"""Step40 complete new finite audit and centered-color certificates.
The radial baseline is reconstructed from the local Step 39 inputs.
"""
if not __debug__:raise RuntimeError('Proof checks must not run with -O')
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,tempfile,copy,sys
from collections import Counter
from verify_baseline import reproduce
from augment_model import construct
from independent_interval_audit import audit,reference_tests,Ball,dot,bernstein
from coupling_checks import extract,scalar_references
from finite_operator_checks import run as finite_checks
from trace_images import word_image
from support_audit import audit as support_audit
ROOT=Path(__file__).resolve().parent
CORE=['baseline_reproduction.json','model_augmentation.json','trace_images.json','coupled_constraints.json',
      'independent_interval_references.json','formal_references.json','finite_N_contacts.json',
      'lower_independent_audit.json','upper_independent_audit.json','lower_new_uniform.json','upper_new_uniform.json',
      'radial_consequences.json','selected_support_audit.json','negative_controls.json','summary.json']

def write(out,name,data):
    out.mkdir(parents=True,exist_ok=True);(out/name).write_text(json.dumps(data,sort_keys=True,ensure_ascii=False,indent=2)+'\n')

def prove(baseline,out,baseline_info):
    cfg=json.loads((ROOT/'inputs/parameters.json').read_text())
    write(out,'baseline_reproduction.json',baseline_info)
    aug,fac,structure,images=construct(baseline)
    if aug!=json.loads((ROOT/'inputs/augmented_pencil.json').read_text()):raise AssertionError('new full pencil not correctly rebuilt')
    if fac!=json.loads((ROOT/'inputs/augmented_faces.json').read_text()):raise AssertionError('new exact discovery kernels do not rebuild')
    write(out,'model_augmentation.json',structure);write(out,'trace_images.json',images)
    coupling=extract(images);write(out,'coupled_constraints.json',coupling)
    write(out,'independent_interval_references.json',reference_tests());write(out,'formal_references.json',scalar_references())
    write(out,'finite_N_contacts.json',finite_checks())
    oldcfg=json.loads((baseline/'inputs/parameters.json').read_text());oldproof={};newproof={};atomcounts={}
    for side in ('lower','upper'):
        print('Independent baseline interval audit:',side,flush=True)
        r=audit(baseline/'inputs/radial_pencil.json',baseline/f'inputs/{side}_refined_continuation.json',oldcfg[side+'_box'])
        oldproof[side]=r;write(out,side+'_independent_audit.json',r)
        print('New centered-color certificate:',side,flush=True)
        c=json.loads((ROOT/f'inputs/{side}_continuation.json').read_text())
        if list(map(F,c['energy_interval']))!=list(map(F,cfg['energy_interval'])):raise AssertionError('uncovered energy interval')
        r=audit(ROOT/'inputs/augmented_pencil.json',ROOT/f'inputs/{side}_continuation.json',cfg[side+'_box'])
        newproof[side]=r;write(out,side+'_new_uniform.json',r)
        counts=Counter()
        for at in c['atoms']:
            name=aug['blocks'][at['block']]['name']
            counts['ordinary' if name.startswith('M_') else 'ground' if name.startswith('N_') else 'centered_color' if name.startswith('C_') else 'scalar_flow_or_covariance']+=1
        atomcounts[side]=dict(sorted(counts.items()))
        print('q<',float(F(r['contraction_upper'])),'objective<',float(F(r['objective_upper'])),flush=True)
    oldl,oldh=map(F,cfg['inherited_radius_interval']);l,h=map(F,cfg['new_radius_interval']);el,eh=map(F,cfg['energy_interval'])
    if not oldl<l<h<oldh:raise AssertionError('not a strict refinement')
    if list(map(F,cfg['lower_box']))!=[oldl,l,el,eh] or list(map(F,cfg['upper_box']))!=[h,oldh,el,eh]:raise AssertionError('parameter coverage gap')
    var=(h-l)**2/4;oldvar=(oldh-oldl)**2/4
    radial={'necessary_limiting_support':list(map(str,[l,h])),'same_energy_interval':list(map(str,[el,eh])),
      'limsup_variance_upper':str(var),'variance_upper_display':format(float(var),'.12g'),
      'old_support':list(map(str,[oldl,oldh])),'relative_variance_bound_reduction':str(1-var/oldvar),
      'full_continuous_rectangles':True,'factorization_assumed':False,'finite_N_hard_support':False,
      'no_primal_feasibility_claim':True,'same_energy_two_radius_primal_test':'not required as fallback because new certified exclusions were obtained; no complete primal pair was constructed'}
    write(out,'radial_consequences.json',radial)
    support = {side: support_audit(ROOT/f'inputs/{side}_continuation.json') for side in ('lower','upper')}
    if not all(x['full_column_rank'] for x in support.values()):
        raise AssertionError('Chosen support rank check failed')
    write(out,'selected_support_audit.json',support)
    rejected=[]
    def fail(label,fn):
        try:fn()
        except (AssertionError,ValueError,ZeroDivisionError,KeyError,IndexError):rejected.append(label);return
        raise AssertionError('false input accepted: '+label)
    bad=copy.deepcopy(images)
    for b in bad['blocks']:
        if b['name']=='C_q1':
            i=next(i for i,p in enumerate(b['basis']) if p==[[[0,0,3],'1']]);b['entries'][i][i]=[r for r in b['entries'][i][i] if r[:3]!=[0,1,0]]
    fail('delete A/2 contact term',lambda:extract(bad))
    bad=copy.deepcopy(images)
    for b in bad['blocks']:
        if b['name']=='C_q0_even':
            b['entries']=[[[r for r in entry if r[0]!=2] for entry in row] for row in b['entries']]
    fail('delete neutral second moment v',lambda:extract(bad))
    bad=copy.deepcopy(images)
    for b in bad['blocks']:
        if b['name']=='C_q1':b['entries']=[[[r for r in entry if r[0]!=3] for entry in row] for row in b['entries']]
    fail('delete current norm j',lambda:extract(bad))
    fail('unproved weight-six trace image',lambda:word_image((0,0,0,1,1,1)))
    fail('negative interval radius',lambda:Ball(0,-1))
    fail('shape mismatch in fused dot',lambda:dot([Ball.scalar(1)],[]))
    # Deliberate actual-certificate mutations rejected before costly interval work.
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)/'bad.json';c=json.loads((ROOT/'inputs/lower_continuation.json').read_text())
        b=copy.deepcopy(c);b['free_weights'][0]='-1';p.write_text(json.dumps(b));fail('negative free positive atom weight',lambda:audit(ROOT/'inputs/augmented_pencil.json',p,cfg['lower_box']))
        b=copy.deepcopy(c);b['M0'][0][0]=str(F(b['M0'][0][0])+1);p.write_text(json.dumps(b));fail('changed nuisance linear system',lambda:audit(ROOT/'inputs/augmented_pencil.json',p,cfg['lower_box']))
        b=copy.deepcopy(c);b['omitted_row_syzygy']=b['omitted_row_syzygy'][:-1];p.write_text(json.dumps(b));fail('false omitted-row identity',lambda:audit(ROOT/'inputs/augmented_pencil.json',p,cfg['lower_box']))
        b=copy.deepcopy(c);b['atoms'][0]['coefficients'][0][1]=str(F(b['atoms'][0]['coefficients'][0][1])+1);p.write_text(json.dumps(b));fail('changed unprojected operator',lambda:audit(ROOT/'inputs/augmented_pencil.json',p,cfg['lower_box']))
    write(out,'negative_controls.json',{'count':len(rejected),'rejected':rejected})
    summary={
        'atom_counts': atomcounts,
        'baseline_finite_reproduction': baseline_info,
        'consequences': radial,
        'negative_controls': len(rejected),
        'new_coupling': 'j <= a/2+p*u+sqrt((R-u^2-sigma)*(T-p^2)); both variances stay nonnegative, not set to zero',
        'new_model': structure,
        'selected_support_audit': support,
    }
    write(out,'summary.json',summary)
    print('CERTIFIED limiting support:',l,h,flush=True);print('CERTIFIED limsup variance <=',var,flush=True);print('All Step40 required finite checks passed.',flush=True)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    baseline = ROOT.parent / 'step39'
    with tempfile.TemporaryDirectory() as tmp:
        info = reproduce(baseline, Path(tmp))
        prove(baseline, args.output_dir, info)
