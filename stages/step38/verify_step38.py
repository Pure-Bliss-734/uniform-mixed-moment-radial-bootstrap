"""Stage38 verification: standard library only; optional pinned predecessor rebuild."""
if not __debug__:
    raise RuntimeError('Do not use -O: exact verification requires assertions enabled.')
from pathlib import Path
import json,argparse,hashlib,copy,sys
from fractions import Fraction as F
from face_reduction import get_faces
from energy_model import construct,pack,ps_enc,ps_dec,unpack
from certificate_checks import (structure,root_enclosure,verify_identity,upward_short,
    tight_simple_cut,prove_cut,covariance_reference,reject_floating_candidate,exact_references)
ROOT=Path(__file__).resolve().parent
CORE=['model_structure.json','energy_interval.json','lower_certificate_checks.json','upper_certificate_checks.json',
      'radial_consequences.json','covariance_reference.json','exact_references.json','negative_controls.json','face_reduction.json','operator_atoms.json','summary.json']

def save(out,name,obj):
    out.mkdir(parents=True,exist_ok=True)
    (out/name).write_text(json.dumps(obj,sort_keys=True,ensure_ascii=False,indent=2)+'\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-dir',type=Path,help='Directory of reconstructed Step 37 results.')
    ap.add_argument('--output-dir',type=Path,default=ROOT/'reproduced');args=ap.parse_args()
    data=json.loads((ROOT/'inputs/energy_pencil.json').read_text());config=json.loads((ROOT/'inputs/parameters.json').read_text())
    if args.source_dir:
        new=pack(construct(args.source_dir,True))
        if new!=data:raise AssertionError('source rebuild differs from frozen new model')
        print('PINNED SOURCE REBUILD: all transformed matrices/maps equal frozen Stage38 data.',flush=True)
    st=structure(data);save(args.output_dir,'model_structure.json',st)
    faces=get_faces(data);save(args.output_dir,'face_reduction.json',faces)
    print('Exact ray/common-kernel projections recomputed (no minimal-face claim).',flush=True)
    en=root_enclosure();save(args.output_dir,'energy_interval.json',en)
    cuts={};certs={};signs={};operators={}
    for side in ('lower','upper'):
        cert=json.loads((ROOT/f'inputs/{side}_certificate.json').read_text());certs[side]=cert
        poly,proof=verify_identity(data,cert)
        if list(map(F,cert['energy_interval']))!=list(map(F,en['certified_energy_interval'])):raise AssertionError('incorrect energy scope')
        short,scale=upward_short(poly,config['rounding_denominator'])
        cut=F(config[side+'_cut']);sign=prove_cut(short,side,cut,en['certified_energy_interval'])
        if ps_enc(short)!=config[side+'_short_polynomial']:raise AssertionError('short polynomial mismatch')
        cuts[side]=cut;signs[side]=sign
        ops=[]
        for atom in cert['atoms']:
            b=data['blocks'][atom['block']];o={}
            for i,c in atom['coefficients']:
                for w,v in b['basis'][i]:
                    key=tuple(w);o[key]=o.get(key,F(0))+F(c)*F(v)
            o={w:c for w,c in o.items() if c}
            max_weight=max((sum(1 if a<2 else 2 for a in w) for w in o),default=0)
            if max_weight>5:raise AssertionError('atom exceeds declared test degree')
            ops.append({'source_block':b['name'],'weight':atom['weight'],
                'interpretation':'scalar trace polynomial' if b['name'].startswith('S_') else 'matrix word polynomial',
                'operator':[[list(w),str(c)] for w,c in sorted(o.items())],'maximum_weight':max_weight})
        operators[side]=ops
        proof.update({'short_integer_polynomial':ps_enc(short),'normalization_scale':str(scale),
            'rounding_denominator':config['rounding_denominator'],'upward_rounding_all_coefficients_checked':True,'uniform_sign_box':sign})
        save(args.output_dir,side+'_certificate_checks.json',proof)
        print(side.upper()+': positive atoms',proof['atoms'],'all nuisance coefficients exactly zero; whole parameter box excluded.',flush=True)
    save(args.output_dir,'operator_atoms.json',operators)
    lo,hi=cuts['lower'],cuts['upper'];var=(hi-lo)**2/4
    if not lo<hi:raise AssertionError('empty result interval')
    # Explicit nesting relative to old integer quadratics, without importing rounded old roots.
    if not 2275001*lo**2+10**9*lo-388844062>0:raise AssertionError('lower bound not improved')
    if not 115063122*hi**2-10**9*hi+372222252>0:raise AssertionError('upper bound not improved')
    if not var<F(293,2000000000):raise AssertionError('summary variance rounding')
    rd={'limiting_radius_support':[str(lo),str(hi)],'width':str(hi-lo),'limsup_variance_upper':str(var),
        'simple_strict_variance_upper':'293/2000000000',
        'no_finite_N_hard_support':True,'no_factorization_or_sharp_endpoint_claim':True,
        'energy_quantifier':'Both exclusion polynomials work for EVERY e in the certified interval, with the SAME positive weights. e remains fixed on each subsequence.',
        'previous_support_strictly_contains_new_interval':True}
    save(args.output_dir,'radial_consequences.json',rd)
    save(args.output_dir,'covariance_reference.json',covariance_reference())
    refs=exact_references();save(args.output_dir,'exact_references.json',refs)
    rejected=[]
    def fail(name,fn):
        try:fn()
        except (AssertionError,ValueError,KeyError):rejected.append(name);return
        raise AssertionError('negative control not rejected: '+name)
    bad=copy.deepcopy(certs['lower']);bad['atoms'][0]['weight']=str(-F(bad['atoms'][0]['weight']))
    fail('negative atom weight',lambda:verify_identity(data,bad))
    bad2=copy.deepcopy(certs['upper']);bad2['radial_polynomial'][0][2]=str(F(bad2['radial_polynomial'][0][2])+1)
    fail('tampered radial coefficient',lambda:verify_identity(data,bad2))
    for side in ('lower','upper'):
        q=ps_dec(config[side+'_short_polynomial'])
        cut=F(39,100) if side=='lower' else F(388,1000)
        fail('unjustified stronger '+side+' domain',lambda q=q,side=side,cut=cut:prove_cut(q,side,cut,en['certified_energy_interval']))
    bad3=copy.deepcopy(data);bad3['blocks'][-1]['entries']=[[[]]]
    fail('delete J norm',lambda:structure(bad3))
    bad4=copy.deepcopy(data);bad4['blocks'][-2]['entries'][2][2]=[[2,0,'1']]
    fail('replace unknown a4 second moment',lambda:structure(bad4))
    bad5=copy.deepcopy(data);bad5['independent_low_rows'][0][0][2]=str(F(bad5['independent_low_rows'][0][0][2])+1)
    fail('wrong energy closure normalization',lambda:structure(bad5))
    floated=json.loads((ROOT/'inputs/infeasible_support.json').read_text())
    badrank=reject_floating_candidate(data,floated);rejected.append('actual spurious floating negative candidate: exact full-column rank')
    neg={'rejected':rejected,'count':len(rejected),'actual_numeric_failure_certificate':badrank}
    save(args.output_dir,'negative_controls.json',neg)
    su={
        'actual_new_positive_certificates': [len(certs['lower']['atoms']),len(certs['upper']['atoms'])],
        'after': 146,
        'constant_rank_minor': 81,
        'homogeneous_coordinates_before': 150,
        'negative_control_count': len(rejected),
        'radial_consequences': rd,
        'rank_new_closure': 4,
        'scalar_G3_and_J_norm_retained': True,
        'source_closure_inserted_into_all_14_ordinary_ground_blocks': True,
    }
    save(args.output_dir,'summary.json',su)
    print('CERTIFIED LIMITING SUPPORT:',lo,'<= A <=',hi,flush=True)
    print('CERTIFIED limsup Var(A) <=',var,flush=True)
    print('All Stage-38 exact checks passed.',flush=True)
if __name__=='__main__':main()
