"""Stage39 independent finite checks; standard library only.
Analytical hypotheses and disintegration remain the explicit written proof.
"""
if not __debug__:raise RuntimeError('Do not disable assertions with -O')
from pathlib import Path
import json,hashlib,argparse,copy,time
from fractions import Fraction as F
from radial_model import build,faces,load,decode,encode,decvec,congruence,psd_pivots,axpy,add
from continuation_certificate import verify_system,prove_box
from exact_tests import conditioning_reference,arithmetic_references
from interval_arithmetic import I
ROOT=Path(__file__).resolve().parent
CORE=['model_structure.json','conditioning_reference.json','arithmetic_references.json','lower_uniform_certificate.json',
 'upper_uniform_certificate.json','radial_consequences.json','negative_controls.json','operator_atoms.json','summary.json']

def save(out,name,obj):
    out.mkdir(parents=True,exist_ok=True);(out/name).write_text(json.dumps(obj,sort_keys=True,ensure_ascii=False,indent=2)+'\n')

def structure(d,bs):
    if d['dimension']!=138 or d['source_dimension']!=146:raise AssertionError('wrong conditional dimension')
    if d['radial_rules']!=[[1,0,1],[3,0,2],[4,2,1],[26,7,1],[27,8,1],[28,9,1],[29,10,1],[30,11,1]]:raise AssertionError('invalid conditioning rules')
    if bs[-2]['entries']!=[[{(0,0,0):F(1)},{(1,0,0):F(1)}],[{(1,0,0):F(1)},{(2,0,0):F(1)}]]:raise AssertionError('neutral quartic fluctuation not retained')
    if bs[-1]['entries']!=[[{(3,0,0):F(1)}]]:raise AssertionError('flow norm missing')
    count=0
    for b in bs:
        for i,row in enumerate(b['entries']):
            for j,r in enumerate(row):
                if r!=b['entries'][j][i]:raise AssertionError('matrix not symmetric')
                if any(k<0 or k>=138 or ap<0 or ep<0 for k,ap,ep in r):raise AssertionError('polynomial key')
                count+=1
    # Prove the stored one-ray/common-kernel projection, although final atoms are unprojected.
    fd=json.loads((ROOT/'inputs/radial_faces.json').read_text());ray=decvec(fd['ray']);dim=[]
    for b,f in zip(bs,fd['blocks']):
        E=b['entries'];R=[]
        for row in E:
            rr=[]
            for r in row:
                v={}
                for (k,a,e),c in r.items():add(v,(a,e),c*ray.get(k,F(0)))
                if any(k!=(0,0) for k in v):raise AssertionError('parameter-dependent ray')
                rr.append(v.get((0,0),F(0)))
            R.append(rr)
        if [[str(v) for v in r] for r in R]!=f['ray_matrix']:raise AssertionError('ray coefficients changed')
        if list(map(str,psd_pivots(R)))!=f['ray_pivots']:raise AssertionError('ray positivity failure')
        P=[decvec(v) for v in f['projection']]
        for v in P:
            if any(sum(row[j]*c for j,c in v.items()) for row in R):raise AssertionError('projection outside exposing kernel')
        C=congruence(E,P)
        if [[encode(r) for r in row] for row in C]!=f['entries']:raise AssertionError('projection mismatch')
        dim.append(len(P))
    prod=json.loads((ROOT/'inputs/radial_product_witnesses.json').read_text());pc=0
    for ob in prod['checks']:
        r=decode(ob['product_expression']);s=decode(ob['single_expression'])
        for (k,ap,ep),v in s.items():add(r,(k,ap+1,ep),-v)
        if r:raise AssertionError('A times observable identity')
        pc+=1
    return {'source_dimension':146,'conditional_homogeneous_dimension':138,'normalized_affine_dimension':137,
       'new_independent_radial_relations':8,'constant_elimination_minor':1,'all_source_radius_products_checked':pc,
       'block_count':len(bs),'full_block_dimensions':[len(b['entries']) for b in bs],
       'exact_symmetric_matrix_entries':count,'discovery_projection_dimensions':dim,
       'neutral_G2_and_flow_retained':True,'conditional_kernel':['-a','1','0'],
       'scope':'Exact eight-row conditional substitution. Scalar G3 compression is equivalent; recession projections used only for discovery are necessary, not asserted sufficient.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-pencils',type=Path);ap.add_argument('--output-dir',type=Path,default=ROOT/'reproduced');args=ap.parse_args()
    p=ROOT/'inputs/radial_pencil.json';d,bs=load(p);out=args.output_dir
    if args.source_pencils:
        v=build(args.source_pencils)
        if v!=d:raise AssertionError('source rebuild not equal')
        f=faces(v,args.source_pencils)
        if f!=json.loads((ROOT/'inputs/radial_faces.json').read_text()):raise AssertionError('source face rebuild not equal')
        print('SOURCE REBUILD: exact radial model and projections match.',flush=True)
    st=structure(d,bs);save(out,'model_structure.json',st)
    save(out,'conditioning_reference.json',conditioning_reference());save(out,'arithmetic_references.json',arithmetic_references())
    config=json.loads((ROOT/'inputs/parameters.json').read_text());systems={};certs={};operators={}
    for side in ('lower','upper'):
        c=json.loads((ROOT/f'inputs/{side}_refined_continuation.json').read_text());certs[side]=c
        if list(map(F,c['energy_interval']))!=list(map(F,config['energy_interval'])):raise AssertionError('energy interval narrowed')
        s=verify_system(p,c);systems[side]=s;proof=prove_box(s,config[side+'_box']);save(out,side+'_uniform_certificate.json',proof)
        ops=[]
        for idx,a in enumerate(c['atoms']):
            b=bs[a['block']];r={}
            for j,v in a['coefficients']:
                for w,q in b['basis'][j]:add(r,tuple(w),F(v)*F(q))
            norm=sum(F(v)**2 for j,v in a['coefficients'])
            if idx in c['free_indices']:weight={'fixed':c['free_weights'][c['free_indices'].index(idx)]}
            else:weight={'exact_solution_row':c['dependent_indices'].index(idx)}
            ops.append({'index':idx,'block':b['name'],'squared_coefficient_norm':str(norm),'normalized_weight':weight,
             'operator':[[list(w),str(v)] for w,v in sorted(r.items())],
             'interpretation':'scalar trace polynomial/norm' if b['name'].startswith('S_') else 'matrix polynomial norm/ground-state form'})
        operators[side]=ops
    save(out,'operator_atoms.json',operators)
    lo,hi=F(config['new_lower']),F(config['new_upper']);oldlo,oldhi=map(F,config['inherited_radius_interval'])
    el,eh=map(F,config['energy_interval'])
    if list(map(F,config['lower_box']))!=[oldlo,lo,el,eh] or list(map(F,config['upper_box']))!=[hi,oldhi,el,eh]:raise AssertionError('incomplete exclusion coverage')
    if not oldlo<lo<hi<oldhi:raise AssertionError('not a strict refinement')
    var=(hi-lo)**2/4
    rd={'limiting_support':[str(lo),str(hi)],'width':str(hi-lo),'limsup_variance_upper':str(var),
      'old_support':[str(oldlo),str(oldhi)],'relative_variance_reduction':str(1-var/((oldhi-oldlo)**2/4)),
      'energy_interval':list(map(str,[el,eh])),
      'quantifier':'Conditioned fields exist on a common mu_A-full set. Each certified closed strip is excluded almost everywhere; together with inherited support this yields the stated closed necessary support interval.',
      'finite_N_hard_support':False,'global_factorization_assumed':False,'sharp_endpoint_claim':False}
    save(out,'radial_consequences.json',rd)
    rejected=[]
    def fail(name,fn):
        try:fn()
        except (AssertionError,ValueError,ZeroDivisionError,KeyError,IndexError):rejected.append(name);return
        raise AssertionError('bad input was accepted: '+name)
    bad=copy.deepcopy(bs);bad[-1]['entries']=[[{}]];fail('remove unresolved flow norm',lambda:structure(d,bad))
    bad=copy.deepcopy(bs);bad[-2]['entries'][1][1]={(1,0,0):F(1)};fail('erase neutral fourth-trace second moment',lambda:structure(d,bad))
    bad=copy.deepcopy(d);bad['radial_rules'][0][2]=0;fail('incorrect conditional radius identity',lambda:structure(bad,bs))
    fail('interval zero divisor',lambda:I(-1,1).inv())
    bad=copy.deepcopy(certs['lower']);bad['free_weights'][0]='-1';fail('negative free atom weight',lambda:verify_system(p,bad,False))
    bad=copy.deepcopy(certs['lower']);bad['atoms'][0]['coefficients'][0][1]=str(F(bad['atoms'][0]['coefficients'][0][1])+1);fail('tampered raw operator',lambda:verify_system(p,bad,False))
    bad=copy.deepcopy(certs['lower']);bad['omitted_row_syzygy']=bad['omitted_row_syzygy'][:-1];fail('false stationarity redundancy',lambda:verify_system(p,bad,False))
    bad=copy.deepcopy(systems['lower']);bad['q0']=F(1);fail('preconditioner without contraction',lambda:prove_box(bad,config['lower_box'],False))
    fail('unjustified lower extension .38905',lambda:prove_box(systems['lower'],[oldlo,F('.38905'),el,eh],False))
    fail('unjustified upper extension .3890',lambda:prove_box(systems['upper'],[F('.3890'),oldhi,el,eh],False))
    bad=copy.deepcopy(systems['upper']);bad['v1'][0]=F(10**9);fail('weight loses positivity on strip',lambda:prove_box(bad,config['upper_box'],False))
    save(out,'negative_controls.json',{'count':len(rejected),'rejected':rejected})
    summary={
        'exact_weight_definition': 'One positive free normalized weight and 138 exact dependent weights obtained by a square Q[a] linear system. Dyadic interval contraction verifies positive solutions on whole radius strips. Weights do not depend on e.',
        'model_structure': st,
        'negative_controls': len(rejected),
        'result': rd,
        'uniform_certificate_atom_counts': [len(certs[s]['atoms']) for s in ('lower','upper')],
    }
    save(out,'summary.json',summary)
    print('CERTIFIED LIMITING SUPPORT:',lo,'<= A <=',hi,flush=True)
    print('CERTIFIED limsup Var(A) <=',var,flush=True)
    print('All Stage-39 finite exact checks passed.',flush=True)
if __name__=='__main__':main()
