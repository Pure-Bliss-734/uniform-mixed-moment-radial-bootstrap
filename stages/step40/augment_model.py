"""Build only the two nonzero color-centering additions to the Stage39 pencil.
Four additional candidates are identically old principal blocks and omitted.
Kernels are solved over Q, never classified by floating eigenvalues.
"""
from pathlib import Path
from fractions import Fraction as F
import json,copy,hashlib
from trace_images import rebuild,add,encoded

def decode(r):return {(k,a,e):F(c) for k,a,e,c in r}
def axpy(r,s,c=F(1)):
    for k,v in s.items():add(r,k,c*v)
    return r

def row_reduce(rows):
    piv={}
    for rr in rows:
        r=dict(rr)
        while r:
            k=min(r)
            if k not in piv:
                v=r[k];piv[k]={j:c/v for j,c in r.items()};break
            axpy(r,piv[k],-r[k])
    return piv

def projection(E):
    n=len(E);equations={}
    for i,row in enumerate(E):
        for j,r in enumerate(row):
            for key,c in r.items():equations.setdefault((key,i),{})[j]=c
    P=row_reduce(equations[k] for k in sorted(equations));ids=sorted(P);null=[]
    for j in range(n):
        if j in P:continue
        z={j:F(1)}
        for i in sorted(P,reverse=True):
            v=-sum(c*z.get(k,F(0)) for k,c in P[i].items() if k!=i)
            if v:z[i]=v
        for row in E:
            rem={}
            for k,c in z.items():axpy(rem,row[k],c)
            if rem:raise AssertionError('false exact common kernel')
        null.append(z)
    return ids,null

def construct(baseline):
    baseline=Path(baseline);p=baseline/'inputs/radial_pencil.json';raw=json.loads(p.read_text());face=json.loads((baseline/'inputs/radial_faces.json').read_text());im=rebuild(p)
    aug=copy.deepcopy(raw);fac=copy.deepcopy(face);added=[];redundant=[]
    for b in im['blocks']:
        if not any(r for row in b['subtracted_scalar_gram'] for r in row):
            old=raw['blocks'][b['source_block']];idx=b['source_indices']
            if b['entries']!=[[old['entries'][i][j] for j in idx] for i in idx]:raise AssertionError('false redundancy')
            redundant.append(b['name']);continue
        E=[[decode(r) for r in row] for row in b['entries']];ids,null=projection(E);n=len(E)
        aug['blocks'].append({'name':b['name'],'basis':b['basis'],'dimension':n,'entries':b['entries'],'color_centered':True})
        fac['blocks'].append({'name':b['name'],'dimension':len(ids),'source_dimension':n,'projection':[[[i,'1']] for i in ids],
          'entries':[[b['entries'][i][j] for j in ids] for i in ids],'common_nullity':len(null),'color_centered':True})
        added.append({'name':b['name'],'raw_dimension':n,'discovery_dimension':len(ids),
          'exact_common_kernel_vectors':[[[i,str(c)] for i,c in sorted(v.items())] for v in null],
          'source_indices':b['source_indices']})
    aug['stage']=40;aug['augmentation']='Finite-N color-centered ordinary PSD; two nonzero subtracted trace Grams.'
    info={'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'new_variables':0,'new_linear_equalities':0,
      'base_test_weight':4,'base_norm_weight':8,'added_blocks':added,'redundant_principal_blocks':redundant,
      'complete_block_count':len(aug['blocks']),'complete_matrix_entries':sum(len(b['entries'])**2 for b in aug['blocks']),
      'explicitly_new_not_logical_independence_claim':'Two additional PSD blocks were not in the prior list. Logical nonredundancy against the ENTIRE old SDP is not proved; exact principal redundancy is proved for the other four.',
      'neutral_variance_and_flow_remain':True}
    return aug,fac,info,im

def write(baseline,output):
    d=Path(output);d.mkdir(parents=True,exist_ok=True);a,f,r,im=construct(baseline)
    (d/'augmented_pencil.json').write_text(json.dumps(a,separators=(',',':'))+'\n');(d/'augmented_faces.json').write_text(json.dumps(f,separators=(',',':'))+'\n')
    return r,im

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--baseline',required=True,type=Path);p.add_argument('--output-dir',required=True,type=Path);a=p.parse_args();r,im=write(a.baseline,a.output_dir);print(json.dumps(r,indent=2))
