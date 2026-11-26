"""Stage 39: exact radial disintegration substitution, Q[a,e] coefficient pencils.
No factorization of the neutral fourth-trace moment or of the flow norm.
The optional builder reads Stage38 mathematical JSON, never executes old code.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib, json

def add(d,k,v):
    if v:
        d[k]=d.get(k,F(0))+v
        if not d[k]: del d[k]

def axpy(d,r,v=F(1)):
    for k,c in r.items():add(d,k,c*v)
    return d

def decode(r):return {(int(k),int(a),int(e)):F(c) for k,a,e,c in r}
def encode(r):return [[k,a,e,str(c)] for (k,a,e),c in sorted(r.items())]
def decvec(r):return {int(i):F(c) for i,c in r}
def encvec(r):return [[i,str(c)] for i,c in sorted(r.items())]

# These coordinates are actual A times other observables, not generic products.
# old -> (old retained coordinate, power of a)
RADIAL={1:(0,1),3:(0,2),4:(2,1),26:(7,1),27:(8,1),28:(9,1),29:(10,1),30:(11,1)}

def substitute(row,kept):
    index={k:j for j,k in enumerate(kept)};out={}
    for k,ep,c in row:
        k=int(k);target,ap=RADIAL.get(k,(k,0))
        add(out,(index[target],ap,int(ep)),F(c))
    return out

def build(source):
    source=Path(source);raw=json.loads(source.read_text());n=raw['dimension']
    if n!=146:raise ValueError('Expected Stage38 dimension146')
    kept=[k for k in range(n) if k not in RADIAL]
    # Verify all named radius products directly in the recorded low-observable basis.
    bn=[raw['basis_names'][i] for i in raw['kept_basis_indices']]
    if bn[:7]!=['m0','m1','d1','m2','dA','d2','JdaggerJ']:raise ValueError('Unexpected observable order')
    for k in (26,27,28,29,30):
        target,_=RADIAL[k];fs=bn[k]['trace_factors'];other=bn[target]['trace_factors']
        if sorted(fs)!=sorted(other+[[0,1]]):raise ValueError('Radius product provenance failed')
    blocks=[]
    for b in raw['blocks']:
        ent=[[substitute(r,kept) for r in row] for row in b['entries']]
        # This is exact full-rank congruence modulo a proven zero kernel, not a
        # small-eigenvalue threshold. Recover G3 = T(a) G2 T(a)^T.
        basis=b['basis'];projection=[{i:F(1)} for i in range(len(ent))]
        if b['name']=='S_G3':
            want=[[{(0,0,0):F(1)},{(0,1,0):F(1)},{(1,0,0):F(1)}],
                  [{(0,1,0):F(1)},{(0,2,0):F(1)},{(1,1,0):F(1)}],
                  [{(1,0,0):F(1)},{(1,1,0):F(1)},{(2,0,0):F(1)}]]
            if ent!=want:raise AssertionError('Incorrect radial Gram kernel')
            ids=[0,2];ent=[[ent[i][j] for j in ids] for i in ids]
            basis=[basis[i] for i in ids];projection=[{i:F(1)} for i in ids]
        for i in range(len(ent)):
            for j in range(len(ent)):
                if ent[i][j]!=ent[j][i]:raise AssertionError('asymmetry')
        blocks.append({'name':'S_G2' if b['name']=='S_G3' else b['name'],
             'basis':basis,'original_dimension':b['dimension'],'dimension':len(ent),
             'radial_projection':[encvec(v) for v in projection],
             'entries':[[encode(r) for r in row] for row in ent]})
    return {'stage':39,'parameters':['a: conditioned radius, not unconditional mean','e: fixed subsequential energy'],
        'source_pencil_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'source_dimension':146,'dimension':len(kept),'normalized_dimension':len(kept)-1,
        'source_kept_coordinates':kept,'coordinate_names':[bn[k] for k in kept],
        'radial_rules':[[k,v[0],v[1]] for k,v in sorted(RADIAL.items())],
        'claim':'Necessary fibers of the all-radially-weighted limit. Not finite-N conditioning; not factorization.',
        'blocks':blocks}

def load(path):
    d=json.loads(Path(path).read_text());b=[dict(x,entries=[[decode(r) for r in row] for row in x['entries']]) for x in d['blocks']]
    return d,b

def congruence(entries,cols):
    out=[]
    for u in cols:
        row=[]
        for v in cols:
            r={}
            for i,c in u.items():
                for j,d in v.items():axpy(r,entries[i][j],c*d)
            row.append(r)
        out.append(row)
    return out

def echelon(rows):
    piv={}
    for rr in rows:
        r=dict(rr)
        while r:
            k=min(r)
            if k not in piv:
                v=r[k];piv[k]={j:c/v for j,c in r.items()};break
            axpy(r,piv[k],-r[k])
    return piv

def kernel(rows,n):
    P=echelon({i:F(c) for i,c in enumerate(r) if c} for r in rows)
    out=[]
    for j in range(n):
        if j in P:continue
        v={j:F(1)}
        for k in sorted(P,reverse=True):
            t=-sum(c*v.get(i,0) for i,c in P[k].items() if i!=k)
            if t:v[k]=t
        out.append(v)
    return out

def psd_pivots(mat):
    a=[list(r) for r in mat];p=[]
    while a:
        if any(a[i][i]<0 for i in range(len(a))):raise AssertionError('negative diagonal')
        k=next((i for i in range(len(a)) if a[i][i]>0),None)
        if k is None:
            if any(v for r in a for v in r):raise AssertionError('nonzero PSD zero diagonal')
            break
        ids=[k]+[i for i in range(len(a)) if i!=k];a=[[a[i][j] for j in ids] for i in ids]
        d=a[0][0];p.append(d)
        a=[[a[i][j]-a[i][0]*a[0][j]/d for j in range(1,len(a))] for i in range(1,len(a))]
    return p

def faces(data,source38):
    # Rebuild the known degree-ten recession direction from observable semantics.
    src=json.loads(Path(source38).read_text());rows=[decvec(r) for r in src['basis_rows']]
    oldray={j:rows[k].get(639,F(0)) for j,k in enumerate(src['kept_basis_indices'])}
    kept=data['source_kept_coordinates'];ray={j:oldray[k] for j,k in enumerate(kept) if oldray[k]}
    if any(oldray[k] for k in RADIAL):raise AssertionError('radial equalities vary along ray')
    result=[]
    for b in data['blocks']:
        E=[[decode(r) for r in row] for row in b['entries']];n=len(E);R=[]
        for rr in E:
            row=[]
            for r in rr:
                v={}
                for (k,a,e),c in r.items():add(v,(a,e),c*ray.get(k,0))
                if any(k!=(0,0) for k in v):raise AssertionError('parameter-dependent recession')
                row.append(v.get((0,0),F(0)))
            R.append(row)
        p=psd_pivots(R);P=kernel(R,n);T=congruence(E,P)
        coeff={}
        for i,row in enumerate(T):
            for j,r in enumerate(row):
                for key,c in r.items():coeff.setdefault((key,i),{})[j]=c
        ec=echelon(coeff[k] for k in sorted(coeff));ids=sorted(ec)
        C=[P[i] for i in ids];K=kernel([[ec[i].get(j,F(0)) for j in range(len(P))] for i in sorted(ec)],len(P))
        for k in K:
            for row in T:
                rem={}
                for j,c in k.items():axpy(rem,row[j],c)
                if rem:raise AssertionError('common kernel does not annihilate coefficient')
        E2=congruence(E,C)
        result.append({'name':b['name'],'dimension':len(C),'source_dimension':n,
           'projection':[encvec(c) for c in C],'ray_pivots':[str(v) for v in p],
           'ray_matrix':[[str(v) for v in r] for r in R],
           'common_nullity':len(K),'common_kernels_in_ray_compression':[encvec(k) for k in K],
           'entries':[[encode(r) for r in row] for row in E2]})
    return {'ray':encvec(ray),'blocks':result,'scope':'Exact common kernels and one recession projection only; no completion/minimal-face claim.'}

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--source',required=True,type=Path);ap.add_argument('--out',required=True,type=Path)
    aa=ap.parse_args();d=build(aa.source);aa.out.mkdir(parents=True,exist_ok=True)
    (aa.out/'radial_pencil.json').write_text(json.dumps(d,sort_keys=True,separators=(',',':'))+'\n')
    f=faces(d,aa.source);(aa.out/'radial_faces.json').write_text(json.dumps(f,sort_keys=True,separators=(',',':'))+'\n')
    print('dimension',d['dimension'],'face dims',[b['dimension'] for b in f['blocks']])
