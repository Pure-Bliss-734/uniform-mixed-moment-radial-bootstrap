"""Exact common-kernel and one-ray compression of Stage38 pencils.
Compressed positivity is a necessary condition. No finite completion is claimed.
"""
import json
from fractions import Fraction as F
from energy_model import add,axpy,enc,ps_enc,ps_dec,unpack

def echelon(rows):
    piv={}
    for row in rows:
        r=dict(row)
        while r:
            i=min(r)
            if i not in piv:
                v=r[i];piv[i]={k:c/v for k,c in r.items()};break
            axpy(r,piv[i],-r[i])
    return piv

def kernel(M):
    d=len(M);E=echelon({j:c for j,c in enumerate(row) if c} for row in M)
    free=[j for j in range(d) if j not in E];out=[]
    for k in free:
        v={k:F(1)}
        for i in sorted(E,reverse=True):
            a=-sum(c*v.get(j,F(0)) for j,c in E[i].items() if j!=i)
            if a:v[i]=a
        out.append(v)
    return out

def psd_pivots(M):
    M=[list(row) for row in M];piv=[]
    while M:
        if any(M[i][i]<0 for i in range(len(M))):raise AssertionError('negative ray diagonal')
        idx=next((i for i in range(len(M)) if M[i][i]>0),None)
        if idx is None:
            if any(x for row in M for x in row):raise AssertionError('nonzero zero-diagonal ray')
            break
        order=[idx]+[i for i in range(len(M)) if i!=idx];M=[[M[i][j] for j in order] for i in order]
        a=M[0][0];piv.append(a)
        M=[[M[i][j]-M[i][0]*M[0][j]/a for j in range(1,len(M))] for i in range(1,len(M))]
    return piv

def congruence(E,P):
    out=[]
    for u in P:
        row=[]
        for v in P:
            r={}
            for i,c in u.items():
                for j,d in v.items():axpy(r,E[i][j],c*d)
            row.append(r)
        out.append(row)
    return out

def combine_column(P,q):
    r={}
    for j,c in q.items():axpy(r,P[j],c)
    return r

def get_faces(data,progress=False):
    m=unpack(data);basis=[dict((k,F(v)) for k,v in row) for row in data['basis_rows']]
    kept=data['kept_basis_indices'];ray={j:basis[old].get(639,F(0)) for j,old in enumerate(kept)}
    ray={j:c for j,c in ray.items() if c}
    if any(j<7 for j in ray):raise AssertionError('ray changes low observables')
    # Independently check lift equals the original pure-coordinate ray.
    for k,v in m['total_map'].items():
        f={}
        for (j,p),c in v.items():add(f,p,c*ray.get(j,F(0)))
        if f!=({0:F(1)} if k==639 else {}):raise AssertionError('ray lift')
    faces=[]
    for b in m['blocks']:
        E=b['entries'];d=len(E);R=[]
        for row in E:
            vals=[]
            for r in row:
                v={}
                for (j,p),c in r.items():add(v,p,c*ray.get(j,F(0)))
                if any(p!=0 for p in v):raise AssertionError('energy-dependent ray')
                vals.append(v.get(0,F(0)))
            R.append(vals)
        piv=psd_pivots(R);P=kernel(R);T=congruence(E,P)
        # Union of all rows of every coefficient matrix, including energy powers.
        mats={}
        for i,row in enumerate(T):
            for j,r in enumerate(row):
                for (k,p),v in r.items():mats.setdefault((k,p,i),{})[j]=v
        ec=echelon(mats[key] for key in sorted(mats));ids=sorted(ec)
        Q=[P[i] for i in ids];C=congruence(E,Q)
        # Verify every common null vector actually annihilates every coefficient.
        matrows=[[ec[i].get(j,F(0)) for j in range(len(P))] for i in ec]
        # Generic rectangular nullspace via padding does not alter the row span.
        matrows+= [[F(0)]*len(P) for _ in range(len(P)-len(matrows))]
        common=kernel(matrows)
        for v in common:
            for row in T:
                z={}
                for j,a in v.items():axpy(z,row[j],a)
                if z:raise AssertionError('common kernel mismatch')
        item={'name':b['name'],'original_dimension':d,'dimension':len(Q),
            'projection':[enc(v) for v in Q],
            'ray_matrix':[[str(v) for v in row] for row in R],
            'ray_psd_pivots':[str(v) for v in piv],
            'ray_kernel_dimension':len(P),'common_kernel_dimension':len(common),
            'entries':[[ps_enc(r) for r in row] for row in C]}
        faces.append(item)
        if progress:print(b['name'],d,'->',len(Q),'rayrank',len(piv),'common',len(common),flush=True)
    return {'ray':enc(ray),'ray_origin':'old Stage37 free column639 (pure coordinate degree10)',
        'claim':'Exact exposing ray and universal common kernels. Projection is necessary; completion and minimal facial reduction are not claimed.',
        'blocks':faces}

if __name__=='__main__':
    import argparse,time
    from pathlib import Path
    ap=argparse.ArgumentParser();ap.add_argument('--pencil',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();t=time.time();v=get_faces(json.loads(args.pencil.read_text()),True)
    args.output.write_text(json.dumps(v,sort_keys=True,indent=2)+'\n');print('face seconds',time.time()-t)
