"""Exact continuation of positive-atom certificates on a radius interval.
A rational preconditioner defines a contraction, proving a unique exact solution
of EVERY stationarity equation, not accepting a small numerical residual.
"""
from fractions import Fraction as F
from math import comb
import json,time
from pathlib import Path
from interval_arithmetic import I,matrix,matmul,mv
from radial_model import add,axpy,load,decode,encode

def bernstein(poly,box):
    lo,hi,el,eh=map(F,box);da=hi-lo;de=eh-el
    na=max((k[0] for k in poly),default=0);ne=max((k[1] for k in poly),default=0);power={}
    for (a,e),c in poly.items():
        for i in range(a+1):
            for j in range(e+1):add(power,(i,j),c*comb(a,i)*lo**(a-i)*da**i*comb(e,j)*el**(e-j)*de**j)
    out=[]
    for i in range(na+1):
        row=[]
        for j in range(ne+1):row.append(sum(v*F(comb(i,k),comb(na,k))*F(comb(j,l),comb(ne,l)) for (k,l),v in power.items() if k<=i and l<=j))
        out.append(row)
    return out

def rebuild_expressions(blocks,atoms):
    out=[]
    for at in atoms:
        b=blocks[int(at['block'])];co={int(k):F(v) for k,v in at['coefficients']}
        if not co:raise ValueError('zero atom')
        norm=sum(v*v for v in co.values());r={}
        for i,c in co.items():
            for j,d in co.items():
                for key,v in b['entries'][i][j].items():add(r,key,v*c*d/norm)
        out.append(r)
    return out

def verify_system(pencil,cert,verbose=True):
    data,blocks=load(pencil);expr=rebuild_expressions(blocks,cert['atoms'])
    if [encode(r) for r in expr]!=cert['expressions']:raise AssertionError('raw atom mismatch')
    a0=F(cert['radius']);keys=[tuple(k) for k in cert['constraint_keys']]
    ids=cert['dependent_indices'];free=cert['free_indices'];fw=list(map(F,cert['free_weights']))
    n=len(ids)
    if sorted(ids+free)!=list(range(len(expr))) or any(v<=0 for v in fw):raise AssertionError('column partition/positive free weights')
    syz=[tuple(k) for k in cert['omitted_row_syzygy']]
    allkeys=sorted({(k,e) for r in expr for k,ap,e in r if k})
    if sorted(keys+[(77,0)])!=allkeys or (77,0) not in syz:raise AssertionError('omitted row provenance')
    for r in expr:
        for ap in (0,1,2):
            if sum(r.get((k,ap,e),F(0)) for k,e in syz):raise AssertionError('polynomial syzygy fails')
    m0=[];m1=[];b0=[];b1=[]
    for k,e in keys:
        if any(kk==k and ee==e and ap>1 for r in expr for kk,ap,ee in r):raise AssertionError('nonlinear nuisance matrix')
        r0=[r.get((k,0,e),F(0))+a0*r.get((k,1,e),F(0)) for r in expr]
        r1=[r.get((k,1,e),F(0)) for r in expr]
        m0.append([r0[i] for i in ids]);m1.append([r1[i] for i in ids]);b0.append(-sum(r0[j]*v for j,v in zip(free,fw)));b1.append(-sum(r1[j]*v for j,v in zip(free,fw)))
    m0.append([F(1)]*n);m1.append([F(0)]*n);b0.append(1-sum(fw));b1.append(F(0))
    for key,M in [('M0',m0),('M1',m1)]:
        if [[F(v) for v in r] for r in cert[key]]!=M:raise AssertionError('coefficient matrix mismatch')
    if list(map(F,cert['b0']))!=b0 or list(map(F,cert['b1']))!=b1:raise AssertionError('rhs mismatch')
    R=matrix(cert['preconditioner']);M0=matrix(m0);M1=matrix(m1)
    v0=list(map(F,cert['center_weights']));v1=list(map(F,cert['linear_weights']))
    if len(R)!=n or any(len(r)!=n for r in R) or min(v0)<=0:raise AssertionError('invalid center/preconditioner')
    v0i=list(map(I,v0));v1i=list(map(I,v1));B0=list(map(I,b0));B1=list(map(I,b1))
    t=time.time();Q0=matmul(R,M0);Q1=matmul(R,M1)
    q0=q1=F(0)
    for i in range(n):
        z0=z1=F(0)
        for j in range(n):
            z0+=((I(int(i==j))-Q0[i][j])*v0[j]/v0[i]).absupper()
            z1+=(Q1[i][j]*v0[j]/v0[i]).absupper()
        q0=max(q0,z0);q1=max(q1,z1)
    ev0=mv(M0,v0i);ev1=mv(M1,v0i);ev2=mv(M0,v1i)
    E0=mv(R,[B0[i]-ev0[i] for i in range(n)])
    E1=mv(R,[B1[i]-ev1[i]-ev2[i] for i in range(n)])
    E2=mv(Q1,v1i)
    d0=max(E0[i].absupper()/v0[i] for i in range(n));d1=max(E1[i].absupper()/v0[i] for i in range(n));d2=max(E2[i].absupper()/v0[i] for i in range(n))
    if verbose:print('preconditioner checked',n,'seconds',time.time()-t,'q0',float(q0),'q1',float(q1),flush=True)
    # Exact nominal objective: c(a,e) dot [v0+(a-a0)v1], plus fixed free weights.
    poly={}
    for j,i in enumerate(ids):
        for (k,ap,ep),c in expr[i].items():
            if k:continue
            add(poly,(ap,ep),c*(v0[j]-a0*v1[j]));add(poly,(ap+1,ep),c*v1[j])
    for i,w in zip(free,fw):
        for (k,ap,ep),c in expr[i].items():
            if not k:add(poly,(ap,ep),w*c)
    return {'expr':expr,'v0':v0,'v1':v1,'ids':ids,'q0':q0,'q1':q1,'d0':d0,'d1':d1,'d2':d2,'poly':poly,'dimension':n,
         'full_nuisance_rows':len(keys)+1,'a0':a0,'elapsed_preconditioner':time.time()-t}

def prove_box(s,box,verbose=True):
    al,ah,el,eh=map(F,box);a0=s['a0'];rad=max(abs(al-a0),abs(ah-a0));q=s['q0']+rad*s['q1']
    if q>=1:raise AssertionError('not a contraction')
    eps=(s['d0']+rad*s['d1']+rad*rad*s['d2'])/(1-q)
    positivity=min(F(1)-rad*abs(d/v)-eps for v,d in zip(s['v0'],s['v1']))
    if positivity<=0:raise AssertionError('weights not certified positive')
    B=bernstein(s['poly'],[al,ah,el,eh]);nom=max(x for row in B for x in row)
    # Weighted error on the objective, separately bounded for each dependent atom.
    mag=F(0)
    for j,i in enumerate(s['ids']):
        poly={(ap,ep):c for (k,ap,ep),c in s['expr'][i].items() if k==0}
        bb=bernstein(poly,[al,ah,el,eh]);mag+=s['v0'][j]*max(abs(x) for row in bb for x in row)
    err=mag*eps;upper=nom+err
    if verbose:print('box',list(map(str,box)),'q',float(q),'relerr',float(eps),'obj nom',float(nom),'error',float(err),'upper',float(upper),flush=True)
    if upper>=0:raise AssertionError('no strict negative uniform objective')
    return {'box_radius_energy':list(map(str,[al,ah,el,eh])),
       'center_radius':str(a0),'system_dimension':s['dimension'],'all_nuisance_rows':s['full_nuisance_rows'],
       'preconditioned_center_residual_norm':str(s['q0']),'parameter_operator_norm':str(s['q1']),
       'contraction_bound':str(q),'solution_relative_error_bound':str(eps),'positive_weight_relative_margin':str(positivity),
       'nominal_bernstein_coefficients':[[str(v) for v in row] for row in B],
       'objective_error_bound':str(err),'strict_negative_margin':str(-upper),
       'interpretation':'Weights are the EXACT solution of a radius-parametric square rational linear system. Preconditioner only verifies existence/invertibility and enclosure. All energies in the box use the same weights at a fixed radius.'}

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--pencil',type=Path,required=True);p.add_argument('--cert',type=Path,required=True);p.add_argument('--box',nargs=4,required=True);p.add_argument('--output',type=Path)
    a=p.parse_args();s=verify_system(a.pencil,json.loads(a.cert.read_text()));r=prove_box(s,a.box)
    if a.output:a.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
