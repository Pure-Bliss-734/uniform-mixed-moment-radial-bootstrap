"""Independent decimal center-radius verification of Stage39 continuations.
No import of the previous interval, continuation, Bernstein or matrix code.
The input preconditioner is NOT assumed to be an exact inverse.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb
from random import Random
import json,hashlib

S=10**110

def ceildiv(n,d):return -((-n)//d)
class Ball:
    __slots__=('c','r')
    def __init__(self,c,r):
        self.c=int(c);self.r=int(r)
        if self.r<0:raise ValueError('negative radius')
    @classmethod
    def scalar(cls,q):
        q=F(q);c=(q.numerator*S)//q.denominator
        return cls(c,int(c*q.denominator!=q.numerator*S))
    def __add__(self,b):return Ball(self.c+b.c,self.r+b.r)
    def __neg__(self):return Ball(-self.c,self.r)
    def __sub__(self,b):return self+-b
    def ab(self):return F(abs(self.c)+self.r,S)
    def bounds(self):return F(self.c-self.r,S),F(self.c+self.r,S)

def dot(x,y):
    if len(x)!=len(y):raise ValueError('dot shape')
    mid=0;err=0
    for a,b in zip(x,y):
        mid+=a.c*b.c
        err+=abs(a.c)*b.r+a.r*abs(b.c)+a.r*b.r
    c=mid//S
    return Ball(c,ceildiv(err,S)+int(c*S!=mid))

def matmul(x,y):
    if not y or any(len(a)!=len(y) for a in x):raise ValueError('matrix shape')
    cols=list(zip(*y));return [[dot(a,b) for b in cols] for a in x]
def mv(m,v):return [dot(r,v) for r in m]
def matrix(x):return [[Ball.scalar(v) for v in r] for r in x]
def add(d,k,v):
    if v:
        d[k]=d.get(k,F(0))+v
        if not d[k]:del d[k]
def decode(r):return {(k,a,e):F(c) for k,a,e,c in r}
def encoded(r):return [[k,a,e,str(c)] for (k,a,e),c in sorted(r.items())]

def mono_bernstein(p,n,i,l,h):
    if p>n:raise ValueError('degree')
    return sum(F(comb(i,j)*comb(n-i,p-j),comb(n,p))*h**j*l**(p-j)
               for j in range(max(0,p-(n-i)),min(p,i)+1))
def bernstein(poly,box):
    l,h,el,eh=map(F,box);n=max((p for p,q in poly),default=0);m=max((q for p,q in poly),default=0)
    return [[sum(c*mono_bernstein(p,n,i,l,h)*mono_bernstein(q,m,j,el,eh)
                 for (p,q),c in poly.items()) for j in range(m+1)] for i in range(n+1)]

def audit(pencil,cert,box):
    raw=json.loads(Path(pencil).read_text());c=json.loads(Path(cert).read_text())
    B=[[[decode(r) for r in row] for row in b['entries']] for b in raw['blocks']]
    rows=[]
    for at in c['atoms']:
        v={int(i):F(q) for i,q in at['coefficients']};norm=sum(q*q for q in v.values())
        if norm<=0:raise AssertionError('zero atom')
        rr={};E=B[at['block']]
        for i,x in v.items():
            for j,y in v.items():
                for k,z in E[i][j].items():add(rr,k,x*y*z/norm)
        rows.append(rr)
    if [encoded(r) for r in rows]!=c['expressions']:raise AssertionError('unprojected atoms changed')
    keys=[tuple(k) for k in c['constraint_keys']];syz=[tuple(k) for k in c['omitted_row_syzygy']]
    keyall=sorted({(k,e) for r in rows for k,a,e in r if k})
    missing=set(keyall)-set(keys)
    if len(missing)!=1 or not missing.issubset(syz) or len(keyall)!=len(keys)+1:raise AssertionError('row omission')
    for r in rows:
        apowers={a for k,a,e in r}
        for ap in apowers:
            if sum(r.get((k,ap,e),F(0)) for k,e in syz):raise AssertionError('syzygy not exact')
    ids=c['dependent_indices'];free=c['free_indices'];fw=list(map(F,c['free_weights']));n=len(ids)
    if sorted(ids+free)!=list(range(len(rows))) or min(fw)<=0:raise AssertionError('weights partition')
    a0=F(c['radius']);M0=[];M1=[];b0=[];b1=[]
    for k,e in keys:
        if any(a>1 and kk==k and ee==e for r in rows for kk,a,ee in r):raise AssertionError('nonaffine nuisance')
        x=[r.get((k,0,e),F(0))+a0*r.get((k,1,e),F(0)) for r in rows]
        y=[r.get((k,1,e),F(0)) for r in rows]
        M0.append([x[j] for j in ids]);M1.append([y[j] for j in ids]);b0.append(-sum(x[j]*q for j,q in zip(free,fw)));b1.append(-sum(y[j]*q for j,q in zip(free,fw)))
    M0.append([F(1)]*n);M1.append([F(0)]*n);b0.append(1-sum(fw));b1.append(F(0))
    if len(M0)!=n:raise AssertionError('not square')
    if M0!=[[F(x) for x in row] for row in c['M0']] or M1!=[[F(x) for x in row] for row in c['M1']]:raise AssertionError('matrix coefficients')
    if b0!=list(map(F,c['b0'])) or b1!=list(map(F,c['b1'])):raise AssertionError('rhs coefficients')
    v0=list(map(F,c['center_weights']));v1=list(map(F,c['linear_weights']))
    if len(v0)!=n or min(v0)<=0:raise AssertionError('center positivity')
    R=matrix(c['preconditioner']);X=matrix(M0);Y=matrix(M1);V=list(map(Ball.scalar,v0));D=list(map(Ball.scalar,v1))
    T=matmul(R,X);U=matmul(R,Y);one=Ball.scalar(1);zero=Ball.scalar(0)
    q0=max(sum(((one if i==j else zero)-T[i][j]).ab()*v0[j]/v0[i] for j in range(n)) for i in range(n))
    q1=max(sum(U[i][j].ab()*v0[j]/v0[i] for j in range(n)) for i in range(n))
    xv=mv(X,V);yv=mv(Y,V);xd=mv(X,D)
    E0=mv(R,[Ball.scalar(b0[i])-xv[i] for i in range(n)])
    E1=mv(R,[Ball.scalar(b1[i])-yv[i]-xd[i] for i in range(n)])
    E2=mv(U,D)
    d0=max(z.ab()/v0[i] for i,z in enumerate(E0));d1=max(z.ab()/v0[i] for i,z in enumerate(E1));d2=max(z.ab()/v0[i] for i,z in enumerate(E2))
    l,h,el,eh=map(F,box);rad=max(abs(l-a0),abs(h-a0));q=q0+rad*q1
    if q>=1:raise AssertionError('not invertible by contraction')
    eps=(d0+rad*d1+rad**2*d2)/(1-q)
    margin=min(1-rad*abs(v1[i]/v0[i])-eps for i in range(n))
    if margin<=0:raise AssertionError('no positive weights')
    nominal={}
    for j,i in enumerate(ids):
        for (k,ap,ep),c0 in rows[i].items():
            if k==0:
                add(nominal,(ap,ep),c0*(v0[j]-a0*v1[j]));add(nominal,(ap+1,ep),c0*v1[j])
    for i,w in zip(free,fw):
        for (k,ap,ep),c0 in rows[i].items():
            if k==0:add(nominal,(ap,ep),c0*w)
    coefficients=bernstein(nominal,box);nom=max(t for r in coefficients for t in r)
    magnitude=F(0)
    for j,i in enumerate(ids):
        p={(a,e):c0 for (k,a,e),c0 in rows[i].items() if k==0};b=bernstein(p,box)
        magnitude+=v0[j]*max(abs(t) for r in b for t in r)
    objective=nom+eps*magnitude
    if objective>=0:raise AssertionError('objective not uniformly negative')
    return {'source_pencil_sha256':hashlib.sha256(Path(pencil).read_bytes()).hexdigest(),
        'certificate_sha256':hashlib.sha256(Path(cert).read_bytes()).hexdigest(),
        'arithmetic':'independent symmetric center-radius decimal grid 10^110; fused dot error bound; blossom Bernstein coefficients',
        'dimension':n,'atom_count':len(rows),'box':list(map(str,[l,h,el,eh])),
        'all_nuisance_rows_checked':len(keyall),'q0':str(q0),'q1':str(q1),'d0':str(d0),'d1':str(d1),'d2':str(d2),
        'contraction_upper':str(q),'weight_relative_error_upper':str(eps),'positive_weight_relative_lower':str(margin),
        'nominal_bernstein_coefficients':[[str(x) for x in r] for r in coefficients],
        'objective_error_upper':str(eps*magnitude),'objective_upper':str(objective),
        'passed':True,'proof_level':'finite exact/validated arithmetic; not an automated proof of physical premises'}

def reference_tests():
    rng=Random(40039);count=0
    for _ in range(500):
        n=rng.randrange(1,8);xx=[F(rng.randint(-10000,10000),rng.randint(1,999)) for _ in range(n)];yy=[F(rng.randint(-10000,10000),rng.randint(1,999)) for _ in range(n)]
        result=dot([Ball.scalar(x) for x in xx],[Ball.scalar(y) for y in yy]);l,h=result.bounds();v=sum(x*y for x,y in zip(xx,yy));assert l<=v<=h;count+=1
    for _ in range(100):
        poly={(p,q):F(rng.randint(-9,9),rng.randint(1,13)) for p in range(4) for q in range(3)}
        box=[F(-2,7),F(5,8),F(2,3),F(7,4)];b=bernstein(poly,box);t=F(rng.randrange(21),20);u=F(rng.randrange(21),20)
        exact=sum(c*(box[0]+t*(box[1]-box[0]))**p*(box[2]+u*(box[3]-box[2]))**q for (p,q),c in poly.items())
        rebuilt=sum(b[i][j]*comb(3,i)*t**i*(1-t)**(3-i)*comb(2,j)*u**j*(1-u)**(2-j) for i in range(4) for j in range(3))
        assert exact==rebuilt;count+=1
    return {'fused_dot_rational_references':500,'independent_bernstein_references':100,'total':count}

if __name__=='__main__':
    import argparse,time
    ap=argparse.ArgumentParser();ap.add_argument('--baseline',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args()
    cfg=json.loads((a.baseline/'inputs/parameters.json').read_text());a.output_dir.mkdir(parents=True,exist_ok=True)
    tests=reference_tests();(a.output_dir/'independent_interval_references.json').write_text(json.dumps(tests,sort_keys=True,indent=2)+'\n')
    for side in ('lower','upper'):
        t=time.time();r=audit(a.baseline/'inputs/radial_pencil.json',a.baseline/f'inputs/{side}_refined_continuation.json',cfg[side+'_box'])
        (a.output_dir/f'{side}_independent_audit.json').write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
        print(side,'q=',float(F(r['contraction_upper'])),'weights>',float(F(r['positive_weight_relative_lower'])),'objective<=',float(F(r['objective_upper'])),'seconds',time.time()-t,flush=True)
