"""Independent rational finite-N differential and color-centering references.
All coefficients lie in Q(i); exact Gaussian integration, not Monte Carlo.
Reference wavefunctions are REAL trial functions, not interacting eigenstates.
"""
from fractions import Fraction as F
from dataclasses import dataclass
from itertools import product
from math import prod
import json

@dataclass(frozen=True)
class C:
    x:F=F(0); y:F=F(0)
    def __post_init__(self):object.__setattr__(self,'x',F(self.x));object.__setattr__(self,'y',F(self.y))
    @staticmethod
    def of(z):return z if isinstance(z,C) else C(z)
    def __add__(self,z):z=C.of(z);return C(self.x+z.x,self.y+z.y)
    __radd__=__add__
    def __neg__(self):return C(-self.x,-self.y)
    def __sub__(self,z):return self+-C.of(z)
    def __rsub__(self,z):return C.of(z)+-self
    def __mul__(self,z):z=C.of(z);return C(self.x*z.x-self.y*z.y,self.x*z.y+self.y*z.x)
    __rmul__=__mul__
    def conj(self):return C(self.x,-self.y)
    def __bool__(self):return bool(self.x or self.y)
I=C(0,1)

def padd(a,b,c=C(1)):
    r=dict(a);c=C.of(c)
    for k,v in b.items():
        r[k]=r.get(k,C())+c*v
        if not r[k]:del r[k]
    return r

def pmul(a,b):
    r={}
    for u,c in a.items():
        for v,d in b.items():
            k=tuple(x+y for x,y in zip(u,v));r[k]=r.get(k,C())+c*d
            if not r[k]:del r[k]
    return r

def scale(p,c):return {k:C.of(c)*v for k,v in p.items() if C.of(c)*v}
def deriv(p,i):
    r={}
    for k,c in p.items():
        if k[i]:v=list(k);v[i]-=1;r[tuple(v)]=c*k[i]
    return r

def conj(p):return {k:c.conj() for k,c in p.items()}
def power(p,n,unit):
    r=unit
    for _ in range(n):r=pmul(r,p)
    return r

def mmul(A,B):
    n=len(A);m=len(B[0]);k=len(B)
    return [[sum_poly(pmul(A[i][l],B[l][j]) for l in range(k)) for j in range(m)] for i in range(n)]
def sum_poly(items):
    r={}
    for p in items:r=padd(r,p)
    return r

def mscale(A,c):return [[scale(p,c) for p in row] for row in A]
def trace(A):return scale(sum_poly(A[i][i] for i in range(len(A))),F(1,len(A)))
def frob(A):return scale(sum_poly(pmul(conj(p),p) for row in A for p in row),F(1,len(A)))

def model(N,M):
    # Real coordinates: diagonal weights1, real/imaginary off-diagonal weights2.
    dim=2*N*N;zero=(0,)*dim;unit={zero:C(1)};variables=[]
    for i in range(dim):e=list(zero);e[i]=1;variables.append({tuple(e):C(1)})
    mats=[];frames=[];metrics=[]
    for flavor in range(2):
        A=[[{} for j in range(N)] for i in range(N)];f=[];off=flavor*N*N;idx=off
        for i in range(N):
            T=[[C() for _ in range(N)] for _ in range(N)];T[i][i]=C(1);A[i][i]=variables[idx];f.append((idx,T,F(1)));metrics.append(1);idx+=1
        for i in range(N):
            for j in range(i+1,N):
                T=[[C() for _ in range(N)] for _ in range(N)];T[i][j]=T[j][i]=C(1)
                A[i][j]=padd(A[i][j],variables[idx]);A[j][i]=padd(A[j][i],variables[idx]);f.append((idx,T,F(1,2)));metrics.append(2);idx+=1
                T=[[C() for _ in range(N)] for _ in range(N)];T[i][j]=I;T[j][i]=-I
                A[i][j]=padd(A[i][j],variables[idx],I);A[j][i]=padd(A[j][i],variables[idx],-I);f.append((idx,T,F(1,2)));metrics.append(2);idx+=1
        mats.append(A);frames.append(f)
    zh=[[padd(mats[0][i][j],mats[1][i][j],I) for j in range(N)] for i in range(N)];bh=[[conj(p) for p in row] for row in zip(*zh)]
    R=scale(trace(mmul(zh,bh)),F(1,2));Chat=mmul(mmul(zh,zh),mmul(bh,bh));a4=scale(trace(Chat),F(1,4))
    variances=[F(1,2*N)*F(1,M)*F(1,g) for g in metrics]
    moment_cache={}
    def integrate(p):
        ans=C()
        for k,c in p.items():
            if k not in moment_cache:
                if any(t%2 for t in k):moment_cache[k]=F(0)
                else:moment_cache[k]=prod(prod(range(1,t,2))*variances[i]**(t//2) for i,t in enumerate(k))
            ans+=c*moment_cache[k]
        return ans
    def momentum(A,bar):
        out=[[{} for _ in range(N)] for _ in range(N)]
        for flavor in range(2):
            phase=C(1) if flavor==0 else (-I if bar else I)
            for idx,T,g in frames[flavor]:
                for k in range(N):
                    for j in range(N):
                        p=padd(deriv(A[k][j],idx),pmul(variables[idx],A[k][j]),-N*M*metrics[idx])
                        for i in range(N):
                            if T[i][k]:out[i][j]=padd(out[i][j],p,phase*T[i][k]*(-g/F(N)))
        return out
    def kinetic(p):
        out={}
        for i,g in enumerate(metrics):
            b=N*M*g
            x=variables[i]
            z=padd(deriv(deriv(p,i),i),pmul(x,deriv(p,i)),-2*b)
            z=padd(z,pmul(pmul(x,x),p),b*b);z=padd(z,p,-b)
            out=padd(out,z,F(-1,2*N*g))
        return out
    return unit,zh,bh,R,Chat,a4,integrate,momentum,kinetic

def run():
    records=[];ncontact=0;ncent=0;ncurrent=0
    for N,M,theta in [(1,F(1),F(0)),(1,F(3,2),F(1,3)),(2,F(1),F(0)),(2,F(3,2),F(0))]:
        print('finite differential case',N,M,theta,flush=True)
        one,Z,Zb,A,Chat,a4,integrate,P,Tkin=model(N,M)
        L=padd(one,A,theta);norm=integrate(pmul(L,L));assert norm.y==0 and norm.x>0
        psi=[[L if i==j else {} for j in range(N)] for i in range(N)]
        Pb=P(psi,True);D=mmul(mmul(Z,Z),Pb);K=mscale(P(Pb,False),F(-1,2))
        S=mscale(P(mmul(mmul(Zb,Zb),D),False),F(-1,8))
        CK=mscale(mmul(Chat,K),F(1,4))
        J=trace(D)
        Fhat=trace(mmul(mmul(Z,Z),Zb));Khat=trace(mmul(mmul(Z,Zb),P(psi,False)))
        comm=padd(Tkin(pmul(Fhat,L)),pmul(Fhat,Tkin(L)),-1)
        residual=padd(padd(padd(scale(Khat,2),J),comm,-1),pmul(trace(Z),L),-4)
        if residual:raise AssertionError('finite-N cubic derivative/gauge current identity')
        ncurrent+=1
        ctr=scale(pmul(trace(Zb),J),F(1,4))
        zz=scale(pmul(trace(Zb),trace(Z)),F(1,2))
        # Exact operator color-centering, norm on the actual test vector.
        center=[[padd(D[i][j],J,-1) if i==j else D[i][j] for j in range(N)] for i in range(N)]
        lhs=frob(center);rhs=padd(frob(D),pmul(conj(J),J),-1)
        assert lhs==rhs;value=integrate(lhs);assert value.y==0 and value.x>=0;ncent+=1
        for m in (range(3) if N==1 else range(2) if M==1 else range(1)):
            print('weighted contact',N,m,flush=True)
            f2=power(A,2*m,one)
            target=integrate(pmul(pmul(f2,L),trace(S)))
            q=integrate(pmul(pmul(f2,L),trace(CK)))
            radius=integrate(pmul(pmul(pmul(f2,L),L),A))*F(1,2)
            centerpair=integrate(pmul(pmul(pmul(f2,L),L),zz))*F(1,2)
            centerflow=integrate(pmul(pmul(f2,L),ctr))
            derivterm=C()
            if m:derivterm=integrate(pmul(pmul(pmul(power(A,2*m-1,one),L),L),a4))*F(m,N*N)
            assert target==q+radius+centerpair+centerflow+derivterm,(N,M,theta,m,target,q+radius+centerpair+centerflow+derivterm)
            records.append({'N':N,'M':str(M),'radial_trial_polynomial':f'1+({theta})A','weight_power':m,
                'left_weighted_correlation_DdagD_not_a_norm':str(target.x/norm.x),
                'contact_A_over2':str(radius.x/norm.x),
                'center_trace_pair':str(centerpair.x/norm.x),
                'center_flow':str(centerflow.x/norm.x),
                'weight_derivative_contact':str(derivterm.x/norm.x)})
            ncontact+=1
    return {'finite_N_weighted_contact_identities':ncontact,'operator_color_projection_identities':ncent,'finite_N_cubic_current_action_identities':ncurrent,
      'records':records,'reference_states':'real radial polynomial times harmonic Gaussian; not interacting ground states',
      'normalization':'Pi=-(grad_X+i grad_Y)/(sqrt(2)N), Tr/N; offdiagonal real coordinates have metric2'}

if __name__=='__main__':
    from pathlib import Path
    import argparse,time
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();t=__import__('time').time();r=run();a.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n');print('finite contacts',r['finite_N_weighted_contact_identities'],'centered operator identities',r['operator_color_projection_identities'],'seconds',__import__('time').time()-t)
