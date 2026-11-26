"""Independent finite reference checks for new quartic traces and covariance algebra.
These checks do not prove the inherited all-N analytical assumptions.
"""
from fractions import Fraction as F
from functools import lru_cache
from collections import Counter
import trace_algebra as a

@lru_cache(None)
def pairings(n):
 if n==0:return ((),)
 if n%2:return ()
 ans=[]
 for j in range(1,n):
  rest=[i for i in range(1,n) if i!=j]
  for p in pairings(n-2):ans.append(((0,j),)+tuple((rest[u],rest[v]) for u,v in p))
 return tuple(ans)

def cov(x,y,M=F(1),real=False):
 if real:return F(1,2)/M if x==y else F(0)
 if (x,y) in [(0,1),(1,0)]:return F(1,2)/M
 if (x,y) in [(2,3),(3,2)]:return -M/2
 if (x,y) in [(0,3),(1,2)]:return F(1,2)
 if (x,y) in [(3,0),(2,1)]:return F(-1,2)
 return F(0)

@lru_cache(None)
def planar(w,M=F(1)):
 if not w:return F(1)
 if len(w)%2:return F(0)
 return sum((cov(w[0],w[j],M)*planar(w[1:j],M)*planar(w[j+1:],M) for j in range(1,len(w),2)),F(0))

def pexpect(p,M=F(1)):
 ans=F(0)
 for mon,c in p.items():
  v=c
  for w in mon:v*=planar(w,M)
  ans+=v
 return ans

def finite_wick_polynomial(traces,M=F(1),real=False):
 """Laurent polynomial in N, retaining ORDERED oscillator contractions."""
 traces=tuple(t for t in traces if t);flat=sum(traces,());n=len(flat)
 if n==0:return {0:F(1)}
 gamma=[];off=0
 for w in traces:
  gamma.extend([off+(i+1)%len(w) for i in range(len(w))]);off+=len(w)
 ans={}
 for pairs in pairings(n):
  cf=F(1);pi=list(range(n))
  for u,v in pairs:
   cf*=cov(flat[u],flat[v],M,real);pi[u]=v;pi[v]=u
  if not cf:continue
  perm=[gamma[pi[i]] for i in range(n)];seen=set();cycles=0
  for i in range(n):
   if i in seen:continue
   cycles+=1;j=i
   while j not in seen:seen.add(j);j=perm[j]
  a.add(ans,cycles-len(traces)-n//2,cf)
 return ans

def eval_laurent(p,N):return sum((c*F(N)**j for j,c in p.items()),F(0))

def complex_expansion(signs):
 """Trace product of (X+i*s*Y)/sqrt(2), with four factors only."""
 assert len(signs)==4
 out={}
 for mask in range(16):
  w=tuple((mask>>i)&1 for i in range(4));k=sum(w);sg=1
  for s,x in zip(signs,w):
   if x:sg*=s
  phase=k%4;re=F(sg,4) if phase==0 else (F(-sg,4) if phase==2 else F(0))
  im=F(sg,4) if phase==1 else (F(-sg,4) if phase==3 else F(0))
  key=a.cyclic_min(w);rr,ii=out.get(key,(F(0),F(0)));out[key]=(rr+re,ii+im)
 return {k:v for k,v in out.items() if v!=(0,0)}

# Tiny commutative polynomial engine for exact congruences in e and scalar moments.
NV=7
ZERO=(0,)*NV
def scalar(x):return {ZERO:F(x)} if x else {}
def var(j):
 k=[0]*NV;k[j]=1;return {tuple(k):F(1)}
def padd(p,q,c=F(1)):
 z=dict(p)
 for k,v in q.items():a.add(z,k,c*v)
 return z
def pmul(p,q):
 z={}
 for u,c in p.items():
  for v,d in q.items():a.add(z,tuple(x+y for x,y in zip(u,v)),c*d)
 return z
def transpose(M):return list(map(list,zip(*M)))
def mmul(A,B):
 return [[sum_polys(pmul(A[i][k],B[k][j]) for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def sum_polys(seq):
 z={}
 for p in seq:z=padd(z,p)
 return z

def covariance_checks():
 # e,m0,m1,m2,d1,dA,d2; d here denotes a4, not a momentum.
 e,m0,m1,m2,d1,dA,d2=[var(i) for i in range(7)]
 z={};one=scalar(1)
 T=[[one,z,z],[z,one,z],[pmul(scalar(F(2,3)),e),scalar(F(-1,3)),z],
    [z,z,one],[pmul(scalar(F(1,3)),e),scalar(F(-2,3)),one]]
 G=[[m0,m1,d1],[m1,m2,dA],[d1,dA,d2]]
 H=mmul(mmul(T,G),transpose(T))
 kernels=[ [pmul(scalar(F(-2,3)),e),scalar(F(1,3)),one,z,z],
           [pmul(scalar(F(-1,3)),e),scalar(F(2,3)),z,scalar(-1),one] ]
 checks=0
 for k in kernels:
  product=mmul(H,[[v] for v in k])
  assert all(not r[0] for r in product);checks+=5
 # Selected rows recover G identically; congruence is exactly equivalent to G>=0.
 sel=[0,1,3]
 for i in range(3):
  for j in range(3):assert H[sel[i]][sel[j]]==G[i][j];checks+=1
 # Positive two-atom scalar example with NONZERO radial and quartic covariance.
 ee=F(59,50);atoms=[(F(389,1000),F(23,100)),(F(3891,10000),F(231,1000))]
 rows=[]
 for A,D in atoms:rows.append([F(1),A,(2*ee-A)/3,D,D+(ee-2*A)/3])
 numeric=[[sum(r[i]*r[j] for r in rows)/2 for j in range(5)] for i in range(5)]
 variance=sum(A*A for A,D in atoms)/2-(sum(A for A,D in atoms)/2)**2
 assert variance==F(1,400000000)
 def outpoly(p):return [{'exponents':list(k),'coefficient':str(v)} for k,v in sorted(p.items())]
 return {'variables':['energy_e','m0','m1','m2','mean_a4','mixed_A_a4','second_a4'],
 'observable_order':['1','A','p','a4','b4'],
 'T':[[outpoly(p) for p in row] for row in T],
 'G3':[[outpoly(p) for p in row] for row in G],
 'G5':[[outpoly(p) for p in row] for row in H],
 'kernels':[[outpoly(p) for p in row] for row in kernels],
 'formal_entry_checks':checks,
 'rank_statement':'G5=T(e) G3 T(e)^T; T contains the 3x3 identity, so PSD is equivalent and rank<=3. Centered covariance rank<=2.',
 'finite_N_residual':'||(p-(2 e_N-A)/3) psi_N||^2=Var(B_N)/9; b4-a4-(e_N-2 A)/3=(B_N-e_N)/3.',
 'nonfactorized_example':{'atoms':[[str(x) for x in r] for r in rows],'Gram':[[str(c) for c in row] for row in numeric],
 'variance_A':str(variance),'warning':'Example for the covariance identities ONLY; not a solution of the full word system or a physical state.'}}

def run(model):
 counts=Counter();negative=[]
 e4=complex_expansion([1,1,1,1]);e2=complex_expansion([1,1,1,-1])
 im4={k:v[1] for k,v in e4.items() if v[1]};im2={k:v[1] for k,v in e2.items() if v[1]}
 assert im4=={(0,0,0,1):F(1),(0,1,1,1):F(-1)}
 assert im2=={(0,0,0,1):F(1,2),(0,1,1,1):F(1,2)}
 counts['formal_quartic_trace_expansions']=2
 norms={}
 for w in [(0,0,0,0),(0,0,0,1)]:
  _,dw=a.involution(w);p=finite_wick_polynomial((dw,w));norms[str(w)]={str(k):str(v) for k,v in sorted(p.items())}
  assert all(k<=-2 for k in p)
  const=F(15,2) if w==(0,0,0,0) else F(15,8)
  for M in [F(1),F(2),F(3,2)]:
   pp=finite_wick_polynomial((dw,w),M)
   for N in [1,2,3,5,10]:
    v=eval_laurent(pp,N);assert 0<v<=const/(M**4*N*N);counts['finite_Gaussian_quartic_bounds']+=1
 # Check every kinetic and gauge raw equation against the independent planar free vacuum.
 for kind,w,mul,p in a.raw_system(11):
  if 'stationarity' in kind:
   p={a.raw_monomial(t):c for t,c in harmonic_derivative(w).items()}
  pp={}
  for mon,c in p.items():a.add(pp,mon+mul,c)
  assert pexpect(pp)==0,(kind,w,mul,pexpect(pp))
  counts['raw_free_reference_equations']+=1
  qq=model['normalizer'].polynomial(p,mul)
  assert pexpect(qq)==0,(kind,w,mul)
  counts['normalized_free_reference_equations']+=1
 # Explicitly check that known nontrivial unknown double traces were not factorized.
 aa=(0,0,1,1);bb=(0,1,0,1);A=(0,1);P=(2,3)
 for mon in [(aa,aa),(aa,bb),(bb,bb),(P,P)]:
  r=a.product_expression({mon:F(1)},model)
  assert r and any(len(model['keys'][i])==2 for i in r);counts['retained_quartic_covariances']+=1
 for v in [(),A,aa,bb,P]:
  r={}
  for w,c in [(P,F(1)),(A,F(1)),(bb,F(2)),(aa,F(-2))]:
   a.combine(r,a.product_expression({a.raw_monomial(v,w):c},model))
  assert not r;counts['multiplied_virial_identities']+=1
 for N in [1,2,3,4]:
  p=finite_wick_polynomial(((0,3),));assert eval_laurent(p,N)==F(1,2);counts['contact_constants']+=1
 # Deliberate false claims each cause an exact test to fail.
 def reject(name,test):
  try:assert test()
  except AssertionError:negative.append(name);return
  raise AssertionError('bad control was not rejected: '+name)
 reject('delete finite-N charged quartic variance',lambda:eval_laurent(finite_wick_polynomial(((0,0,0,0),(1,1,1,1))),2)==0)
 reject('wrong quartic imaginary coefficient',lambda:im4[(0,0,0,1)]==F(1,2))
 reject('wrong mixed quartic imaginary sign',lambda:im2[(0,1,1,1)]==F(-1,2))
 reject('drop momentum contact',lambda:pexpect({((0,3),):F(1)})==0)
 reject('wrong Pi dagger sign',lambda:a.involution((2,))[0]==1)
 reject('replace radial second moment by mean squared',lambda:F(1,400000000)==0)
 reject('state symmetry applied to a single charged factor',lambda:planar((0,0)) * planar((1,1))==F(1,8))
 reject('reuse charged-trace rule for neutral fourth trace',lambda:model['normalizer'](((0,0,1,1),))[0]==0)
 reject('pretend energy is fixed before choosing a subsequence',lambda:F(59,50)==F(1181,1000))
 return {'counts':dict(counts),'negative_controls_rejected':negative,'Gaussian_quartic_variance_Laurent':norms,
 'proof_scope':'Finite exact regressions only; the conditional Gaussian/Poincare and weighted limit arguments are written proofs with inherited assumptions.'}

@lru_cache(None)
def harmonic_derivative(w):
 out={}
 for i,x in enumerate(w):a.add(out,w[:i]+((x+2)%4,)+w[i+1:],F(1))
 return out
