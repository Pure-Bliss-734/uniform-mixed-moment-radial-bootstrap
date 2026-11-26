"""Unfactorized trace algebra through weight eleven.
Whole-product symmetries; no product of unknown traces is factorized.
All exact arithmetic uses fractions.Fraction.
"""
from fractions import Fraction as F
from functools import lru_cache
from collections import Counter
WEIGHTS=(1,1,2,2); CHARGES=(1,-1,1,-1); CONJ=(3,2,1,0)
LABELS=('Z','Zbar','Pi','Pibar')
def weight(w):return sum(WEIGHTS[x] for x in w)
def charge(w):return sum(CHARGES[x] for x in w)
def add(p,k,c):
 if c:
  p[k]=p.get(k,F(0))+c
  if not p[k]:del p[k]
def combine(p,q,c=1):
 for k,v in q.items():add(p,k,c*v)
 return p
@lru_cache(None)
def words(L,q=None):
 if L<0:return ()
 ans=[]
 def visit(w,k,c):
  if q is None or c==q:ans.append(w)
  for x in range(4):
   if k+WEIGHTS[x]<=L:visit(w+(x,),k+WEIGHTS[x],c+CHARGES[x])
 visit((),0,0)
 return tuple(sorted(ans,key=lambda w:(weight(w),len(w),w)))
def coordinate(w):return all(x<2 for x in w)
def cyclic_min(w):return min(w[i:]+w[:i] for i in range(len(w))) if w else w
def raw_monomial(*ws):return tuple(w for w in ws if w)
@lru_cache(None)
def differentiate(w):
 force=[{(2,):F(1)},{(3,):F(1)},
   {(0,):F(1),(0,1,0):F(2),(0,0,1):F(-1),(1,0,0):F(-1)},
   {(1,):F(1),(1,0,1):F(2),(1,1,0):F(-1),(0,1,1):F(-1)}]
 ans={}
 for i,x in enumerate(w):
  for t,c in force[x].items():add(ans,w[:i]+t+w[i+1:],c)
 return ans
def involution(w):return (-1)**sum(x>=2 for x in w),tuple(x^1 for x in w[::-1])
def raw_cyclic(w):
 if not w:return {}
 x,v=w[0],w[1:];p={}
 add(p,raw_monomial(w),F(1));add(p,raw_monomial(v+(x,)),F(-1))
 for i,y in enumerate(v):
  if y==CONJ[x]:add(p,raw_monomial(v[:i],v[i+1:]),F(-1 if x<2 else 1))
 return p
def raw_gauge(w):
 p={}
 for tail,c in [((0,3),1),((3,0),-1),((1,2),1),((2,1),-1)]:add(p,raw_monomial(w+tail),F(c))
 add(p,raw_monomial(w),F(-2));return p
class TraceNormalizer:
 def __init__(self,charged_four=False):self.charged_four=charged_four;self.cached={}
 def __call__(self,m):
  m=tuple(tuple(w) for w in m if w)
  if m in self.cached:return self.cached[m]
  c=F(1);ws=[];reasons=[]
  for w in m:
   if len(w)==1:
    out=(F(0),(),('center',));self.cached[m]=out;return out
   if coordinate(w):
    if len(w)%2:
     out=(F(0),(),('odd_coordinate_trace',));self.cached[m]=out;return out
    if len(w)==2 and charge(w):
     out=(F(0),(),('charged_coordinate2',));self.cached[m]=out;return out
    if self.charged_four and len(w)==4 and charge(w):
     out=(F(0),(),('charged_coordinate4',));self.cached[m]=out;return out
    w=cyclic_min(w)
   if len(w)==2 and sum(x>=2 for x in w)==1:
    if charge(w):
     out=(F(0),(),('charged_bilinear_current',));self.cached[m]=out;return out
    c*=F(1,2) if w[0]<2 else F(-1,2);reasons.append('neutral_current_contact');continue
   ws.append(w)
  ws=tuple(sorted(ws));options=[]
  # Reflection/time reversal applied to ALL factors together.
  for flip in (False,True):
   for rev in (False,True):
    sign=(-1)**sum(x>=2 for w in ws for x in w) if rev else 1;out=[]
    for w in ws:
     v=w[::-1] if rev else w
     if flip:v=tuple(x^1 for x in v)
     if coordinate(v):v=cyclic_min(v)
     out.append(v)
    options.append((tuple(sorted(out)),sign))
  key=min(a for a,s in options);signs={s for a,s in options if a==key}
  out=(F(0),(),('odd_whole_symmetry',)) if len(signs)==2 else (c*next(iter(signs)),key,tuple(reasons))
  self.cached[m]=out;return out
 def polynomial(self,p,multiplier=()):
  out={}
  for m,c in p.items():
   d,k,_=self(m+multiplier);add(out,k,c*d)
  return out
class Echelon:
 def __init__(self):self.pivots={}
 def insert(self,r):
  r=dict(r)
  while r:
   j=min(r)
   if j not in self.pivots:
    a=r[j];self.pivots[j]={k:c/a for k,c in r.items()};return j
   combine(r,self.pivots[j],-r[j])
  return None
 def reduce(self,r):
  r=dict(r);out={}
  while r:
   j=min(r)
   if j in self.pivots:combine(r,self.pivots[j],-r[j])
   else:out[j]=r.pop(j)
  return out

def raw_system(L):
 for w in words(L,0):
  p=raw_cyclic(w)
  if p:yield 'single_cut',w,(),p
 for w in words(L-3,0):yield 'single_gauge',w,(),raw_gauge(w)
 for w in words(L-1,0):
  p={raw_monomial(t):c for t,c in differentiate(w).items()}
  if p:yield 'single_stationarity',w,(),p
 B=L-3
 for v in words(B-1):
  if not v:continue
  top=B-weight(v)
  for w in words(top,-charge(v)):
   if not w:continue
   p=raw_cyclic(w)
   if p:yield 'product_cut',w,(v,),p
  for w in words(top-3,-charge(v)):yield 'product_gauge',w,(v,),raw_gauge(w)
  for w in words(top-1,-charge(v)):
   p={raw_monomial(t):c for t,c in differentiate(w).items()}
   if p:yield 'product_stationarity',w,(v,),p

def mon_key(m):return (1 if all(w==(0,1) for w in m) else 0,len(m),sum(weight(w) for w in m),m)
def build(L=11,charged_four=False,progress=False,verify_equations=True):
 import time
 t0=time.time();norm=TraceNormalizer(charged_four);eq=[];keys=set();counts=Counter();cuts=[]
 for w in words(L,0):
  c,m,_=norm(raw_monomial(w))
  if c:keys.add(m)
 for kind,w,mul,p in raw_system(L):
  r=norm.polynomial(p,mul);counts[kind]+=1
  if kind=='single_cut':
   for mon,c in p.items():
    if len(mon)==2:
     d,key,reasons=norm(mon)
     cuts.append({'word':list(w),'factors':[list(u) for u in mon],'coefficient':str(c),
       'multiplier':str(d),'reduced_factors':[list(u) for u in key], 'reason':list(reasons) if reasons else ['retained']})
  keys.update(r)
  if r:eq.append((kind,w,mul,r))
 eq.sort(key=lambda row: (len(row[3]), max((sum(weight(w) for w in mon) for mon in row[3]),default=0)))
 keys=sorted(keys,key=mon_key);index={k:i for i,k in enumerate(keys)};E=Echelon()
 for it,(_,_,_,p) in enumerate(eq):
  E.insert({index[k]:v for k,v in p.items()})
  if progress and it%5000==0:print('rows',it,'rank',len(E.pivots),'seconds',time.time()-t0,flush=True)
 free=[i for i in range(len(keys)) if i not in E.pivots]
 if progress: print('forward complete',len(E.pivots),'seconds',time.time()-t0,flush=True)
 expr={}
 for i in reversed(range(len(keys))):
  if i not in E.pivots:expr[i]={i:F(1)}
  else:
   r={}
   for j,c in E.pivots[i].items():
    if j!=i:combine(r,expr[j],-c)
   expr[i]=r
 if progress: print('back substitution complete','seconds',time.time()-t0,flush=True)
 if verify_equations:
  for kind,w,mul,p in eq:
   r={}
   for k,c in p.items():combine(r,expr[index[k]],c)
   if r:raise AssertionError(('nonzero equation',kind,w,mul,r))
 return {'L':L,'normalizer':norm,'keys':keys,'index':index,'E':E,'free':free,'expressions':expr,
   'equations':eq,'counts':dict(counts),'cuts':cuts,'elapsed':time.time()-t0}
def trace_expression(p,m):
 r={}
 for w,c in p.items():
  d,k,_=m['normalizer'](raw_monomial(w))
  if d:combine(r,m['expressions'][m['index'][k]],c*d)
 return r
def product_expression(p,m):
 r={}
 for w,c in p.items():
  d,k,_=m['normalizer'](w)
  if d:combine(r,m['expressions'][m['index'][k]],c*d)
 return r
def gram(u,v,ground=False):
 out={}
 for w,c in u.items():
  s,dw=involution(w)
  for z,e in v.items():
   for t,g in (differentiate(z) if ground else {z:F(1)}).items():add(out,dw+t,s*c*e*g)
 return out
def bases(order):
 ans=[];seen=set();even=[];odd=[]
 for w in words(order,0):
  if w in seen:continue
  v=tuple(a^1 for a in w);seen|={w,v}
  if v==w:even.append({w:F(1)})
  else:even.append({w:F(1,2),v:F(1,2)});odd.append({w:F(1,2),v:F(-1,2)})
 ans.extend([('q0_even',even),('q0_odd',odd)])
 for q in range(1,order+1):
  base=[{w:F(1)} for w in words(order,q)]
  if base:ans.append(('q'+str(q),base))
 return ans
def blocks(m,order=5):
 ans=[]
 for ground in (False,True):
  for name,basis in bases(order):
   ent=[[None for v in basis] for u in basis]
   for i,u in enumerate(basis):
    for j in range(i,len(basis)):
     r=trace_expression(gram(u,basis[j],ground),m)
     rev=trace_expression(gram(basis[j],u,ground),m)
     if r!=rev:raise AssertionError(('asymmetric block',ground,name,i,j))
     ent[i][j]=r;ent[j][i]=r
   ans.append({'name':('N_' if ground else 'M_')+name,'basis':basis,'entries':ent})
 return ans
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--level',type=int,default=11);ap.add_argument('--charged-four',action='store_true')
 a=ap.parse_args();m=build(a.level,a.charged_four,True)
 print('COUNTS',len(words(a.level,0))-1,len(m['keys']),len(m['equations']),len(m['E'].pivots),len(m['free']),m['counts'],m['elapsed'])
 print('FREE',[(i,m['keys'][i]) for i in m['free']])
 bs=blocks(m,(a.level-1)//2);print('BLOCKS',[(b['name'],len(b['basis'])) for b in bs])
