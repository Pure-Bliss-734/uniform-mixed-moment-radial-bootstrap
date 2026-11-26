"""Outward rounded fixed dyadic intervals. Standard library proof arithmetic.
No empirical floating point tolerance enters accepted inequalities.
"""
from fractions import Fraction as F
PRECISION=320
SCALE=1<<PRECISION

def updiv(a,b):return -((-a)//b)
class I:
    __slots__=('lo','hi')
    def __init__(self,x=0,hi=None,raw=False):
        if raw:self.lo=int(x);self.hi=int(hi)
        elif isinstance(x,I):self.lo=x.lo;self.hi=x.hi
        else:
            a=F(x);b=a if hi is None else F(hi)
            self.lo=(a.numerator*SCALE)//a.denominator
            self.hi=updiv(b.numerator*SCALE,b.denominator)
        if self.lo>self.hi:raise ValueError('empty interval')
    def __add__(self,b):
        b=b if isinstance(b,I) else I(b);return I(self.lo+b.lo,self.hi+b.hi,raw=True)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo,raw=True)
    def __sub__(self,b):return self+-I(b)
    def __rsub__(self,b):return I(b)+-self
    def __mul__(self,b):
        b=b if isinstance(b,I) else I(b)
        p=(self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi)
        return I(min(p)//SCALE,updiv(max(p),SCALE),raw=True)
    __rmul__=__mul__
    def inv(self):
        if self.lo<=0<=self.hi:raise ZeroDivisionError('interval contains zero')
        return I((SCALE*SCALE)//self.hi,updiv(SCALE*SCALE,self.lo),raw=True)
    def __truediv__(self,b):return self*I(b).inv()
    def __rtruediv__(self,b):return I(b)*self.inv()
    def absupper(self):return F(max(abs(self.lo),abs(self.hi)),SCALE)
    def lower(self):return F(self.lo,SCALE)
    def upper(self):return F(self.hi,SCALE)
    def to_json(self):return [str(self.lower()),str(self.upper())]
    def __repr__(self):return repr(self.to_json())

def matmul(A,B):
    n=len(A);k=len(B);m=len(B[0]);out=[]
    # Fuse the sum: do only a single outward division per dot product.
    cols=list(zip(*B))
    for row in A:
        rr=[]
        for col in cols:
            l=h=0
            for a,b in zip(row,col):
                pp=(a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi)
                l+=min(pp);h+=max(pp)
            rr.append(I(l//SCALE,updiv(h,SCALE),raw=True))
        out.append(rr)
    return out

def mv(A,v):return [r[0] for r in matmul(A,[[x] for x in v])]
def matrix(raw):return [[I(F(v)) for v in row] for row in raw]
