"""Finite exact references for Stage39. These are not a proof assistant for disintegration."""
from fractions import Fraction as F
from random import Random
from math import comb
from interval_arithmetic import I,matrix,matmul,mv
from continuation_certificate import bernstein
from radial_model import psd_pivots

def expect(atoms,fun):return sum(w*fun(a,b,j) for w,a,b,j in atoms)
def gram(atoms,fun):
    vec=[(w,fun(a,b,j)) for w,a,b,j in atoms];n=len(vec[0][1])
    return [[sum(w*v[i]*v[j] for w,v in vec) for j in range(n)] for i in range(n)]

def conditioning_reference():
    # Two distinct radius fibers and nonzero neutral fluctuations inside EACH fiber.
    atoms=[]
    for a,b0 in [(F(389,1000),F(23,100)),(F(3891,10000),F(231,1000))]:
        for d in [-F(1,1000),F(1,1000)]:atoms.append((F(1,4),a,b0+d,F(1,10)))
    G=gram(atoms,lambda a,b,j:[F(1),a,b]);piv=psd_pivots(G)
    if len(piv)!=3:raise AssertionError('unconditional neutral Gram was factored')
    ma=expect(atoms,lambda a,b,j:a);va=expect(atoms,lambda a,b,j:a*a)-ma*ma
    mb=expect(atoms,lambda a,b,j:b);vb=expect(atoms,lambda a,b,j:b*b)-mb*mb
    cov=expect(atoms,lambda a,b,j:a*b)-ma*mb
    if va!=F(1,400000000) or cov!=F(1,40000000):raise AssertionError('unconditional covariance mismatch')
    fibers=[];tower=0
    for a in sorted(set(a for w,a,b,j in atoms)):
        sub=[(2*w,aa,b,j) for w,aa,b,j in atoms if aa==a]
        g=gram(sub,lambda a,b,j:[F(1),a,b]);C=gram(sub,lambda a,b,j:[F(1),b]);K=[-a,F(1),F(0)]
        if any(sum(g[i][j]*K[j] for j in range(3)) for i in range(3)):raise AssertionError('fiber kernel failed')
        if len(psd_pivots(C))!=2:raise AssertionError('conditional a4 variance deleted')
        T=[[F(1),F(0)],[a,F(0)],[F(0),F(1)]]
        back=[[sum(T[i][k]*C[k][l]*T[j][l] for k in range(2) for l in range(2)) for j in range(3)] for i in range(3)]
        if g!=back:raise AssertionError('exact Gram congruence failed')
        for n in range(7):
            for m in range(5):
                lhs=expect(sub,lambda aa,b,j:aa**n*b**m)
                rhs=a**n*expect(sub,lambda aa,b,j:b**m)
                if lhs!=rhs:raise AssertionError('fiber multiplication')
                tower+=1
        fibers.append({'a':str(a),'neutral_conditional_mean':str(C[0][1]),'neutral_conditional_variance':str(C[1][1]-C[0][1]**2),
          'flow_norm':'1/10','G2':[[str(v) for v in row] for row in C],
          'G2_positive_pivots':[str(v) for v in psd_pivots(C)],'G3_kernel':[str(v) for v in K]})
    # Conditional coefficient vectors may depend on radius without becoming a global constant.
    identities=0
    for p in range(8):
        f=lambda a:a**p+F(1,3)
        lhs=expect(atoms,lambda a,b,j:f(a)**2*b)
        rhs=sum(F(1,2)*f(F(z['a']))**2*F(z['neutral_conditional_mean']) for z in fibers)
        if lhs!=rhs:raise AssertionError('weighted recombination failed')
        identities+=1
    return {'atoms':[[str(x) for x in row] for row in atoms],
       'unconditional_G3':[[str(v) for v in row] for row in G],'unconditional_G3_positive_pivots':[str(v) for v in piv],
       'Var_radius':str(va),'Var_neutral_quartic':str(vb),'Cov_radius_quartic':str(cov),'fibers':fibers,
       'fiber_multiplication_checks':tower,'weighted_recombination_checks':identities,
       'scope':'Only the conditional/energy Gram algebra is tested. This distribution is NOT claimed feasible for the full matrix-word system or to be a quantum state.'}

def arithmetic_references():
    rng=Random(390039);counts={'interval_operations':0,'matrix_products':0,'bernstein_evaluations':0,'variable_weight_continuation_examples':0}
    for _ in range(200):
        a=F(rng.randrange(-10000,10001),rng.randrange(1,999));b=F(rng.randrange(-10000,10001),rng.randrange(1,999))
        A=I(a);B=I(b)
        for z,v in [(A+B,a+b),(A-B,a-b),(A*B,a*b)]+([(A/B,a/b)] if b else []):
            if not z.lower()<=v<=z.upper():raise AssertionError('outward interval failed')
            counts['interval_operations']+=1
    for n in range(1,6):
        for _ in range(6):
            A=[[F(rng.randrange(-9,10),rng.randrange(1,20)) for j in range(n)] for i in range(n)]
            B=[[F(rng.randrange(-9,10),rng.randrange(1,20)) for j in range(n)] for i in range(n)]
            Z=matmul(matrix(A),matrix(B))
            for i in range(n):
                for j in range(n):
                    v=sum(A[i][k]*B[k][j] for k in range(n))
                    if not Z[i][j].lower()<=v<=Z[i][j].upper():raise AssertionError('interval matrix product')
                    counts['matrix_products']+=1
    for _ in range(25):
        p={(i,j):F(rng.randrange(-20,21),rng.randrange(1,17)) for i in range(4) for j in range(3)}
        lo=F(-2,5);hi=F(4,7);el=F(9,8);eh=F(6,5);B=bernstein(p,[lo,hi,el,eh])
        for u,v in [(F(0),F(0)),(F(1),F(1)),(F(1,3),F(2,5)),(F(7,10),F(4,9))]:
            exact=sum(c*(lo+(hi-lo)*u)**i*(el+(eh-el)*v)**j for (i,j),c in p.items())
            got=sum(B[i][j]*comb(3,i)*u**i*(1-u)**(3-i)*comb(2,j)*v**j*(1-v)**(2-j) for i in range(4) for j in range(3))
            if exact!=got:raise AssertionError('Bernstein formula')
            counts['bernstein_evaluations']+=1
    # A(a) has varying exact inverse: never replace a-dependent weights by their center.
    for n in range(1,51):
        a=F(n,1000);M=[[F(1),a],[a,F(1)]];b=[F(1),F(2)];det=1-a*a
        w=[(1-2*a)/det,(2-a)/det]
        if min(w)<=0 or any(sum(M[i][j]*w[j] for j in range(2))!=b[i] for i in range(2)):raise AssertionError('parametric exact weights')
        counts['variable_weight_continuation_examples']+=1
    return counts
