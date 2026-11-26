"""Extract the directed sigma/current budgets from actual conditional matrices."""
from fractions import Fraction as F
from collections import Counter
from random import Random
from trace_images import add,encoded
from augment_model import decode,axpy

def parity_part(word,c,parity):
    out={};add(out,tuple(word),F(c,2));add(out,tuple(x^1 for x in word),F(c*parity,2));return out

def find_vector(basis,wanted):
    if not wanted:return {}
    for i,poly in enumerate(basis):
        p={tuple(w):F(c) for w,c in poly}
        if set(p)!=set(wanted):continue
        k=min(p);r=wanted[k]/p[k]
        if all(wanted[j]==r*p[j] for j in p):return {i:r}
    raise AssertionError(('exact parity representation absent',wanted))

def extract(images):
    byname={b['name']:b for b in images['blocks']};G=[[{} for _ in range(2)] for _ in range(2)]
    for name,sg in [('C_q0_even',1),('C_q0_odd',-1)]:
        b=byname[name];E=[[decode(r) for r in row] for row in b['entries']]
        vectors=[find_vector(b['basis'],parity_part((0,0,1,1),1,sg)),find_vector(b['basis'],parity_part((2,3),-1,sg))]
        for i in range(2):
            for j in range(2):
                for u,c in vectors[i].items():
                    for v,d in vectors[j].items():axpy(G[i][j],E[u][v],c*d)
    if G[0][0]!={(18,0,0):F(1),(2,0,0):F(-1)}:raise AssertionError('neutral budget')
    expected={(19,0,0):F(-1),(1,0,1):F(-2,3),(1,1,0):F(1,3)}
    if G[0][1]!=expected or G[1][0]!=expected:raise AssertionError('shared q/p/u coefficient')
    b=byname['C_q1'];idx=next(i for i,p in enumerate(b['basis']) if p==[[[0,0,3],'1']])
    flow=decode(b['entries'][idx][idx]);expected_flow={(0,1,0):F(1,2),(19,0,0):F(-1),(3,0,0):F(-1)}
    if flow!=expected_flow:raise AssertionError('J contact budget')
    # A nonfactorized local reference. Not claimed to satisfy the full SDP.
    a,e,u,sigma,j=F(389,1000),F(117,100),F(23,100),F(1,10**6),F(1,10)
    p=(2*e-a)/3;R=u*u+sigma+1;K=p*p+1;q=p*u+F(1,10);v=u*u+sigma
    assert sigma>0 and j>0 and R-v>0 and K-p*p>0
    det=(R-v)*(K-p*p)-(q-p*u)**2
    assert det==F(99,100) and a/2+q-j>0
    return {'notation':{'C':'Z^2 Zbar^2','K':'-Pi Pibar','D':'Z^2 Pibar','p':'(2e-a)/3','u':'h_tr(C)',
             'v':'h_(tr(C))^2','sigma':'v-u^2','R':'h_tr(C^2)','T':'h_tr(K^2)','q':'Re h_tr(C K)','j':'h_(tr(D))^dagger tr(D)'},
        'centered_2x2':[[encoded(r) for r in row] for row in G],'flow_budget':encoded(flow),
        'exact_relations':['[[R-u^2-sigma,q-p*u],[q-p*u,T-p^2]] is PSD',
           '0 <= j <= a/2+q',
           'j <= a/2+p*u+sqrt((R-u^2-sigma)*(T-p^2))'],
        'finite_N_contact':'<f^2 tr DdagD>=<f^2 tr C K>+(1/2)<f^2 A>+(1/2)<f^2 tr Zbar tr Z>+<f^2 tr Zbar tr D>+(1/(2N^2))<(f^2)prime tr C>',
        'finite_N_contact_scope':'real reflection-invariant wavefunction and fixed radial polynomial. Center terms and derivative contact retained; only their proved weighted limit is removed.',
        'new_variables':0,'independence':'new ordinary trace-polynomial test blocks, no new stationarity rows; no claim of lowest degree or logical independence against every old inequality',
        'positive_fluctuation_reference':{'a':str(a),'e':str(e),'u':str(u),'v':str(v),'sigma':str(sigma),'j':str(j),'R':str(R),'T':str(K),'q':str(q),
           'centered_minor':str(det),'flow_margin':str(a/2+q-j),'full_system_feasibility':False}}

def free_algebra_centering():
    """Universal noncommuting coefficient identity, with ordering retained."""
    checks=0
    for N in range(1,6):
        for left,right in [('O','O'),('O','V'),('V','O'),('V','V')]:
            # Each letter records an operator entry and its adjoint bit.
            def word(letter,i,j,star):return (letter,i,j,star)
            def center(letter,i,j,star=False):
                p={word(letter,i,j,star):F(1)}
                if i==j:
                    for k in range(N):add(p,word(letter,k,k,star),F(-1,N))
                return p
            lhs={};rhs={}
            for i in range(N):
                for j in range(N):
                    for x,c in center(left,i,j,True).items():
                        for y,d in center(right,i,j).items():add(lhs,(x,y),c*d/N)
                    add(rhs,(word(left,i,j,True),word(right,i,j,False)),F(1,N))
            for i in range(N):
                for j in range(N):add(rhs,(word(left,i,i,True),word(right,j,j,False)),F(-1,N*N))
            if lhs!=rhs:raise AssertionError('noncommuting color identity')
            checks+=1
    return checks

def scalar_references():
    rng=Random(40040);checks=0
    for _ in range(100):
        s=F(rng.randint(1,20),17);t=F(rng.randint(1,20),19);r=F(rng.randint(-10,10),20)
        # [[s^2,s*t*r],[s*t*r,t^2]] PSD without a square root calculation.
        x=s*s;y=t*t;z=s*t*r
        assert x*y-z*z>=0;checks+=1
    return {'universal_noncommuting_color_identities':free_algebra_centering(),'coupling_positive_minor_references':checks,
       'tests_do_not_prove_all_N_analysis':True}
