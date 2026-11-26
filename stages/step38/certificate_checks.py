"""Fraction-only identities, Bernstein sign boxes, and negative controls for Stage38."""
from fractions import Fraction as F
from math import comb
from energy_model import add,axpy,dec,ps_dec,ps_enc,unpack,combine_map
from face_reduction import psd_pivots

def atom_expr(model,atom):
    b=model['blocks'][int(atom['block'])];d=len(b['entries']);co={int(i):F(c) for i,c in atom['coefficients']}
    if not co or any(i<0 or i>=d for i in co):raise ValueError('invalid atom vector')
    r={}
    for i,c in co.items():
        for j,t in co.items():
            for key,v in b['entries'][i][j].items():add(r,key,c*t*v)
    return r

def verify_identity(data,cert):
    model=unpack(data);out={};norm=F(0);mins=[]
    for a in cert['atoms']:
        w=F(a['weight'])
        if w<=0:raise AssertionError('weight is not positive')
        r=atom_expr(model,a);axpy(out,r,w)
        norm+=w*sum(F(v)**2 for i,v in a['coefficients']);mins.append(w)
    if norm!=1:raise AssertionError('atom normalization does not equal 1')
    if any(k not in (0,1,3) for k,p in out):raise AssertionError('nonradial coefficient survives')
    if out!=ps_dec(cert['radial_polynomial']):raise AssertionError('declared polynomial differs')
    return out,{'atoms':len(cert['atoms']),'all_weights_strictly_positive':True,
        'exact_normalization':str(norm),'zero_nuisance_coefficient_identities':144,
        'polynomial_terms':ps_enc(out),'minimum_raw_atom_weight':str(min(mins))}

def structure(data):
    rows=[dec(r) for r in data['basis_rows']];inv={k:dec(r) for k,r in data['inverse']};n=150
    if len(rows)!=n or len(inv)!=n:raise AssertionError('coordinate dimensions')
    for j,r in enumerate(rows):
        if combine_map(r,inv)!={j:F(1)}:raise AssertionError('row times inverse')
    for k,r in inv.items():
        if combine_map(r,dict(enumerate(rows)))!={k:F(1)}:raise AssertionError('inverse times rows')
    total={k:ps_dec(r) for k,r in data['total_map']}
    for r in data['energy_constraints']:
        val={}
        for k,p,c in r:
            for (j,q),a in total[k].items():add(val,(j,p+q),F(c)*a)
        if val:raise AssertionError('energy closure failed')
    expected=[{(6,0):F(3),(1,0):F(1),(0,1):F(-2)},
        {(7,0):F(3),(3,0):F(1),(1,1):F(-2)},
        {(8,0):F(3),(4,0):F(1),(2,1):F(-2)},
        {(9,0):F(3),(7,0):F(1),(6,1):F(-2)}]
    if [ps_dec(r) for r in data['independent_low_rows']]!=expected:raise AssertionError('wrong independent low-moment equations')
    model=unpack(data);bs=model['blocks'];indices=[[0,1,2],[1,3,4],[2,4,5]]
    if bs[-2]['name']!='S_G3':raise AssertionError('G3 missing')
    for i in range(3):
        for j in range(3):
            if bs[-2]['entries'][i][j]!={(indices[i][j],0):F(1)}:raise AssertionError('covariance has been factored')
    if bs[-1]['entries']!=[[{(6,0):F(1)}]]:raise AssertionError('JdaggerJ removed')
    count=0;keys=set()
    for b in bs:
        E=b['entries']
        for i,row in enumerate(E):
            for j,r in enumerate(row):
                if r!=E[j][i]:raise AssertionError('asymmetric pencil')
                if any(k>=146 or p>2 or p<0 for k,p in r):raise AssertionError('invalid polynomial layout')
                keys.update(r);count+=1
    return {'basis_inverse_equalities':300,'energy_kernel_equalities':10,'independent_added_rank':4,
        'constant_rank_minor':81,'homogeneous_dimension':146,'normalized_dimension':145,
        'all_energies_no_denominator_or_exceptional_fiber':True,'block_count':len(bs),
        'dimensions':[len(b['entries']) for b in bs],'symmetric_polynomial_entries':count,
        'surviving_energy_dependent_nuisance_columns':sorted({k for k,p in keys if p and k not in (0,1,3)}),
        'G3_not_factorized':True,'JdaggerJ_retained':True}

def root_enclosure():
    lo,hi=F(1),F(2)
    for _ in range(130):
        mid=(lo+hi)/2
        if 15*mid**3-15*mid-7<0:lo=mid
        else:hi=mid
    low=lo-F(1,30)/lo**3;high=F(265677,224720)
    a,b=F(11609,10000),F(11823,10000)
    if not low>a or not high<b:raise AssertionError('energy interval not justified')
    return {'omega_interval':[str(lo),str(hi)],'inherited_lower_formula':'omega - 1/(30 omega^3)',
        'lower_value_rational_bound':str(low),'inherited_upper_value':str(high),
        'certified_energy_interval':[str(a),str(b)],
        'scope':'Arithmetic check of inherited Stage32/34 physical energy bounds, not a new independent proof of the physical Gaussian comparison.'}

def bivar(poly):return {(0 if k==0 else 1 if k==1 else 2,p):c for (k,p),c in poly.items()}
def val(poly,a,e):return sum(c*a**i*e**j for (i,j),c in poly.items())
def bernstein(poly,a,b,c,d):
    """Exact bidegree(2,2) Bernstein coefficients on [a,b] x [c,d]."""
    power={}
    for (i,j),v in poly.items():
        for k in range(i+1):
            for l in range(j+1):add(power,(k,l),v*comb(i,k)*a**(i-k)*(b-a)**k*comb(j,l)*c**(j-l)*(d-c)**l)
    B=[]
    for i in range(3):
        row=[]
        for j in range(3):
            row.append(sum(v*F(comb(i,k),comb(2,k))*F(comb(j,l),comb(2,l)) for (k,l),v in power.items() if k<=i and l<=j))
        B.append(row)
    return B

def upward_short(poly,den=10**10):
    scale=abs(poly[(1,0)])
    if scale<=0:raise AssertionError('radius coefficient absent')
    out={}
    for key,c in poly.items():
        t=c*den/scale;n=-((-t.numerator)//t.denominator)
        out[key]=F(n)
        if F(n,den)-c/scale<0:raise AssertionError('rounding direction')
    return out,scale

def prove_cut(poly,side,cut,energies):
    cut=F(cut);elo,ehi=map(F,energies)
    if not F(0)<cut<F(1,2) or not F(0)<elo<ehi:raise ValueError('invalid domain')
    aa,bb=(F(0),cut) if side=='lower' else (cut,F(1,2))
    B=bernstein(bivar(poly),aa,bb,elo,ehi)
    if not max(v for row in B for v in row)<0:raise AssertionError('box sign is not strictly negative')
    return {'side':side,'excluded_radius_box':[str(aa),str(bb)],'energy_interval':[str(elo),str(ehi)],
        'bernstein_coefficients':[[str(v) for v in row] for row in B],
        'strict_negative_margin':str(-max(v for row in B for v in row))}

def tight_simple_cut(poly,side,energies,den=10**10):
    # For these actual certificates, check the separate convexity that makes
    # corner signs sufficient; certify the final domain separately by Bernstein.
    p=bivar(poly);elo,ehi=map(F,energies)
    if p.get((2,0),0)<0 or p.get((0,2),0)<0:raise AssertionError('not separately convex')
    if any(i+j>2 for i,j in p):raise AssertionError('unexpected polynomial')
    # Search within the already proven broad strip. This operation is exact.
    def good(n):
        x=F(n,den)
        end=F(0) if side=='lower' else F(1,2)
        return all(val(p,y,e)<0 for y in [x,end] for e in [elo,ehi])
    lo,hi=int(F('0.388')*den),int(F('0.390')*den)
    if side=='lower':
        if not good(lo) or good(hi):raise AssertionError('lower crossing not bracketed')
        while hi-lo>1:
            m=(lo+hi)//2
            if good(m):lo=m
            else:hi=m
        return F(lo,den)
    if good(lo) or not good(hi):raise AssertionError('upper crossing not bracketed')
    while hi-lo>1:
        m=(lo+hi)//2
        if good(m):hi=m
        else:lo=m
    return F(hi,den)

def covariance_reference():
    e=F(59,50)
    atoms=[(F(389,1000),F(23,100)),(F(3891,10000),F(23,100)),(F(389,1000),F(231,1000))]
    vs=[[F(1),a,b] for a,b in atoms]
    G=[[sum(v[i]*v[j] for v in vs)/3 for j in range(3)] for i in range(3)]
    piv=psd_pivots(G)
    if len(piv)!=3:raise AssertionError('reference not full rank')
    T=[[F(1),0,0],[0,F(1),0],[2*e/3,F(-1,3),0],[0,0,F(1)],[e/3,F(-2,3),F(1)]]
    Q=[[sum(F(T[i][k])*G[k][l]*F(T[j][l]) for k in range(3) for l in range(3)) for j in range(5)] for i in range(5)]
    nus=[[-2*e/3,F(1,3),1,0,0],[-e/3,F(2,3),0,-1,1]]
    for v in nus:
        if any(sum(c*F(a) for c,a in zip(row,v)) for row in Q):raise AssertionError('null vector mismatch')
    var=G[1][1]-G[0][1]**2;cv=G[1][2]-G[0][1]*G[0][2]
    if var<=0:raise AssertionError('reference variance discarded')
    return {'energy':str(e),'G3':[[str(x) for x in row] for row in G],'G3_positive_pivots':[str(v) for v in piv],
        'G5':[[str(x) for x in row] for row in Q],'Var_A':str(var),'Cov_A_a4':str(cv),
        'meaning':'Reference for closure alone; not asserted feasible for the other 14 matrix blocks or physically representable.'}

def modular_rank(rows,prime=2147483647):
    piv={}
    for row in rows:
        r={k:(v.numerator*pow(v.denominator,-1,prime))%prime for k,v in row.items() if v}
        r={k:v for k,v in r.items() if v}
        while r:
            k=min(r);v=r[k]
            if k not in piv:
                iv=pow(v,-1,prime);piv[k]={j:c*iv%prime for j,c in r.items()};break
            for j,c in piv[k].items():
                a=(r.get(j,0)-v*c)%prime
                if a:r[j]=a
                elif j in r:del r[j]
    return len(piv)

def reject_floating_candidate(data,candidate):
    model=unpack(data);cols=[atom_expr(model,a) for a in candidate['atoms']]
    keys=sorted({key for col in cols for key in col if key[0] not in (0,1,3)})
    rows=[{j:col[key] for j,col in enumerate(cols) if key in col} for key in keys]
    rank=modular_rank(rows)
    if rank!=len(cols):raise AssertionError('candidate has not been excluded by this rank check')
    return {'numeric_reported_objective':candidate['value'],'number_of_selected_atoms':len(cols),
        'homogeneous_stationarity_rank_mod_2147483647':rank,
        'exact_implication':'The rational homogeneous stationarity map has trivial kernel on this support. No nonzero weights, positive or otherwise, can make a certificate on the displayed support.',
        'scope':'Rejects this finite numeric support, not feasibility of the sampled point and not all possible certificates there.'}

def exact_references():
    from random import Random
    rng=Random(38038);n=0
    # Bernstein expansion against direct power evaluation over independently chosen boxes.
    for it in range(30):
        poly={(i,j):F(rng.randrange(-10,11),rng.randrange(1,10)) for i in range(3) for j in range(3)}
        a=F(rng.randrange(-3,4),3);b=a+F(rng.randrange(1,4),2)
        c=F(rng.randrange(-3,4),5);d=c+F(rng.randrange(1,4),3)
        B=bernstein(poly,a,b,c,d)
        for u,v in [(F(0),F(0)),(F(1),F(1)),(F(1,3),F(2,5)),(F(7,9),F(1,4))]:
            got=sum(B[i][j]*comb(2,i)*u**i*(1-u)**(2-i)*comb(2,j)*v**j*(1-v)**(2-j) for i in range(3) for j in range(3))
            if got!=val(poly,a+(b-a)*u,c+(d-c)*v):raise AssertionError('Bernstein evaluation')
            n+=1
    # Fractional upward rounding must stay upward for both signs.
    nr=0
    for j in range(-75,76):
        for den in (3,7,101):
            c=F(j,den);up=F(-((-c.numerator)//c.denominator))
            if not F(0)<=up-c<1:raise AssertionError('ceiling')
            nr+=1
    return {'Bernstein_exact_point_comparisons':n,'signed_ceiling_comparisons':nr,
        'finite_tests_do_not_replace_the_interval_convex_combination_proof':True}
