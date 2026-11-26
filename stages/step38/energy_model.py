"""Stage 38: low-observable coordinates and an exact Q[e] energy-conditioned pencil.
No predecessor Python modules are imported; only pinned JSON model data are read.
The predecessor's physical meaning is an explicit inherited premise.
"""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib

def add(d,k,c):
    if c:
        d[k]=d.get(k,F(0))+c
        if not d[k]:del d[k]
def axpy(d,r,c=F(1)):
    for k,v in r.items():add(d,k,c*v)
    return d

def dec(r):return {int(k):F(v) for k,v in r}
def enc(r):return [[int(k),str(v)] for k,v in sorted(r.items())]
def ps_enc(r):return [[int(k),int(p),str(v)] for (k,p),v in sorted(r.items())]
def ps_dec(r):return {(int(k),int(p)):F(v) for k,p,v in r}
def combine_map(v,maps):
    out={}
    for k,c in v.items():axpy(out,maps[k],c)
    return out

class Basis:
    """Select independent observable rows; record exact row operations."""
    def __init__(self):self.piv={};self.rows=[];self.names=[];self.pivot_values=[]
    def offer(self,v,name):
        j=len(self.rows);r=dict(v);rhs={j:F(1)}
        while r:
            p=min(r)
            if p not in self.piv:
                a=r[p]; self.pivot_values.append(a)
                self.piv[p]=({k:c/a for k,c in r.items()},{k:c/a for k,c in rhs.items()})
                self.rows.append(dict(v));self.names.append(name);return True
            row,b=self.piv[p];a=r[p];axpy(r,row,-a);axpy(rhs,b,-a)
        return False
    def inverse(self):
        out={}
        for p in sorted(self.piv,reverse=True):
            row,rhs=self.piv[p];r=dict(rhs)
            for k,c in row.items():
                if k!=p:axpy(r,out[k],-c)
            out[p]=r
        return out

def read_source(directory):
    d=Path(directory)
    paths={x:d/x for x in ['quotient11.json','gram_blocks.json','covariance_structure.json']}
    raw={k:json.loads(p.read_text()) for k,p in paths.items()}
    hashes={k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in paths.items()}
    return raw,hashes

def construct(directory,progress=False):
    raw,hashes=read_source(directory);q=raw['quotient11.json'];cov=raw['covariance_structure.json']
    oldfree=q['free_columns'];expr=[dec(r) for r in q['expressions']]
    if len(oldfree)!=150:raise ValueError('Expected the pinned 150-coordinate Stage37 model')
    G=[[dec(r) for r in row] for row in cov['Gram5']]
    names=['m0','m1','d1','m2','dA','d2','p','Ap','pa4','pp','JdaggerJ']
    special=[G[0][0],G[0][1],G[0][3],G[1][1],G[1][3],G[3][3],G[0][2],G[1][2],G[2][3],G[2][2],dec(cov['current_variance_expression'])]
    B=Basis()
    for v,name in zip(special,names):
        if not B.offer(v,name):raise AssertionError(('dependent special coordinate',name))
    def key(i):
        m=q['monomials'][i]
        return (sum((1 if x<2 else 2) for w in m for x in w),len(m),m)
    chosen=[]
    for i in sorted(range(len(expr)),key=key):
        if B.offer(expr[i],{'source_column':i,'trace_factors':q['monomials'][i]}):chosen.append(i)
        if len(B.rows)==150:break
    if len(B.rows)!=150:raise AssertionError('low-observable basis not full')
    inv=B.inverse()
    # Bijective coordinate change, tested both ways rather than by a floating determinant.
    for i,r in enumerate(B.rows):
        if combine_map(r,inv)!={i:F(1)}:raise AssertionError('basis/inverse identity')
    for i,r in inv.items():
        out={}
        for j,c in r.items():axpy(out,B.rows[j],c)
        if out!={i:F(1)}:raise AssertionError('inverse/basis identity')
    # Ordered new six scalar moments then J norm etc. Remove p,Ap,pa4,pp only.
    kept=[i for i in range(150) if i not in (6,7,8,9)];ki={j:i for i,j in enumerate(kept)}
    energy={i:{(ki[i],0):F(1)} for i in kept}
    def S(items):return {(ki[k],p):F(v) for k,p,v in items}
    energy[6]=S([(0,1,F(2,3)),(1,0,F(-1,3))])
    energy[7]=S([(1,1,F(2,3)),(3,0,F(-1,3))])
    energy[8]=S([(2,1,F(2,3)),(4,0,F(-1,3))])
    energy[9]=S([(0,2,F(4,9)),(1,1,F(-4,9)),(3,0,F(1,9))])
    total={k:combine_map(v,energy) for k,v in inv.items()}
    def transform(r):return combine_map(r,total)
    Gnew=[[transform(r) for r in row] for row in G]
    # Check all ten energy null-vector equations as identities in e.
    constraints=[]
    for coeff in [[(0,1,F(-2)),(1,0,F(1)),(2,0,F(3))],[(0,1,F(-1)),(1,0,F(2)),(3,0,F(-3)),(4,0,F(3))]]:
        for j in range(5):
            unred={};red={}
            for l,p,c in coeff:
                for k,v in G[l][j].items():add(unred,(k,p),c*v)
                for (k,p2),v in Gnew[l][j].items():add(red,(k,p+p2),c*v)
            if red:raise AssertionError(('energy kernel',red))
            constraints.append(unred)
    # Before energy elimination, the four equations determine those four coordinates
    # with constant pivots 3,3,3,3. pp uses the first two identities.
    indep=[]
    for j in [0,1,3,2]:
        rr={}
        for l,p,c in [(0,1,F(-2)),(1,0,F(1)),(2,0,F(3))]:
            v=combine_map(G[l][j],inv)
            for k,t in v.items():add(rr,(k,p),c*t)
        indep.append(rr)
    def multiply_e(r,p=1):return {(k,n+p):v for (k,n),v in r.items()}
    # This triangular 4x4 minor has constant nonzero determinant81, for every real e.
    targets=[6,7,8,9]
    for i,r in enumerate(indep):
        for j in range(i+1,4):
            if any(k==targets[j] for (k,p) in r):raise AssertionError('not triangular')
        if r.get((targets[i],0))!=3:raise AssertionError('pivot normalization')
    bs=[]
    for b in raw['gram_blocks.json']['ordinary_and_ground_blocks']:
        ent=[[transform(dec(r)) for r in row] for row in b['entries']]
        if any(ent[i][j]!=ent[j][i] for i in range(len(ent)) for j in range(len(ent))):raise AssertionError('not symmetric')
        bs.append({'name':b['name'],'basis':b['operator_basis'],'entries':ent})
    # Retain the full scalar Gram3. Do not center by setting variances to zero.
    ids=[0,1,3]
    bs.append({'name':'S_G3','basis':[[[[], '1']],[[[0,1],'1']],[[[0,0,1,1],'1']]],'entries':[[Gnew[i][j] for j in ids] for i in ids]})
    jexpr=transform(special[10]);bs.append({'name':'S_JdaggerJ','basis':[[[[0,0,3],'1']]],'entries':[[jexpr]]})
    if jexpr!={(6,0):F(1)}:raise AssertionError('J norm was lost')
    if progress:print('energy map ready: 150 ->146, 10 kernel equations, 4 independent; Gram3 and J retained',flush=True)
    return {'source_hashes':hashes,'source_dimensions':{'columns':len(q['monomials']),'rank':q['rank'],'free':150},
        'basis_names':B.names,'basis_rows':B.rows,'basis_pivots':B.pivot_values,'inverse':inv,'energy_map':energy,
        'total_map':total,'old_free':oldfree,'kept':kept,'Gnew':Gnew,'energy_rows':constraints,'independent_rows':indep,
        'blocks':bs,'source_raw':raw,'new_dimension':len(kept)}

def evaluated_matrix(block,e,n=146):
    out=[[[F(0) for _ in range(len(block['entries']))] for _ in block['entries']] for _ in range(n)]
    for i,row in enumerate(block['entries']):
        for j,r in enumerate(row):
            for (k,p),c in r.items():out[k][i][j]+=c*e**p
    return out

def pack(m):
    return {'stage':38,'energy_parameter':'e (fixed subsequential total energy density, not a radius-branch random variable)',
    'source_hashes':m['source_hashes'],'old_dimension':150,'dimension':m['new_dimension'],
    'low_observable_order':['m0','m1','d1','m2','dA','d2','JdaggerJ'],
    'basis_names':m['basis_names'],'kept_basis_indices':m['kept'],
    'basis_rows':[enc(r) for r in m['basis_rows']],
    'basis_pivots':[str(x) for x in m['basis_pivots']],
    'inverse':[[int(k),enc(r)] for k,r in sorted(m['inverse'].items())],
    'total_map':[[int(k),ps_enc(r)] for k,r in sorted(m['total_map'].items())],
    'energy_constraints':[ps_enc(r) for r in m['energy_rows']],
    'independent_low_rows':[ps_enc(r) for r in m['independent_rows']],
    'all_energies_constant_rank_minor':'81','independent_new_rank':4,
    'blocks':[{'name':b['name'],'basis':b['basis'],'dimension':len(b['entries']),
        'entries':[[ps_enc(r) for r in row] for row in b['entries']]} for b in m['blocks']]}

def unpack(d):
    return {'new_dimension':d['dimension'],'blocks':[{'name':b['name'],'basis':b['basis'],
        'entries':[[ps_dec(r) for r in row] for row in b['entries']]} for b in d['blocks']],
        'total_map':{k:ps_dec(r) for k,r in d['total_map']},'source_hashes':d['source_hashes']}

if __name__=='__main__':
    import argparse,time
    ap=argparse.ArgumentParser();ap.add_argument('--source-dir',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();t=__import__('time').time();m=construct(args.source_dir,True)
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(pack(m),sort_keys=True,indent=2)+'\n')
    print('Reconstructed the energy pencil.')
