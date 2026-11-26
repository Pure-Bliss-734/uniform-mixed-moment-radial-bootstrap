"""Weight<=4 trace images derived from finite-N contacts/gauge/time derivatives.
Independent of the old word normalizer or its sparse elimination. These are
limiting images, not asserted as exact finite-N operator identities.
"""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
W=(1,1,2,2);CH=(1,-1,1,-1)
def add(d,k,v):
    if v:
        d[k]=d.get(k,F(0))+v
        if not d[k]:del d[k]
def cyc(w):return min(w[i:]+w[:i] for i in range(len(w))) if w else ()
def pmul(p,q):
    r={}
    for (a,e),c in p.items():
        for (b,f),d in q.items():add(r,(a+b,e+f),c*d)
    return r

def word_image(w):
    """Return (charge, constant polynomial(a,e), a4 coefficient, J coefficient).
    charge -1 uses J^dagger as the reference. Other nonzero charges map to zero.
    """
    w=tuple(w);degree=sum(W[x] for x in w);q=sum(CH[x] for x in w)
    if degree>4:raise ValueError('not a proved weight<=4 image')
    if not w:return 0,{(0,0):F(1)},F(0),F(0)
    if len(w)==1:return q,{},F(0),F(0)
    np=sum(x>=2 for x in w)
    if np==0:
        if len(w)%2 or q:return q,{},F(0),F(0)
        if len(w)==2:return 0,{(1,0):F(1)},F(0),F(0)
        if len(w)==4:
            c=cyc(w)
            if c==(0,0,1,1):return 0,{},F(1),F(0)
            if c==(0,1,0,1):return 0,{(0,1):F(1,3),(1,0):F(-2,3)},F(1),F(0)
        raise ValueError(('unexpected coordinate image',w))
    if len(w)==2 and np==1:
        return (q,{(0,0):F(1,2) if w[0]<2 else F(-1,2)},F(0),F(0)) if not q else (q,{},F(0),F(0))
    if len(w)==2 and np==2:
        return (q,{(0,1):F(-2,3),(1,0):F(1,3)},F(0),F(0)) if not q else (q,{},F(0),F(0))
    if len(w)==3 and np==1:
        if abs(q)==3:return q,{},F(0),F(0)
        if q==1:return q,{},F(0), F(1) if 3 in w else F(-1,2)
        if q==-1:return q,{},F(0), F(-1) if 2 in w else F(1,2)
    raise ValueError(('unclassified image',w))

def image(poly):
    out=None
    for w,c in poly:
        q,p,k,r=word_image(w);c=F(c)
        if out is None:out=[q,{},F(0),F(0)]
        if q!=out[0]:raise ValueError('basis mixes charges')
        for key,v in p.items():add(out[1],key,c*v)
        out[2]+=c*k;out[3]+=c*r
    return out or [0,{},F(0),F(0)]

def scalar_gram(x,y):
    if x[0]!=y[0]:raise ValueError('different charge block')
    q=x[0];r={}
    if q==0:
        for (a,e),c in pmul(x[1],y[1]).items():add(r,(0,a,e),c)
        for key,c in x[1].items():add(r,(1,*key),c*y[2])
        for key,c in y[1].items():add(r,(1,*key),c*x[2])
        add(r,(2,0,0),x[2]*y[2])
    elif abs(q)==1:add(r,(3,0,0),x[3]*y[3])
    elif x[1] or y[1] or any([x[2],x[3],y[2],y[3]]):raise AssertionError('charged field not proved small')
    return r

def encoded(r):return [[k,a,e,str(c)] for (k,a,e),c in sorted(r.items())]
def rebuild(pencil):
    data=json.loads(Path(pencil).read_text());out=[];images=[]
    for ib,b in enumerate(data['blocks']):
        if not b['name'].startswith('M_'):continue
        ids=[i for i,p in enumerate(b['basis']) if max((sum(W[x] for x in w) for w,c in p),default=0)<=4]
        if not ids:continue
        im=[image(b['basis'][i]) for i in ids];S=[[scalar_gram(x,y) for y in im] for x in im];C=[]
        for ii,i in enumerate(ids):
            row=[]
            for jj,j in enumerate(ids):
                r={(k,a,e):F(c) for k,a,e,c in b['entries'][i][j]}
                for key,v in S[ii][jj].items():add(r,key,-v)
                row.append(encoded(r))
            C.append(row)
        out.append({'name':'C_'+b['name'][2:],'source_block':ib,'source_indices':ids,'basis':[b['basis'][i] for i in ids],
          'entries':C,'subtracted_scalar_gram':[[encoded(r) for r in row] for row in S]})
        images.append({'block':b['name'],'images':[{'charge':v[0],'constant_polynomial':[[a,e,str(c)] for (a,e),c in sorted(v[1].items())],
          'neutral_a4_coefficient':str(v[2]),'current_coefficient':str(v[3])} for v in im]})
    return {'stage':40,'pencil_sha256':hashlib.sha256(Path(pencil).read_bytes()).hexdigest(),'new_variables':0,
            'blocks':out,'trace_images':images,'trace_count':sum(len(x['images']) for x in images),
            'entry_count':sum(len(x['entries'])**2 for x in out)}

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--pencil',required=True,type=Path);p.add_argument('--compare-quotient-build',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=rebuild(a.pencil)
    if a.compare_quotient_build:
        old=json.loads(a.compare_quotient_build.read_text())
        for x,y in zip(old['blocks'],r['blocks']):
            for k in ('name','source_indices','basis','entries','subtracted_scalar_gram'):
                if x[k]!=y[k]:raise AssertionError(('independent trace-image discrepancy',x['name'],k))
        r['all_entries_match_independent_quotient_build']=True
    a.output.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n');print('Independent trace images:',r['trace_count'],'traces,',r['entry_count'],'matrix elements; exact agreement.')
