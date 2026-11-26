"""Exact modular audit: can each chosen support work after centered atoms are removed?
This is a support-specific check, NOT a logical-independence proof for whole cones.
"""
from fractions import Fraction as F
from pathlib import Path
import json
P=2147483647

def rank_mod(rows):
    piv={}
    for row in rows:
        r=[int(x)%P for x in row]
        for j in range(len(r)):
            if not r[j]:continue
            if j in piv:
                c=r[j];b=piv[j]
                for k in range(j,len(r)):r[k]=(r[k]-c*b[k])%P
            else:
                c=pow(r[j],-1,P);piv[j]=[(x*c)%P for x in r];break
    return len(piv)
def audit(cert):
    c=json.loads(Path(cert).read_text());x=F(c['radius']);expr=[{(k,a,e):F(v) for k,a,e,v in r} for r in c['expressions']]
    ids=[i for i,at in enumerate(c['atoms']) if at['block']<16];removed=[i for i,at in enumerate(c['atoms']) if at['block']>=16]
    keys=sorted({(k,e) for r in expr for k,a,e in r if k});mat=[]
    for k,e in keys:
        row=[]
        for i in ids:
            v=sum(v*x**a for (kk,a,ee),v in expr[i].items() if kk==k and ee==e)
            if v.denominator%P==0:raise ValueError('bad modular prime')
            row.append(v.numerator*pow(v.denominator,-1,P)%P)
        mat.append(row)
    rk=rank_mod(mat)
    return {'radius':str(x),'prime':P,'remaining_original_atoms':len(ids),'removed_centered_atoms':len(removed),'rank_mod_prime':rk,
       'full_column_rank':rk==len(ids),
       'claim':'At this fixed rational radius, full column rank certifies that the chosen remaining old atoms have no nonzero rational homogeneous-stationarity weights. It does not exclude other old supports or prove cone-level independence.'}
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('cert',type=Path);a=p.parse_args();print(json.dumps(audit(a.cert),indent=2))
