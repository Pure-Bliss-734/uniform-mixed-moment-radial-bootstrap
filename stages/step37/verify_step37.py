"""Rebuild and check Stage 37. No predecessor code or numerical solver is imported.
The new all-N analytical statements are written proofs, not inferred from these tests.
"""
if not __debug__:
 raise RuntimeError('Run without -O: exact checks use assertions.')
import argparse,json,subprocess,sys,time
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import trace_algebra as a
ROOT=Path(__file__).resolve().parent
NAMES=['quotient11.json','cut_inventory.json','gram_blocks.json','covariance_structure.json',
       'energy_covariance_reduction.json','charged_quartic_bounds.json','exact_checks.json','summary.json']
def save(out,name,x):
 (out/name).write_text(json.dumps(x,ensure_ascii=False,sort_keys=True,indent=2)+'\n')
def enc(r):return [[int(i),str(c)] for i,c in sorted(r.items())]
def dec(r):return {int(i):F(c) for i,c in r}
def pack(m):
 return {'level':11,'charged_quartic_small_trace_lemma_used':True,
  'original_nonempty_neutral_words':len(a.words(11,0))-1,
  'monomials':[[list(w) for w in k] for k in m['keys']],
  'rank':len(m['E'].pivots),'free_columns':m['free'],
  'expressions':[enc(m['expressions'][i]) for i in range(len(m['keys']))],
  'raw_equation_counts':m['counts'],'nonzero_normalized_equations':len(m['equations']),
  'ordering':'single-trace orbits, double-trace orbits, then radial columns; whole-product reflection/reversal',
  'proof_scope':'Homogeneous radially weighted limiting functionals; not finite-N exact factorization.'}
def load(out):
 d=json.loads((out/NAMES[0]).read_text());keys=[tuple(tuple(w) for w in k) for k in d['monomials']]
 return {'L':11,'keys':keys,'index':{k:i for i,k in enumerate(keys)},'free':d['free_columns'],
  'expressions':{i:dec(r) for i,r in enumerate(d['expressions'])},'normalizer':a.TraceNormalizer(True)},d

def build_phase(out):
 m=a.build(11,True,True,verify_equations=False)
 assert len(m['keys'])==2701 and len(m['E'].pivots)==2551 and len(m['free'])==150
 assert len(m['equations'])==15304
 save(out,NAMES[0],pack(m))
 counts=Counter(x['reason'][0] for x in m['cuts'])
 save(out,NAMES[1],{'raw_ordered_nonempty_cuts':len(m['cuts']),'classified_counts':dict(counts),
  'all_cuts':m['cuts'],'free_double_trace_columns':[{'column':k,'factors':[list(w) for w in m['keys'][k]]}
  for k in m['free'] if len(m['keys'][k])==2],
  'new_uncertainties_not_factorized':['neutral quartic trace products','radius times six-weight traces','one-momentum cubic-current norm'],
  'zero_rules':'Only separately proved small traces and whole-moment symmetries; no neutral-trace covariance is silently deleted.'})
 print('BUILT: 10646 neutral words; 2701 symmetry-reduced columns; rank2551; free150.',flush=True)

def modular_rank(rows,prime=2147483647):
 piv={}
 for row in rows:
  r={i:(c.numerator*pow(c.denominator,-1,prime))%prime for i,c in row.items()}
  r={i:c for i,c in r.items() if c}
  while r:
   i=min(r);v=r[i]
   if i not in piv:
    iv=pow(v,-1,prime);piv[i]={j:c*iv%prime for j,c in r.items()};break
   for j,c in piv[i].items():
    z=(r.get(j,0)-v*c)%prime
    if z:r[j]=z
    elif j in r:del r[j]
 return len(piv)

def equations_phase(out):
 m,d=load(out);norm=m['normalizer'];rows=[];nnon=0;counts=Counter()
 for kind,w,mul,p in a.raw_system(11):
  pp=norm.polynomial(p,mul);counts[kind]+=1
  r={}
  for k,c in pp.items():a.combine(r,m['expressions'][m['index'][k]],c)
  assert not r,(kind,w,mul,r)
  if pp:
   rows.append({m['index'][k]:c for k,c in pp.items()});nnon+=1
 assert nnon==15304 and dict(counts)==d['raw_equation_counts']
 for i in m['free']:assert m['expressions'][i]=={i:F(1)}
 rank=modular_rank(rows);assert rank==d['rank']==2551
 # The explicit map with 150 independently free coordinates bounds rank from above;
 # a nonzero modular minor bounds it from below. Thus the rational rank is exactly2551.
 save(out,'rank_check.json',{'prime':2147483647,'rank_mod_prime':rank,
  'rational_parametrization_dimension':len(m['free']),'raw_equations_rechecked':sum(counts.values()),
  'nonzero_equations_rechecked':nnon,'rational_rank_certified_by_two_bounds':True})
 print('EQUATIONS: 18665 generated rows checked; independent modular rank2551.',flush=True)

def gram_phase(out):
 m,d=load(out);blocks=a.blocks(m,5)
 expected=[20,19,33,24,12,5,1]*2
 assert [len(b['basis']) for b in blocks]==expected
 data=[]
 for b in blocks:
  data.append({'name':b['name'],'dimension':len(b['basis']),
   'operator_basis':[[[list(w),str(c)] for w,c in sorted(p.items())] for p in b['basis']],
   'entries':[[enc(r) for r in row] for row in b['entries']]})
 save(out,NAMES[2],{'ordinary_and_ground_blocks':data,
  'dimension_order':expected,'symmetric_matrix_entries_checked':sum(n*n for n in expected),
  'warning':'These are symbolic PSD constraints, not numerically feasible matrices.'})
 A=(0,1);aa=(0,0,1,1);bb=(0,1,0,1);pp=(2,3)
 pols=[{():F(1)},{A:F(1)},{pp:F(-1)},{aa:F(1)},{bb:F(1)}]
 G=[]
 for u in pols:
  row=[]
  for v in pols:
   p={}
   for w,c in u.items():
    sg,dw=a.involution(w)
    for t,e in v.items():a.add(p,a.raw_monomial(dw,t),sg*c*e)
   row.append(a.product_expression(p,m))
  G.append(row)
 kernel=[0,-1,1,2,-2]
 for row in G:
  r={}
  for c,p in zip(kernel,row):a.combine(r,p,F(c))
  assert not r
 assert all(G[i][j]==G[j][i] for i in range(5) for j in range(5))
 J=(0,0,3);sg,dJ=a.involution(J);jvar=a.product_expression({(dJ,J):F(sg)},m)
 assert jvar
 save(out,NAMES[3],{'observable_order':['1','A','p','a4','b4'],
  'Gram5':[[enc(r) for r in row] for row in G],
  'exact_virial_null_vector':kernel,'four_by_four_core_indices':[0,1,3,4],
  'current_J':list(J),'current_variance_expression':enc(jvar),'current_variance_not_set_to_zero':True,
  'interpretation':'G5 has the exact limiting kernel p-A+2a4-2b4=0. Core moments include E[a4^2], E[a4 b4], E[b4^2], not products of their means.'})
 print('GRAM: 14 order-five matrix blocks and unfactorized scalar covariance checked.',flush=True)

def tests_phase(out):
 import exact_checks as t
 m,d=load(out);r=t.run(m);save(out,NAMES[6],r)
 save(out,NAMES[4],t.covariance_checks())
 save(out,NAMES[5],{
 'U4':'tr Z^4','V4':'tr Z^3 Zbar',
 'bounds':{'E_abs_U4_squared':'15/(2 M^4 N^2)','E_abs_V4_squared':'15/(8 M^4 N^2)'},
 'imaginary_identities':{'Im_U4':'tr X^3 Y - tr X Y^3','Im_V4':'(tr X^3 Y + tr X Y^3)/2'},
 'conditional_Poincare_input':'E |tr X^3 Y|^2, E |tr X Y^3|^2 <= 15/(16 M^4 N^2)',
 'O2_input':'E U4^2=E V4^2=0 by charges4 and2; E|T|^2=2E(Im T)^2',
 'limitations':'Written implication from inherited physical conditional Gaussian comparison; finite-N variances are not zero; neutral fourth traces not covered.',
 'finite_free_reference_polynomials':r['Gaussian_quartic_variance_Laurent']})
 print('REFERENCES:',r['counts'],'negative controls',len(r['negative_controls_rejected']),flush=True)

def summary_phase(out):
 d=json.loads((out/NAMES[0]).read_text());r=json.loads((out/NAMES[6]).read_text());rank=json.loads((out/'rank_check.json').read_text())
 save(out,NAMES[-1],{
        'finite_reference_counts': r['counts'],
        'free_coordinates': 150,
        'independent_modular_rank_check': rank,
        'negative_controls': len(r['negative_controls_rejected']),
        'original_neutral_word_count': 10646,
        'rational_rank': 2551,
        'symmetry_reduced_columns': 2701,
    })
 print('All Stage-37 exact finite checks passed. No new radius bound is claimed.',flush=True)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,default=ROOT/'results/verified')
 ap.add_argument('--phase',choices=['build','equations','gram','tests','summary']);args=ap.parse_args();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
 phases={'build':build_phase,'equations':equations_phase,'gram':gram_phase,'tests':tests_phase,'summary':summary_phase}
 if args.phase:phases[args.phase](out)
 else:
  # Separate processes avoid retaining the large intermediate elimination graph.
  for name in phases:
   p=subprocess.run([sys.executable,'-B','-S',str(Path(__file__).resolve()),'--phase',name,'--output-dir',str(out)],cwd=ROOT)
   if p.returncode:raise SystemExit(p.returncode)
if __name__=='__main__':main()
