"""Collect our four real u-sector primitives in their algebraic coefficient field.
GPL words keep their physical path and r=sqrt(s) basepoint. This is algebraic
collection, not a replacement of any sector by the published total.
"""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
from collections import defaultdict
import json,time
import sympy as S
B=Path(__file__).resolve().parent
r,y,ell,z,s,H=S.symbols('r y ell z s H');G=S.Function('G')
files={'R0':'compton_real_analytic_rational_gpl.json','W':'compton_real_analytic_single_gpl.json','logs':'compton_real_analytic_logs_gpl.json','D':'compton_box_analytic_bose_gpl.json'}
relation=S.Poly(z*z-(1-r*r),z)
def reduce_root(p):
 return S.Poly(S.expand(p),z).rem(relation).as_expr()
def powers_to_s(p):
 # Coefficients of words beginning at +r and -r need not be even.
 # Keep their odd factor r, with r=sqrt(s), until these words are combined.
 poly=S.Poly(p,r)
 value=sum(c*r**(k[0]%2)*s**(k[0]//2) for k,c in poly.terms())
 assert S.expand(value.subs(s,r*r)-p)==0
 return value
def rationalize(e):
 e=e.subs(H,r*r-1)
 root_powers={p:z**(2*p.exp) for p in e.atoms(S.Pow)
              if p.exp.is_Rational and p.exp.q==2 and S.expand(p.base)==1-r*r}
 e=e.xreplace(root_powers)
 n,d=S.fraction(S.together(e))
 n,d=reduce_root(n),reduce_root(d)
 if d.has(z):
  conjugate=d.subs(z,-z)
  n,d=reduce_root(n*conjugate),reduce_root(d*conjugate)
 assert not d.has(z)
 # Odd factors can occur in numerator and denominator even when their ratio
 # is invariant under r -> -r. Cancel them before changing to s=r^2.
 n,d=S.fraction(S.cancel(n/d))
 value=S.cancel(powers_to_s(n)/powers_to_s(d))
 residual=S.fraction(S.together(value.subs(s,r*r)-e))[0]
 assert S.expand(reduce_root(residual))==0
 return value
co=defaultdict(lambda:S.S.Zero);counts={};nonzero_counts={};sources={}
for sector,name in files.items():
 start=time.time();data=json.loads(read_artifact(B/name))['coefficients'];counts[sector]=len(data)-1
 sources[sector]={}
 nonzero_counts[sector]=sum(key!='constant' and S.sympify(value)!=0 for key,value in data.items())
 for key,c in data.items():
  w=() if key=='constant' else tuple(S.expand(a.subs(H,r*r-1)).subs(r*r,s) for a in S.sympify(key,locals={'G':G}).args[0])
  value=rationalize(S.sympify(c))*(ell if sector=='logs' else 1)
  sources[sector][str(w)]=str(value);co[w]+=value
 print(sector,'algebraic coefficients',len(data),'seconds',round(time.time()-start,2),flush=True)
def clean(e):
 p=S.Poly(e,ell)
 return sum(S.factor(S.cancel(c))*ell**k[0] for k,c in p.terms())
out={}
for i,(w,c) in enumerate(co.items()):
 c=clean(c)
 if c!=0:out[w]=c
cancelled=[str(w) for w,c in co.items() if w and clean(c)==0 and any(S.sympify(d.get(str(w),'0'))!=0 for d in sources.values())]
print('collected',len(out)-1,'cancelled nonzero words',len(cancelled),flush=True)
before=dict(out);products={}
for w in list(out):
 if w not in out or len(w)!=2:continue
 a,b=w
 if a==b:products[w]=out.pop(w)/2
 elif str(a)>str(b) and (b,a) in out:
  c=out.pop(w);out[(b,a)]=clean(out[(b,a)]-c);products[(b,a)]=c
out={w:c for w,c in out.items() if c!=0}
res=defaultdict(lambda:S.S.Zero,out)
for (a,b),c in products.items():res[(a,b)]+=c;res[(b,a)]+=c
for w,c in before.items():res[w]-=c
assert all(S.cancel(c)==0 for c in res.values())
report={'source':'independent R0,W,logs,D primitives','basepoint':'r=sqrt(s)','variable':'y=u','ell':'ln(s)',
 'input_dictionary_slots':counts,'input_nonzero_words':nonzero_counts,'cancelled_nonzero_words':cancelled,'collected_words':len(before)-1,'after_shuffle_words':len(out)-1,
 'products':len(products),'coefficients':{str(w):str(c) for w,c in out.items()},
 'log_products':{str(w):str(c) for w,c in products.items()},'rationalized_sources':sources,'coefficient_roots':'r=sqrt(s), z=i*sqrt(s-1); odd root coefficients are kept until their corresponding words are combined',
 'reconstruction':'Exact algebraic identities after z^2=1-r^2 and r^2=s; exact shuffle reconstruction',
 'scope':'Primitive collection only; finite endpoint and regulator conversions still required','physical_roots':'s>1, r=positive sqrt(s), sqrt(1-s)=i sqrt(s-1); identical letters are combined before continuation to endpoints'}
write_artifact(B.parent/'real_equivalence_allfour.json' , json.dumps(report,indent=2)+'\n')
print('after shuffle',len(out)-1,'GPL and',len(products),'products',flush=True)
