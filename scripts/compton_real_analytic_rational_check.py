"""Normalize equivalent algebraic GPL letters and verify their coefficients."""
from compton_artifacts import read_artifact
from pathlib import Path
from collections import defaultdict
import sympy as S,json,time
from compton_quadratic_field import reduce_expression
r,y,H,u=S.symbols('r y H u');G=S.Function('G');base=Path(__file__).parent;p=base/'compton_real_analytic_rational_gpl.json';data=json.loads(read_artifact(p));acc=defaultdict(lambda:S.S.Zero)
def letter(a):
 return a.xreplace({q:S.sqrt(S.expand(q.base)) for q in a.atoms(S.Pow) if q.exp==S.Rational(1,2)})
for key,c in data['coefficients'].items():
 w=() if key=='constant' else tuple(letter(a) for a in S.sympify(key,locals={'G':G}).args[0]);acc[w]+=S.sympify(c)
co={}
for w,c in acc.items():co['constant' if not w else str(G(S.Tuple(*w),y))]=str(S.factor(c))
print('normalized words',len(co)-1,flush=True)
res=defaultdict(lambda:S.S.Zero)
for key,c in co.items():
 c=S.sympify(c);w=() if key=='constant' else tuple(S.sympify(key,locals={'G':G}).args[0]);res[w]+=S.diff(c,y)
 if w:res[w[1:]]+=c/(y-w[0])
K=json.loads(read_artifact(base/'compton_real_analytic_rational_final_kernels.json'));ks={k:S.factor(S.sympify(v).subs({H:r*r-1,u:y})) for k,v in K.items()};s=r*r
words={():ks['rational'],(0,):2*ks['L1']+ks['L2']-ks['L3'],(1,):ks['L2']+ks['L3'],(s,):-ks['L2']-ks['L3'],(-s,):2*ks['L3'],(-1,):-2*ks['L3']}
for w,c in words.items():res[w]-=c
for w,c in res.items():
 t=time.time();q=reduce_expression(c)
 assert not any(q),(w,q)
 print('verified',w,time.time()-t,flush=True)
print('Exact rational-sector GPL derivative verified.',flush=True)
