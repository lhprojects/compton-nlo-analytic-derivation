"""Exact derivative check in the two quadratic coefficient extensions."""
from compton_artifacts import read_artifact
from pathlib import Path
from collections import defaultdict
import sympy as S,json,time
from compton_quadratic_field import reduce_expression
r,y,H,u,j,k=S.symbols('r y H u j k');G=S.Function('G');base=Path(__file__).parent
co=json.loads(read_artifact(base/'compton_real_analytic_rational_gpl.json'))['coefficients'];res=defaultdict(lambda:S.S.Zero)
for key,c in co.items():
 c=S.sympify(c);w=() if key=='constant' else tuple(S.sympify(key,locals={'G':G}).args[0]);res[w]+=S.diff(c,y)
 if w:res[w[1:]]+=c/(y-w[0])
K=json.loads(read_artifact(base/'compton_real_analytic_rational_final_kernels.json'));ks={a:S.sympify(v).subs({H:r*r-1,u:y}) for a,v in K.items()};s=r*r
words={():ks['rational'],(0,):2*ks['L1']+ks['L2']-ks['L3'],(1,):ks['L2']+ks['L3'],(s,):-ks['L2']-ks['L3'],(-s,):2*ks['L3'],(-1,):-2*ks['L3']}
for w,c in words.items():res[w]-=c
for w,c in res.items():
 start=time.time()
 reduced=reduce_expression(c)
 assert not any(reduced),(w,reduced)
 if not w:
  assert any(reduce_expression(c+1)), "negative control must reject an incorrect constant derivative"
 print('verified',w,'seconds',time.time()-start,flush=True)
print('All rational-sector derivatives exactly verified in the quadratic field.',flush=True)
