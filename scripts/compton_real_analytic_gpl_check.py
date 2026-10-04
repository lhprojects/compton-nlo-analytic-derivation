"""Exact derivative verification of factored GPL coefficients."""
from compton_artifacts import read_artifact
from pathlib import Path
import sympy as S,json,time
from collections import defaultdict
H,u,r,y=S.symbols('H u r y');G=S.Function('G');base=Path(__file__).parent
D=json.loads(read_artifact(base/'compton_real_analytic_single_gpl.json'))['coefficients']
deriv=defaultdict(lambda:S.S.Zero)
for key,c in D.items():
 c=S.sympify(c)
 if key=='constant':word=()
 else:word=tuple(S.sympify(key,locals={'G':G}).args[0])
 deriv[word]+=S.diff(c,y)
 if word:deriv[word[1:]]+=c/(y-word[0])
src=json.loads(read_artifact(base/'compton_real_analytic_single_sector.json'))
A=S.factor(S.sympify(src['u_rational_single_log']).subs({H:r*r-1,u:y}))
B=S.factor(S.sympify(src['u_rational_two_logs']).subs({H:r*r-1,u:y}))
for c,R,w in [(2,A,(0,)),(4,B,(0,0)),(2,B,(0,1)),(2,B,(1,0)),(-2,B,(0,r*r)),(-2,B,(r*r,0))]:deriv[w]-=c*R
for w,c in deriv.items():
 t=time.time();v=S.factor(c);assert v==0,(w,v)
 print('verified',w,'seconds',time.time()-t,flush=True)
print('All exact GPL derivative identities hold.',flush=True)
