"""Exact elimination of the artificial s=9 root pair from our R0 primitive.
No literature formula or numeric fitted coefficient is used.
"""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
from collections import defaultdict
import json
import sympy as S
BASE=Path(__file__).resolve().parent
r,y=S.symbols('r y');G=S.Function('G');ss=r*r
qminus=(ss-3-S.sqrt(r**4-10*r*r+9))/2
qplus=(ss-3+S.sqrt(r**4-10*r*r+9))/2
C=128*ss*(ss+5)/(3*(ss-1))
k={S.Integer(0):-1,S.Integer(1):1,ss:-1,-ss:2,S.Integer(-1):-2}
data=json.loads(read_artifact(BASE/'compton_real_analytic_rational_gpl.json'))
co=defaultdict(lambda:S.S.Zero)
removed=[]
for key,c in data['coefficients'].items():
 w=() if key=='constant' else tuple(S.sympify(key,locals={'G':G}).args[0])
 c=S.sympify(c)
 if w and (S.simplify(w[0]-qminus)==0 or S.simplify(w[0]-qplus)==0):
  assert len(w)==2 and w[1] in k
  assert S.cancel(c-C*k[w[1]])==0
  removed.append(w)
 else:co[w]+=c
assert len(removed)==10
# Q'/Q=dlog(1-A)-dlog(y^2-s)+dlog[y(s-y)(y+1)^2].
A=(y-1)*(y+ss)**2/(y*(y+1)**2*(ss-y))
Q=y*y-(ss-3)*y+ss
correction={r:-1,-r:-1,S.Integer(0):1,ss:1,S.Integer(-1):2}
assert S.cancel(S.diff(Q,y)/Q-S.diff(1-A,y)/(1-A)-sum(v/(y-a) for a,v in correction.items()))==0
assert S.cancel(A.subs(y,r)-1)==0
for a,v in correction.items():
 for b,w in k.items():co[(a,b)]+=C*v*w
out={('constant' if not w else str(G(S.Tuple(*w),y))):str(S.factor(c)) for w,c in co.items() if c!=0}
result={'basepoint':'r=sqrt(s)>1','variable':'y=u','coefficients':out,'dilogarithm':{'argument':str(S.factor(1-A)),'coefficient':str(-C)},'identity':'sum_{q=qminus,qplus} sum_b k_b G(q,b;y) = -Li2(1-A(y)) + sum_{a,b} correction_a k_b G(a,b;y)','certificates':{'ten_input_coefficients':'exact zero residual','rational_differential_identity':'exact zero residual','basepoint_A':'exactly one; both primitives vanish'},'physical_branch':'A>0 for sqrt(s)<y<s, hence 1-A<1 and Li2 is real; all GPL combinations use the inherited real lower-endpoint continuation. No s=9 branch remains.','note':'This removes qplus/qminus but introduces the already existing Gram letters +/-sqrt(s). It is not a claim that all auxiliary letters have disappeared.'}
path=BASE.parent/'real_equivalence_qpair.json';write_artifact(path , json.dumps(result,indent=2)+'\n')
print('Removed',len(removed),'qpair words; all exact differential/basepoint checks passed.',flush=True)
