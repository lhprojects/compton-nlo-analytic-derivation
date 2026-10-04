"""Formal exact high-energy UV expansion with explicit GPL boundary periods.
Unreduced periods remain identified rather than fitted or silently discarded.
"""
from compton_artifacts import read_artifact, write_artifact
import sympy as S,json
from pathlib import Path
from functools import lru_cache
q,z,L=S.symbols('q z L');ROOT=Path(__file__).resolve().parents[1];N=4
periods={};theta=S.Symbol('Cl2_pi_over_3',real=True)
@lru_cache(None)
def constant(word):
 if not word:return S.S.One
 if all(a==0 for a in word) or all(a==1 for a in word):return S.S.Zero
 if len(word)==1:
  if word[0]==-1:return -S.log(2)
  return -S.sign(S.im(word[0]))*S.I*S.pi/3
 if len(word)==2:
  a,b=word
  if a==b:return constant((a,))**2/2
  if a==0 or b==0:
   t=b if a==0 else a
   if t==1:v=S.pi**2/6
   elif t==-1:v=-S.pi**2/12
   elif S.simplify(t*S.conjugate(t)-1)==0:v=S.pi**2/36-S.sign(S.im(t))*S.I*theta
   else:v=S.polylog(2,1/t)
   return v if a==0 else -v
 if len(word)==3 and word.count(0)==2:
  a=next(a for a in word if a!=0)
  if a==1:v=S.zeta(3)
  elif a==-1:v=-3*S.zeta(3)/4
  else:v=S.zeta(3)/3-S.sign(S.im(a))*S.I*5*S.pi**3/162
  return (-2 if word[1]!=0 else 1)*v
 key='P'+str(len(periods));symbol=S.Symbol(key);periods[key]=list(map(str,word));return symbol

def anti(poly,n):
 if n==-1:return S.integrate(poly,L)
 return sum(co*sum((-1)**j*S.factorial(k)/S.factorial(k-j)*L**(k-j)/S.Integer(n+1)**(j+1) for j in range(k+1)) for (k,),co in S.Poly(poly,L).terms())
@lru_cache(None)
def expansion(word):
 if not word:return {0:S.S.One}
 out={0:constant(word)};tail=expansion(word[1:]);a=word[0]
 if a==0:
  for n,p in tail.items():out[n]=out.get(n,0)+anti(p,n-1)
 else:
  for n,p in tail.items():
   for j in range(N-n):out[n+j+1]=out.get(n+j+1,0)-anti(p,n+j)/a**(j+1)
 return {n:S.expand(c) for n,c in out.items() if n<=N}
rows=json.loads(read_artifact(ROOT/'soft_uv_constant_basis_finite.json'))['coefficients']
X=(1-q)**2/q;total=0
for row in rows:
 word=tuple(S.sympify(a) for a in row['word']);co=S.series(S.sympify(row['coefficient'])/(8*(X+1)*X),q,0,2).removeO()
 for n,p in expansion(word).items():total+=co*p*q**n
qt=S.series(total,q,0,2).removeO().expand()
qs=S.series((1+2*z-S.sqrt(1+4*z))/(2*z),z,0,6).removeO()
lqs=S.series(S.log(qs/z),z,0,5).removeO()-L
out=S.series(qt.subs({q:qs,L:lqs}),z,0,2).removeO().expand()
coefficients={str(n):str(S.simplify(out.coeff(z,n))) for n in range(-4,2) if out.coeff(z,n)!=0}
report=dict(normalization='single UV finite contribution INCLUDING 1/(8*s*x)',variable='z=1/x,L=ln x',coefficients=coefficients,boundary_periods=periods,period_definition='P_w = regularized G based at 1 (w;0) = (-1)^len(w) regularized G based at 0(reverse(w);1)',status='Exact formal expansion; listed boundary periods still require analytic reduction',remainder='O(log(x)^3/x^2)')
write_artifact(ROOT/'soft_asymptotic_UV_high_periods.json' , json.dumps(report,indent=2)+'\n')
print(report,flush=True)

# Analytic reduction of the boundary periods.  The relations are shuffle
# identities and integrals of log(1-t+t^2), not numerically inferred fits.
a=(1-S.sqrt(3)*S.I)/2;b=S.conjugate(a)
symbols={tuple(S.sympify(a) for a in w):S.Symbol(name) for name,w in periods.items()}
P=lambda *w:symbols[tuple(w)]
U=sum(P(0,c,d) for c in [a,b] for d in [a,b])
relations={P(a,1):-P(1,a),P(b,1):-P(1,b),
 P(0,a,1):-P(1,0,a)-P(0,1,a),P(0,b,1):-P(1,0,b)-P(0,1,b),
 P(b,a):S.pi**2/9-P(a,b),P(1,b):-S.pi**2/18-P(1,a),
 P(-1,0,b):U-S.zeta(3)/12-P(-1,0,a)}
reduced={}
for n,co in coefficients.items():
 f=S.sympify(co)
 for _ in range(3):f=f.xreplace(relations)
 f=S.simplify(S.expand(f))
 assert not any(f.has(p) for p in symbols.values()),f
 reduced[n]=str(f)
report['coefficients']=reduced
report['status']='Boundary periods reduced analytically by shuffle identities and D(t)(1+t)=1+t^3'
report['period_relations']={str(k):str(v) for k,v in relations.items()}
write_artifact(ROOT/'soft_asymptotic_UV_high.json' , json.dumps(report,indent=2)+'\n')
print('REDUCED',reduced,flush=True)
