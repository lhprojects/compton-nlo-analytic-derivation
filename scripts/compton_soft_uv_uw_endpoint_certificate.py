"""Exact all-H certificate for removal of UW finite t=1 apparent poles.
Only the actual thirteen-term inner finite expression is used. No numerical
sampling, reference cross section, or numerical quadrature enters.
"""
from compton_artifacts import read_artifact, write_artifact
import json,time
from pathlib import Path
import sympy as S
H,t,d,E=S.symbols('H t d E');ROOT=Path(__file__).resolve().parents[1]
source=ROOT/'scripts/compton_real_analytic_UW_finite.json'
records=json.loads(read_artifact(source))['finite']
N=4
trunc=lambda f:S.series(f,d,0,N+1).removeO().expand()
ell=E-sum((H/(1+H))**n*d**n/S.Integer(n) for n in range(1,N+1))
lt=-sum(d**n/S.Integer(n) for n in range(1,N+1))
lp=trunc(sum((-1)**(n+1)*(H*d-H*d*d)**n/S.Integer(n) for n in range(1,N+1)))
z1=-d*(1+H-H*d)/(1-d);z2=-d/(1-d)
d3=trunc(sum((z1**n-z2**n)/S.Integer(n*n) for n in range(1,N+1)))
atoms=[ell,S.Symbol('LH'),S.Symbol('L2'),lp,lt,S.Symbol('D0'),S.Symbol('D2'),d3]
contributions={};sumco={k:S.S.Zero for k in range(-5,0)};regular=[]
for j,row in enumerate(records):
 start=time.time();R=S.sympify(row['coefficient']);num,den=S.fraction(R)
 m=sum(power for factor,power in S.factor_list(den,t)[1] if factor.subs(t,1)==0)
 if not m:regular.append(j);continue
 # Extract the known (t-1)^m factor BEFORE substituting t=1-d. No global
 # common denominator or full transcendental expression is expanded.
 regular_den=S.cancel(den/(t-1)**m)
 pn=S.Poly(S.expand(num.subs(t,1-d)),d)
 pd=S.Poly(S.expand(regular_den.subs(t,1-d)),d)
 assert pd.nth(0)!=0
 nc=[S.factor(pn.nth(k)) for k in range(m)]
 dc=[S.factor(pd.nth(k)) for k in range(m)]
 rational=[]
 for k in range(m):rational.append(S.cancel((nc[k]-sum(dc[l]*rational[k-l] for l in range(1,k+1)))/dc[0]))
 atom=S.S.One
 for power,f in zip(row['powers'],atoms):
  if power:atom=trunc(atom*f**power)
 atom=S.Poly(atom,d)
 out={}
 for k in range(m):
  co=S.factor((-1)**m*sum(rational[l]*atom.nth(k-l) for l in range(k+1)))
  out[str(k-m)]=str(co);sumco[k-m]+=co
 contributions[str(j)]={'powers':row['powers'],'rational_pole_order':int(m),'negative_coefficients':out}
 print('row',j,'pole',m,'seconds',round(time.time()-start,2),flush=True)
certificate={}
for k in range(-5,0):
 # E is independent. Identity in E is stronger than substituting log(1+H).
 poly=S.Poly(S.expand(sumco[k]),E)
 terms={str(power[0]):str(S.factor(co)) for power,co in poly.terms()}
 assert all(S.cancel(co)==0 for _,co in poly.terms()),(k,terms)
 certificate[str(k)]=terms
 print('d power',k,'all E coefficients zero',flush=True)
report={'domain':'H>0, d=1-t ->0+', 'source':str(source.relative_to(ROOT)),
 'independent_symbol':'E=ln(1+H), treated as algebraically independent',
 'method':'Finite Taylor convolution through degree four; exact rational arithmetic',
 'regular_rows':regular,'singular_row_contributions':contributions,
 'atom_series':{'ell':str(ell),'lt':str(lt),'lp':str(lp),'d3':str(d3)},
 'certificate':certificate,'conclusion':'All coefficients of d^-5,...,d^-1 vanish identically for every H>0. The complete UW finite integrand is analytic at t=1.',
 'scope':'Exact cancellation certificate, not a rigorous numerical remainder bound for the production Taylor evaluator.'}
write_artifact(ROOT/'uw_endpoint_certificate.json' , json.dumps(report,indent=2)+'\n')
