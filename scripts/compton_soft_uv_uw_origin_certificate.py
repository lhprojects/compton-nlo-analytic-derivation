"""All-H exact lower-endpoint certificate for the finite UW inner kernel.
The only log(t) cancels before series expansion.  Its dilogarithm boundary
constant is fixed on the negative real axis, as described in the TeX proof.
"""
from compton_artifacts import read_artifact, write_artifact
import sympy as S,json,time
from pathlib import Path
H,t,LH,L2=S.symbols('H t LH L2');ROOT=Path(__file__).resolve().parents[1]
data=json.loads(read_artifact(ROOT/'scripts/compton_real_analytic_UW_finite.json'))
records=data['finite'];N=6
powers_d3=(0,0,0,0,0,0,0,1);powers_ellt=(1,0,0,0,1,0,0,0)
coeffs={tuple(r['powers']):S.sympify(r['coefficient']).subs(S.log(2),L2) for r in records}
assert S.cancel(coeffs[powers_d3]+coeffs[powers_ellt])==0
print('Exact coefficient of ell*log(t): zero',flush=True)
trunc=lambda f:S.series(f,t,0,N+1).removeO().expand()
ell=sum((-1)**(n+1)*(H*t)**n/S.Integer(n) for n in range(1,N+1))
lp=trunc(sum((-1)**(n+1)*(H*t-H*t*t)**n/S.Integer(n) for n in range(1,N+1)))
l2=L2+sum((-1)**(n+1)*(H*t/2)**n/S.Integer(n) for n in range(1,N+1))
# D=d3-ell*log(t): D'=lp/(1-t)-H*lp/(1+Ht)+(lp-ell)/t, D(0)=0.
deriv=S.series(lp/(1-t)-H*lp/(1+H*t)+(lp-ell)/t,t,0,N).removeO().expand()
D=sum(co*t**(p[0]+1)/S.Integer(p[0]+1) for p,co in S.Poly(deriv,t).terms())
d0=sum((-H*t)**n/S.Integer(n*n) for n in range(1,N+1))
d2deriv=S.series(-H*l2/(1+H*t),t,0,N).removeO().expand()
d2=sum(co*t**(p[0]+1)/S.Integer(p[0]+1) for p,co in S.Poly(d2deriv,t).terms())
atoms=[ell,LH,l2,lp,S.Symbol('UNUSED_LOGT'),d0,d2,D]
acc={n:S.S.Zero for n in range(-5,2)};rows={}
for i,row in enumerate(records):
 powers=tuple(row['powers'])
 if powers==powers_ellt:continue
 start=time.time();R=coeffs[powers];nu,de=S.fraction(R);pn=S.Poly(nu,t);pd=S.Poly(de,t)
 m=min(p[0] for p,c in pd.terms());count=m+2
 nc=[S.factor(pn.nth(k)) for k in range(count)]
 dc=[S.factor(pd.nth(k+m)) for k in range(count)]
 quot=[]
 for k in range(count):quot.append(S.cancel((nc[k]-sum(dc[j]*quot[k-j] for j in range(1,k+1)))/dc[0]))
 atom=S.S.One
 for power,f in zip(powers,atoms):
  if power:atom=trunc(atom*f**power)
 atom=S.Poly(atom,t);out={}
 for k in range(count):
  co=S.factor(sum(quot[j]*atom.nth(k-j) for j in range(k+1)))
  out[str(k-m)]=str(co);acc[k-m]+=co
 rows[str(i)]={'powers':list(powers),'rational_pole_order':m,'coefficients':out}
 print('row',i,'pole',m,'seconds',round(time.time()-start,2),flush=True)
certificate={}
for n in range(-5,2):
 poly=S.Poly(S.expand(acc[n]),LH,L2)
 reduced=sum(S.factor(co)*LH**power[0]*L2**power[1] for power,co in poly.terms())
 reduced=S.factor(reduced);certificate[str(n)]=str(reduced)
 if n<0:assert reduced==0,(n,reduced)
 print('coefficient t^',n,':',reduced,flush=True)
report={'domain':'All fixed H>0; real t->0+', 'log_cancellation':'R_d3+R_ell_logt=0 exactly',
 'analytic_combination':'D=d3-log(1+Ht)*log(t)',
 'D_derivative':'log(P)/(1-t)-H*log(P)/(1+Ht)+(log(P)-log(1+Ht))/t',
 'D_boundary':'D(0)=0, established with negative-real-axis Li2 inversion',
 'D_series':str(S.collect(D,t)),'row_coefficients':rows,'combined_coefficients':certificate,
 'result':'Every negative power vanishes exactly; no log(t) remains in the combined finite kernel.',
 'scope':'Lower finite endpoint analytic certificate. It does not bound the production numerical Taylor remainder.'}
write_artifact(ROOT/'uw_origin_certificate.json' , json.dumps(report,indent=2)+'\n')
