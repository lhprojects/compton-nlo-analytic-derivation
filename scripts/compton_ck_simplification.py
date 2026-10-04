"""Remove both artificial imaginary poles from our finite soft angular term."""
from compton_artifacts import write_artifact
from pathlib import Path
from collections import defaultdict
import json
import sympy as S
import mpmath as mp
from compton_soft_analytic import soft_log_letters
from compton_soft_analytic_gpl import primitive,combine
x,u,w=S.symbols('x u w',positive=True)
s=x+1;beta=x/(x+2);E2=(x+2)**2/(4*s);den=1-beta**2*u**2
p=2*(x+2)*u*(1-u*u)/x
b1=(x**3-2*x*x-16*x-8)/x**2
h0=u*((x*x+2*x+4)/(2*s)-x**3/(4*s*s)-3*x**3*(1-u*u)/(8*s*s))
delta=x*(1-u*u)/s
h1=-(x**3-6*x*x-16*x-8)/x**2+delta*(x*x+2*x+4)/(2*x)-3*x*delta**2/8

def certificate():
 A,F,k=S.symbols('A F k');K=(1+x*u*u)/s;W=A/den
 Wprime=beta/den**2+2*beta**2*u*A/den**2
 # Coefficients of F/k: use k^2=K after differentiation.
 old=W*(b1+2*(x+2)*(1-u*u)/(x*s*K))
 boundary_derivative=(Wprime*p+W*S.diff(p,u))-W*p*x*u/(s*K)
 new=A*(b1/den+16*s*u*u/(x*(x+2)*den**2))-2*u*(1-u*u)/den**2
 assert S.cancel(old-boundary_derivative-new)==0
 assert S.cancel(W*2*(x+2)*u/(x*K)-W*p/(K*(1-u*u)))==0
 kernels=symbolic_kernels();R0,R1=kernels[:2]
 assert S.residue(R1,w,0)==0
 P1=-sum(S.residue(R1*w**(j-1),w,0)/(j-1)/w**(j-1) for j in range(2,5))
 U=(w-1/w)/(2*S.sqrt(x))
 assert S.factor(S.residue(R0-P1*S.diff(U,w)/(1-U*U),w,0))==0
 return True

def symbolic_kernels():
 U=(w-1/w)/(2*S.sqrt(x));jac=S.diff(U,w);J=S.sqrt(s/x)/w
 common=4/(beta*E2)
 expr=[jac*h0/den,jac*h1/den,J*(b1/den+16*s*u*u/(x*(x+2)*den**2)),J*(-2*u*(1-u*u)/den**2)]
 return [S.factor(common*f.subs(u,U)) for f in expr]

def partials(ss):
 # Rational residues are finite derivatives, not quadrature.
 xx=mp.mpf(ss)-1;b=xx/(xx+2);a=mp.sqrt(xx)/b;r=1/(mp.sqrt(1+a*a)+a)
 result=[]
 for index,expr in enumerate(symbolic_kernels()):
  R=expr.subs(x,S.Rational(str(xx)));num,denom=map(lambda z:S.Poly(z,w),S.fraction(S.cancel(R)))
  poly,_=S.div(num,denom);Z=S.Symbol('z');poly=S.Poly(poly.as_expr().subs(w,Z+1),Z)
  cv=lambda a:mp.mpf(str(S.N(a,mp.mp.dps+5)))
  terms=[('p',int(power[0]),cv(c)) for power,c in poly.terms() if c]
  zeros=min(n[0] for n,c in denom.terms());factors=S.factor_list(denom.as_expr(),w)[1]
  nonzero=[(f,n) for f,n in factors if S.degree(f,w)>0 and f!=w]
  orders={n for f,n in nonzero};assert len(orders)==1
  order=next(iter(orders));roots=[(a,order) for a in [r,-r,1/r,-1/r]]
  if zeros:roots.append((mp.mpf(0),zeros))
  assert sum(n for a,n in roots)==denom.degree()
  nc=[cv(a) for a in num.all_coeffs()];lead=cv(denom.LC())
  for a,n in roots:
   def reduced(v):
    val=mp.polyval(nc,v)/lead
    for b,m in roots:
     if b!=a:val/=(v-b)**m
    return val
   for j in range(1,n+1):
    if index==1 and a==0 and j==1:continue # identically zero residue, proved in certificate
    terms.append(('f',(a-1,j),mp.diff(reduced,a,n-j)/mp.factorial(n-j)))
  result.append(terms)
 return result

def expression(ss,Z):
 parts=partials(ss);B,L,K=soft_log_letters(mp.mpf(ss));out=defaultdict(mp.mpc)
 for a,ca in B:
  for kind,data,c in parts[0]:combine(out,primitive(kind,data,(a,),Z),ca*c)
  for j,letters in [(1,L),(2,K)]:
   for b,cb in letters:
    for kind,data,c in parts[j]:
     combine(out,primitive(kind,data,(a,b),Z),ca*cb*c)
     combine(out,primitive(kind,data,(b,a),Z),ca*cb*c)
 for a,ca in K:
  for kind,data,c in parts[3]:combine(out,primitive(kind,data,(a,),Z),ca*c)
 # The second exact residue identity in certificate removes these four words.
 for a,ca in B:
  out.pop((mp.mpf(-1),a),None)
 return dict(out)

if __name__=='__main__':
 certificate();kernels=symbolic_kernels()
 for f in kernels:
  assert S.rem(S.denom(f),w*w+1,w)!=0
  print(S.factor(S.denom(f)),flush=True)
 folder=Path(__file__).resolve().parents[1]
 write_artifact(folder/'ck_simplification.json' , json.dumps(dict(exact_integration_by_parts=True,kernels=[str(f) for f in kernels],boundary='W*p*atanh(U/k)/k vanishes at U=0 and U=1',removed_poles=['+i','-i']),indent=2)+'\n')

def closed(s,dps=40):
 """The simplified full CK endpoint, evaluated with the existing tail series.

 The finite-cutoff boundary is added exactly before reusing the original
 endpoint tail. It vanishes at the mathematical endpoint U=1.
 """
 from compton_soft_analytic_gpl import evaluate_gpl_taylor
 from compton_soft_uv_interface import endpoint_series
 from compton_soft_analytic import born_angular
 with mp.workdps(dps+40):
  ss=mp.mpf(s);xx=ss-1;bb=xx/(xx+2);ee=(xx+2)**2/(4*ss)
  degree=8;delta=mp.power(10,-mp.mpf(dps+15)/(degree+1))/(ss+1)
  U=1-delta;k=mp.sqrt((1+xx*U*U)/ss);Z=mp.sqrt(xx)*U+mp.sqrt(1+xx*U*U)-1
  value,_,_=evaluate_gpl_taylor(expression(ss,Z),Z)
  p=2*(xx+2)*U*(1-U*U)/xx
  boundary=4/(bb*ee)*mp.atanh(bb*U)/(1-bb*bb*U*U)*p*mp.atanh(U/k)/k
  aa,cc=endpoint_series(ss,degree+2)
  def tail(n):return sum(delta**(j+1)/(j+1)*(aa[j]*(mp.log(delta)-1/mp.mpf(j+1))+cc[j]) for j in range(n+1))
  v=value+boundary+tail(degree)-mp.atanh(bb)/bb*born_angular(ss)
  assert abs(tail(degree+2)-tail(degree))<mp.mpf(10)**(-dps-2)*max(1,abs(v))
  assert abs(mp.im(v))<mp.mpf(10)**(-dps-2)*max(1,abs(v))
  return +mp.re(v)
