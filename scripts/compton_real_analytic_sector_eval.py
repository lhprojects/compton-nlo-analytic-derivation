"""Analytic GPL endpoint values and z-cut finite parts of real sectors.

No quadrature is used: GPLs are evaluated by convergent power series and
analytic continuation, with an analytic logarithmic endpoint tail.
"""
from compton_artifacts import read_artifact
from pathlib import Path
from functools import lru_cache
import json,math
import sympy as S
import mpmath as mp
from compton_soft_uv_check import word_series
H,u,r,y,b=S.symbols('H u r y b');G=S.Function('G');base=Path(__file__).parent
names={'W':'single_gpl','logs':'logs_gpl','rational':'rational_gpl'}
@lru_cache(None)
def records(sector):
 d=json.loads(read_artifact(base/f'compton_real_analytic_{names[sector]}.json'));out=[]
 for key,c in d['coefficients'].items():
  word=() if key=='constant' else tuple(S.sympify(key,locals={'G':G}).args[0]);expr=S.sympify(c)
  out.append((word,expr,S.lambdify((H,r,y),expr,'mpmath')))
 return out

def gpl_continue(expr,Z,side=1):
 words=set()
 for w in expr:
  while w:words.add(w);w=w[1:]
 words=sorted(words,key=lambda w:(len(w),str(w)));values={w:mp.mpc(0) for w in words};letters=set(w[0] for w in words)
 order=int(mp.ceil((mp.mp.dps+15)/(-mp.log10(mp.mpf('.35')))))
 interior=any(abs(mp.im(a))<mp.eps and 0<mp.re(a)<Z for a in letters)
 path=[side*mp.j*Z/3,Z+side*mp.j*Z/3,mp.mpc(Z)] if interior else [mp.mpc(Z)]
 center=mp.mpc(0);steps=0
 for target in path:
  while abs(target-center)>mp.eps*max(1,abs(target))*16:
   distance=min(abs(center-a) for a in letters if center!=a)
   length=min(abs(target-center),mp.mpf('.35')*distance)
   step=length*(target-center)/abs(target-center)
   series={():[mp.mpf(1)]+[mp.mpf(0)]*order}
   for w in words:
    tail=series[w[1:]];a=w[0];co=[values[w]]
    if center==a:
     assert center==0 and tail[0]==0
     co.extend(tail[n]/n for n in range(1,order+1))
    else:
     den=center-a
     for n in range(order):co.append((tail[n]-n*co[n])/(den*(n+1)))
    series[w]=co
   for w in words:values[w]=mp.polyval(list(reversed(series[w])),step)
   center+=step;steps+=1
   if steps>2000:raise RuntimeError('GPL series continuation did not finish')
 return sum(c*(values[w] if w else 1) for w,c in expr.items()),steps

def numeric_letter(a,h,rr):
 return S.lambdify((H,r),a,'mpmath')(h,rr)-rr

def primitive_value(h,sector,point,side=1):
 rr=mp.sqrt(1+h);expr={}
 for word,c,fn in records(sector):
  w=tuple(mp.mpc(numeric_letter(a,h,rr)) for a in word)
  expr[w]=expr.get(w,0)+fn(h,rr,point)
 value,steps=gpl_continue(expr,point-rr,side)
 return value*(mp.log(1+h) if sector=='logs' else 1),steps

def base_constant(h,sector):
 # Exact local GPL-series constants. Weight-two words vanish to second
 # order and hence give no constant at the at-most-simple basepoint poles.
 if sector in ('W','logs'):return mp.mpf(0)
 if sector=='rational':
  s=h+1
  return 2*s*(-s**5-5*s**4+18*s**3+122*s**2+151*s-45)/(3*(s**2+2*s-3))
 raise ValueError(sector)

def base_constant_from_series(h,sector):
 rr=mp.sqrt(1+h);value=mp.mpc(0)
 for word,c,fn in records(sector):
  if not word:
   value+=fn(h,rr,rr);continue
  den=S.factor(S.denom(c));shifted=S.Poly(den.subs(y,r+b),b);m=min(k[0] for k,coef in shifted.terms())
  if not m:continue
  smooth=S.factor(c*(y-r)**m);f=S.lambdify((H,r,y),smooth,'mpmath')
  w=tuple(mp.mpc(numeric_letter(a,h,rr)) for a in word);gs=word_series(w,m)
  for n in range(1,m+1):
   if gs[n]:value+=gs[n]*mp.diff(lambda z:f(h,rr,rr+z),mp.mpf(0),m-n)/mp.factorial(m-n)
 return value*(mp.log(1+h) if sector=='logs' else 1)

@lru_cache(None)
def kernel_expressions(sector):
 if sector=='W':
  d=json.loads(read_artifact(base/'compton_real_analytic_single_sector.json'))
  return [('L1',S.sympify(d['u_rational_single_log'])),('L1L2',S.sympify(d['u_rational_two_logs']))]
 if sector=='logs':
  d=json.loads(read_artifact(base/'compton_real_analytic_logs_inner.json'));return [(k,S.sympify(d[k])) for k in ('rational','L1','L2')]
 d=json.loads(read_artifact(base/'compton_real_analytic_rational_final_kernels.json'));return [(k,S.sympify(v)) for k,v in d.items()]

def rational_series(expr,h,degree):
 hh=S.Rational(str(h));ss=hh+1
 num,den=S.fraction(S.factor(expr.subs({H:hh,u:ss-b})))
 def poly(p):
  return {int(k[0]):mp.mpf(str(S.N(c,mp.mp.dps+5))) for k,c in S.Poly(p,b).terms()}
 N,D=poly(num),poly(den);n0=min(N);d0=min(D);offset=n0-d0
 nn={k-n0:v for k,v in N.items()};dd={k-d0:v for k,v in D.items()};values=[]
 for n in range(degree-offset+1):values.append((nn.get(n,0)-sum(dd.get(k,0)*values[n-k] for k in range(1,n+1)))/dd[0])
 return {k+offset:v for k,v in enumerate(values)}

def end_series(h,sector,degree):
 s=h+1
 L1={0:(mp.log(s),0)};L2={0:(mp.log(s*h),-1)};L3={0:(mp.log(4*s*h/(s+1)**2),-1)}
 for n in range(1,degree+3):
  L1[n]=(-2/(n*s**n),0);L2[n]=(-1/(n*s**n)-1/(n*h**n),0)
  L3[n]=(-1/(n*h**n)-2/(n*(2*s)**n)+1/(n*s**n)+2/(n*(s+1)**n),0)
 prod={}
 for i,(a0,a1) in L1.items():
  for j,(c0,c1) in L2.items():
   if i+j>degree+2:continue
   old=prod.get(i+j,(0,0));prod[i+j]=(old[0]+a0*c0,old[1]+a0*c1)
 logs={'rational':{0:(1,0)},'L1':L1,'L2':L2,'L3':L3,'L1L2':prod};out={}
 for kind,expr in kernel_expressions(sector):
  for power,coef in rational_series(expr,h,degree).items():
   for n,(c,d) in logs[kind].items():
    if power+n>degree:continue
    old=out.get(power+n,(0,0));out[power+n]=(old[0]+coef*c,old[1]+coef*d)
 if sector=='logs':out={k:(a*mp.log(s),b*mp.log(s)) for k,(a,b) in out.items()}
 return out

def real_sector_zcut(h,sector='W',dps=35):
 with mp.workdps(dps+60):
  h=mp.mpf(h);s=h+1;rr=mp.sqrt(s);delta=mp.mpf('1e-5')*min(1,h,s-rr);degree=int(mp.ceil(dps/5))+6
  expansion=end_series(h,sector,degree);pole,logpole=expansion.get(-1,(0,0))
  assert abs(logpole)<mp.mpf(10)**(-dps-20)*max(1,abs(pole)),logpole
  for n,(c,d) in expansion.items():
   if n < -1:assert max(abs(c),abs(d)) < mp.mpf(10)**(-dps-20)*max(1,abs(pole)), (n,c,d)
  point=s-delta;v,steps=primitive_value(h,sector,point);origin=base_constant(h,sector)
  tail=mp.mpf(0)
  for n,(c,d) in expansion.items():
   if n<0:continue
   tail+=delta**(n+1)*(c/(n+1)+d*(mp.log(delta)/(n+1)-1/mp.mpf(n+1)**2))
  finite=v-origin+pole*mp.log(delta)+tail+pole*mp.log(h/s)
  return {'pole':pole,'finite_zcut':finite,'gpl_steps':steps,'base_constant':origin,'imaginary_residual':mp.im(finite)}

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--H',default='9');p.add_argument('--sector',default='W',choices=names);p.add_argument('--dps',type=int,default=30);args=p.parse_args()
 print({k:mp.nstr(v,args.dps) for k,v in real_sector_zcut(args.H,args.sector,args.dps).items()},flush=True)
