"""Evaluate the analytic GPL answer by a convergent power series.

At the chosen low-energy point the full angular endpoint is strictly inside
all GPL-letter radii. Independent 1D quadrature is a cross-check only.
"""
from compton_artifacts import read_artifact
import json
from pathlib import Path
from functools import lru_cache
import mpmath as mp
import sympy as S
from compton_box_analytic_closed import box_angular_components

from compton_box_analytic_continuation import GPLSeries

def value(channel,xx,lam,order=240):
    return box_angular_components(1+xx,lam,dps=mp.mp.dps)[0 if channel=='s' else 1]

def quadrature(channel,xx,lam):
 yy=S.Symbol('y');xxsym=S.Symbol('x')
 data=json.loads(read_artifact(Path(__file__).with_name('compton_box_analytic_trace.json')))[channel]
 fs={key:S.lambdify((xxsym,yy),S.sympify(c),'mpmath') for key,c in data.items()}
 def integrand(y):
  ss=1+xx if channel=='s' else 1-y;D=ss-1;T=xx-y
  a=mp.sqrt(T)/2;h=mp.sqrt(1+a*a);r=mp.asinh(a);k=a/h
  g0=r/(a*h);ga=(mp.polylog(2,k)-mp.polylog(2,-k))/(2*a*h)
  hlog=(2*r*mp.log(2*h)-r*r-mp.polylog(2,-mp.exp(-2*r))+mp.polylog(2,-1))/(a*h)
  ms={'A':1,'Ba':2,'Bs':2+(1-ss)/ss*mp.log(abs(1-ss)), 'Bt':2-2*h*r/a,'Bz':0,
      'Cs':-mp.re(mp.polylog(2,ss)-mp.pi**2/6)/(ss-1), 'Ct':g0*mp.log(lam)-hlog/2,'Cm':-2*r*r/T,
      'Db':(g0*mp.log(lam/abs(D))+ga)/D,'constant':1}
  return sum(fs[key](xx,y)*ms[key] for key in fs)
 return mp.quad(integrand,[xx/(1+xx),xx])

if __name__=='__main__':
 mp.mp.dps=55;xx=mp.mpf('.6');lam=mp.mpf('1e-5')
 for channel in ['s','u']:
  analytic=value(channel,xx,lam)
  independent=quadrature(channel,xx,lam)
  print(channel,'GPL',mp.nstr(analytic,48),flush=True)
  print(channel,'check',mp.nstr(independent,48),'difference',mp.nstr(analytic-independent,10),flush=True)
  assert abs(analytic-independent)<mp.mpf('1e-40'), (channel,analytic-independent)
