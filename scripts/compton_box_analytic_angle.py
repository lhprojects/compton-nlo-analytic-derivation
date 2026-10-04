"""Build the exact GPL input words for the spanning angular interference.

There is no numerical integration. The closed endpoint evaluation and finite
Laurent reduction are implemented in compton_box_analytic_closed.py.
"""
from compton_artifacts import read_artifact, write_artifact
from collections import defaultdict
from functools import lru_cache
import itertools,json
from pathlib import Path
import sympy as S

v=S.Symbol('v')
x=S.Symbol('x',positive=True)
ll,lx,ds,du=S.symbols('ll lx ds du',real=True)
y=(x-(x+4)*v*v)/(1-v*v)
z=S.sqrt(x/(x+4));eta=S.sqrt((x-1)/(x+3))

def add(*args):
 out=defaultdict(lambda:S.S.Zero)
 for d in args:
  for w,c in d.items():out[w]=S.cancel(out[w]+c)
 return {w:S.cancel(c) for w,c in out.items() if c!=0}
def merge(*args):
 out=defaultdict(lambda:S.S.Zero)
 for d in args:
  for w,c in d.items():out[w]+=c
 return {w:c for w,c in out.items() if c!=0}
def scale(d,r):return {w:S.cancel(r*c) for w,c in d.items() if c!=0}
@lru_cache(None)
def shuffle(a,b):
 if not a:return (b,)
 if not b:return (a,)
 return tuple((a[0],)+w for w in shuffle(a[1:],b))+tuple((b[0],)+w for w in shuffle(a,b[1:]))
def mul(a,b):
 out=defaultdict(lambda:S.S.Zero)
 for w,c in a.items():
  for u,d in b.items():
   for merged in shuffle(w,u):out[merged]+=c*d
 return {w:S.cancel(c) for w,c in out.items() if c!=0}
def build():
 one={():S.S.One};r={(-1,):S.Rational(1,2),(1,):-S.Rational(1,2)}
 logy={():lx,(z,):1,(-z,):1,(1,):-1,(-1,):-1}
 g0=scale(r,(1-v*v)/v)
 ga={(0,-1):(1-v*v)/(2*v),(0,1):-(1-v*v)/(2*v)}
 h={(1,1):(1-v*v)/(2*v),(1,-1):-(1-v*v)/(2*v),(-1,1):(1-v*v)/(2*v),(-1,-1):-(1-v*v)/(2*v)}
 ct=add(scale(g0,ll),scale(h,-S.Rational(1,2)))
 cm=scale(mul(r,r),-(1-v*v)/(2*v*v))
 bt=add(scale(one,2),scale(r,-2/v))
 # du=Li2(1-x)-pi^2/6. The poles eta follow 1-y+i0.
 fu={():du}
 for a,c in [(eta,-1),(-eta,-1),(1,1),(-1,1)]:
  for word,coef in logy.items():fu=add(fu,{(a,)+word:c*coef})
 masters_s={'A':one,'Ba':scale(one,2),'Bs':{():2-x/(1+x)*lx},'Bt':bt,'Bz':{},'Cs':{():ds},'Ct':ct,'Cm':cm,'Db':scale(add(scale(g0,ll-lx),ga),1/x),'constant':one}
 masters_u=dict(masters_s)
 masters_u.update(Bs=add(scale(one,2),scale(logy,y/(1-y))),Cs=scale(fu,1/y),Db=scale(add(scale(g0,ll),scale(mul(g0,logy),-1),ga),-1/y))
 raw=json.loads(read_artifact(Path(__file__).with_name('compton_box_analytic_trace.json')))
 result={}
 yy=S.Symbol('y');xx=S.Symbol('x')
 for channel,ms in [('s',masters_s),('u',masters_u)]:
  terms=[]
  for name,c in raw[channel].items():
   coeff=S.sympify(c).subs({xx:x,yy:y})
   # dcos(theta)=2*s/(s-1)^2 dy. The expression below integrates dy.
   terms.append(scale(ms[name],coeff*8*v/(1-v*v)**2))
  result[channel]=add(*terms)
 return result

def main():
 result={channel:{'representation':'finite Laurent reduction to GPL endpoints',
                   'integrand':[{'word':[str(a) for a in w],'coefficient':str(c)} for w,c in terms.items()]}
         for channel,terms in build().items()}
 write_artifact(Path(__file__).with_suffix('.json') , json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
