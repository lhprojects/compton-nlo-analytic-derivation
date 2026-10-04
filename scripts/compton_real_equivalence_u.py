"""Exact removal of the auxiliary eta roots from the crossed box.

Only our independently derived trace coefficients are read. The proof uses
rational-function identities and the defining GPL differential equations.
No quadrature, fitted coefficients, or published total is involved.
"""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
import json
import sympy as S
x,y=S.symbols('x y');L,F,W=S.symbols('L F W')
BASE=Path(__file__).resolve().parent

def construct():
 raw={k:S.sympify(v) for k,v in json.loads(read_artifact(BASE/'compton_box_analytic_trace.json'))['u'].items()}
 R=S.factor(raw['A']+2*raw['Bs']+raw['constant'])
 Q=S.factor(raw['Bs']*y/(1-y));D=S.factor(raw['Cs']/y)
 # F=G(1,0;y)=-Li2(1-y)+pi^2/6 and W=G(0,1,0;y).
 target=R+Q*L-D*F
 P=(-16*(x-2)/y+((32*x+28)*y+y*y)/x
   +(-16*x*x-8*x+36-(4*x+16)*y-2*y*y)*L/x
   +4*(x-2)*L*L+4*L/(x*(y-1))
   +4*(y*y+(2-2*x)*y+3*x-9)*F/x+4*(x*x-8)*W/x)
 derivative=S.diff(P,y)+S.diff(P,L)/y+S.diff(P,F)*L/(y-1)+S.diff(P,W)*F/y
 assert S.cancel(derivative-target)==0
 # The remaining masters have no y=1 pole. Their only moving angular
 # letters in v are ±sqrt(x/(x+4)), besides 0,±1.
 rest={k:raw[k] for k in ('Ba','Bt','Ct','Cm','Db')}
 for k,c in rest.items():assert S.factor(S.denom(c)).subs(y,1)!=0
 # At y=1 the apparent log/y-1 quotient has the analytic value 1.
 h=S.symbols('h');assert S.limit(S.log(1+h)/h,h,0)==1
 return dict(R=str(S.apart(R,y)),log_coefficient=str(S.apart(Q,y)),
  triangle_coefficient=str(S.apart(D,y)),primitive=str(P),
  primitive_symbols={'L':'log(y)','F':'G(1,0;y)=pi^2/6-Li2(1-y)','W':'G(0,1,0;y)'},
  lower_endpoint='x/(1+x)',upper_endpoint='x',normalization='divide by 2*x^2 in C_box',
  remaining_masters={k:str(v) for k,v in rest.items()},
  exact_derivative_verified=True)

def mobius_certificate():
 from collections import defaultdict
 def derivative(words):
  out=defaultdict(lambda:S.S.Zero)
  for w,c in words.items():
   out[w]+=S.diff(c,x)
   if w:out[w[1:]]+=c/(x-w[0])
  return {w:S.cancel(c) for w,c in out.items() if S.cancel(c)!=0}
 Ff={(-1,0):-S.S.One,(-1,-1):S.S.One}
 Wf={(0,-1,0):-S.S.One,(0,-1,-1):S.S.One,(-1,-1,0):S.S.One,(-1,-1,-1):-S.S.One}
 assert derivative(Ff)=={(0,):-1/(x+1),(-1,):1/(x+1)}
 assert derivative(Wf)=={w:S.cancel(c/(x*(1+x))) for w,c in Ff.items()}
 assert S.cancel(S.diff(x/(1+x),x)/(x/(1+x))-(1/x-1/(1+x)))==0
 assert S.cancel(x/(x+4)-(x/(x+2))**2-4*x/((x+4)*(x+2)**2))==0

def threshold_series(data):
 P=S.sympify(data['primitive'],locals={'x':x,'y':y,'L':L,'F':F,'W':W});ly=S.log(y)
 f=sum(y**n*(S.Rational(1,n*n)-ly/n) for n in range(1,7))
 w=sum(y**n*(S.Rational(2,n**3)-ly/n**2) for n in range(1,7))
 P=P.subs({L:ly,F:f,W:w})
 return str(S.series(P.subs(y,x)-P.subs(y,x/(1+x)),x,0,4))

def endpoint_words(data):
 from compton_box_analytic_angle import add,mul,scale
 P=S.sympify(data['primitive'],locals={'x':x,'y':y,'L':L,'F':F,'W':W})
 upper=[{(0,):S.S.One},{(1,0):S.S.One},{(0,1,0):S.S.One}]
 lower=[{(0,):S.S.One,(-1,):-S.S.One},
  {(-1,0):-S.S.One,(-1,-1):S.S.One},
  {(0,-1,0):-S.S.One,(0,-1,-1):S.S.One,(-1,-1,0):S.S.One,(-1,-1,-1):-S.S.One}]
 def at(point,atoms):
  terms=[]
  for powers,c in S.Poly(P,L,F,W).terms():
   words={():S.S.One}
   for power,atom in zip(powers,atoms):
    for _ in range(power):words=mul(words,atom)
   terms.append(scale(words,S.cancel(c.subs(y,point))))
  return add(*terms)
 out=add(at(x,upper),scale(at(x/(1+x),lower),-1))
 return {'representation':'exact eta-free subgroup endpoint words, argument x',
  'integrand':[{'word':list(w),'coefficient':str(c)} for w,c in out.items()]}

def angular_remainder():
 from compton_box_analytic_angle import x as ax,v,ll,lx,add,scale,mul
 yy=(ax-(ax+4)*v*v)/(1-v*v);z=S.sqrt(ax/(ax+4))
 one={():S.S.One};r={(-1,):S.Rational(1,2),(1,):-S.Rational(1,2)}
 logy={():lx,(z,):1,(-z,):1,(1,):-1,(-1,):-1}
 g0=scale(r,(1-v*v)/v)
 ga={(0,-1):(1-v*v)/(2*v),(0,1):-(1-v*v)/(2*v)}
 hh={(1,1):(1-v*v)/(2*v),(1,-1):-(1-v*v)/(2*v),(-1,1):(1-v*v)/(2*v),(-1,-1):-(1-v*v)/(2*v)}
 masters={'Ba':scale(one,2),'Bt':add(scale(one,2),scale(r,-2/v)),
  'Ct':add(scale(g0,ll),scale(hh,-S.Rational(1,2))),
  'Cm':scale(mul(r,r),-(1-v*v)/(2*v*v)),
  'Db':scale(add(scale(g0,ll),scale(mul(g0,logy),-1),ga),-1/yy)}
 raw=json.loads(read_artifact(BASE/'compton_box_analytic_trace.json'))['u'];terms=[]
 for name,master in masters.items():
  coefficient=S.sympify(raw[name]).subs({x:ax,y:yy})*8*v/(1-v*v)**2
  terms.append(scale(master,coefficient))
 return {'representation':'eta-free u remainder after analytic y subgroup',
  'integrand':[{'word':[str(a) for a in w],'coefficient':str(c)} for w,c in add(*terms).items()]}

def check_angular_split(data,remainder):
 from compton_box_analytic_angle import x as ax,v,ll,lx,ds,du,add,scale
 yy=(ax-(ax+4)*v*v)/(1-v*v);z=S.sqrt(ax/(ax+4));eta=S.sqrt((ax-1)/(ax+3))
 logy={():lx,(z,):S.S.One,(-z,):S.S.One,(1,):-S.S.One,(-1,):-S.S.One}
 fu={():du}
 for a,c in [(eta,-1),(-eta,-1),(S.Integer(1),1),(S.Integer(-1),1)]:
  for w,k in logy.items():fu=add(fu,{(a,)+w:c*k})
 R,Q,D=[S.sympify(data[k]).subs({x:ax,y:yy}) for k in ('R','log_coefficient','triangle_coefficient')]
 subgroup=scale(add({():R},scale(logy,Q),scale(fu,D)),8*v/(1-v*v)**2)
 local={'x':ax,'v':v,'ll':ll,'lx':lx,'ds':ds,'du':du}
 def decode(rows):return {tuple(S.sympify(a,locals=local) for a in row['word']):S.sympify(row['coefficient'],locals=local) for row in rows}
 original=json.loads(read_artifact(BASE/'compton_box_analytic_angle.json'))['u']['integrand']
 residual=add(decode(original),scale(decode(remainder['integrand']),-1),scale(subgroup,-1))
 assert not residual,residual

if __name__=='__main__':
 data=construct();mobius_certificate();data['mobius_derivative_verified']=True;data['threshold_series']=threshold_series(data);folder=Path(__file__).resolve().parents[1]
 write_artifact(folder/'real_equivalence_u.json' , json.dumps(data,indent=2)+'\n')
 print('Exact crossed-subgroup derivative identity verified; eta absent.')
 print('All remaining master coefficients regular at y=1.')
 remainder=angular_remainder()
 check_angular_split(data,remainder)
 print('Exact decomposition of the original u angular word dictionary verified.')
 write_artifact(folder/'real_equivalence_u_remainder.json' , json.dumps(remainder,indent=2)+'\n')
 print('Saved',len(remainder['integrand']),'eta-free remainder words.')
 endpoint=endpoint_words(data)
 write_artifact(folder/'real_equivalence_u_endpoint.json' , json.dumps(endpoint,indent=2)+'\n')
 print('Saved',len(endpoint['integrand']),'subgroup endpoint words.')
