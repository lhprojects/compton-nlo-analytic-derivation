"""Differentiate the finite GPL primitives, checking their exact kernels.
No angular quadrature or published result is used.
"""
import mpmath as mp
import sympy as S
from collections import defaultdict
from compton_box_analytic_closed import _Primitive,_specifications
from compton_box_analytic_continuation import GPLContinuation

def run(point,channel):
 with mp.workdps(80):
  xx=mp.mpf(point);lam=mp.mpf('1e-5')
  args=(xx,mp.log(lam),mp.log(xx),-mp.re(mp.polylog(2,1+xx)-mp.pi**2/6)/xx,mp.polylog(2,1-xx)-mp.pi**2/6)
  evaluator=GPLContinuation(xx,180);primitive=_Primitive(evaluator);inputs=[]
  for word,pfun,qfun,poles in _specifications(channel):
   word=tuple(S.S.Zero if evaluator.letter(a)==0 else a for a in word)
   p,q=pfun(*args),qfun(*args);primitive.integrate_rational(word,p,q,poles)
   inputs.append((word,p,q))
  worst=mp.mpf(0)
  for fraction in [mp.mpf('.2'),mp.mpf('.7'),mp.mpf('.93')]:
   z=fraction*xx/(xx+2);derivative=defaultdict(lambda:mp.mpc(0));target=defaultdict(lambda:mp.mpc(0))
   for (word,kind,power,a),c in primitive.terms.items():
    if kind=='poly':r=z**power;dr=power*z**(power-1) if power else 0
    else:
     aa=evaluator.letter(a);r=(z-aa)**(-power);dr=-power*(z-aa)**(-power-1)
    derivative[word]+=c*dr
    if word:derivative[word[1:]]+=c*r/(z-evaluator.letter(word[0]))
   for word,p,q in inputs:target[word]+=mp.polyval(list(reversed(p)),z)/mp.polyval(list(reversed(q)),z)
   scale=max([mp.mpf(1)]+[abs(c) for c in target.values()])
   worst=max(worst,max(abs(derivative[w]-target[w])/scale for w in set(derivative)|set(target)))
  assert worst<mp.mpf('1e-65'),(point,channel,worst)
  return worst
if __name__=='__main__':
 for point in ['.6','1','9','99']:
  for channel in ['s','u']:print(point,channel,mp.nstr(run(point,channel),8),flush=True)
