"""GPL continuation along the causal contour above positive real letters.

This evaluates a closed special-function expression by Taylor continuation.
It does not integrate the original angular or phase-space integrand.
"""
import mpmath as mp
from functools import lru_cache
import sympy as S
from compton_box_analytic_angle import x

class GPLSeries:
 def __init__(self,xx,order=240):self.xx=xx;self.order=order
 @lru_cache(None)
 def letter(self,a):
  z=S.sympify(a).subs(x,S.Float(str(self.xx),mp.mp.dps))
  return mp.mpc(str(S.re(z).evalf(mp.mp.dps)),str(S.im(z).evalf(mp.mp.dps)))
 @lru_cache(None)
 def coefficients(self,word):
  if not word:return (mp.mpf(1),)+(mp.mpf(0),)*self.order
  assert word[-1]!=0
  tail=self.coefficients(word[1:]);a=self.letter(word[0]);out=[mp.mpf(0)]
  if a==0:out.extend(tail[n]/n for n in range(1,self.order+1))
  else:
   derivative=mp.mpf(0)
   for n in range(self.order):
    derivative=(derivative-tail[n])/a;out.append(derivative/(n+1))
  return tuple(out)
 def value(self,word,z):
  if not word:return mp.mpf(1)
  if word==(0,):return mp.log(z)
  return mp.polyval(tuple(reversed(self.coefficients(word))),z)


class GPLContinuation(GPLSeries):
 def values(self,words,endpoint):
  family={()}
  for word in words:
   for j in range(len(word)+1):family.add(word[j:])
  family=sorted(family,key=lambda w:(len(w),str(w)))
  letters={a:self.letter(a) for word in family for a in word}
  start=min(endpoint/20,min(abs(a) for a in letters.values() if a!=0)/10)
  center=start*(1+1j)
  values={word:self.value(word,center) for word in family}
  height=min(endpoint/4,mp.mpf('.2'))
  waypoints=[endpoint/2+1j*height,endpoint+1j*height,mp.mpc(endpoint)]
  steps=0
  for destination in waypoints:
   while abs(destination-center)>mp.mpf(10)**(-mp.mp.dps+5):
    radius=min(abs(center-a) for a in letters.values())
    distance=abs(destination-center)
    delta=(destination-center)*min(mp.mpf(1),mp.mpf('.32')*radius/distance)
    series={():[mp.mpf(1)]+[mp.mpf(0)]*self.order}
    for word in family:
     if not word:continue
     tail=series[word[1:]];inverse=1/(center-letters[word[0]])
     coefficients=[values[word]];derivative=mp.mpf(0)
     for n in range(self.order):
      derivative=(tail[n]-derivative)*inverse
      coefficients.append(derivative/(n+1))
     series[word]=coefficients
    values={word:mp.polyval(list(reversed(series[word])),delta) for word in family}
    center+=delta;steps+=1
  return values,steps
