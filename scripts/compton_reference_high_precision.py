"""High-precision evaluation of LSZ Eq. (6), used ONLY as a reference.

The independent amplitude integration never imports this module.
"""
import mpmath as mp
from compton_nlo import WORDS


def coefficient(s,dps=40,log_lower=-80):
    with mp.workdps(dps):
        return +_coefficient(mp.mpf(str(s))-1,log_lower)


def _coefficient(x,log_lower=-80):
    y=mp.sqrt(x/(x+4));l=mp.log1p(x);h=mp.log(x)
    li2=mp.polylog(2,-x);li2_1=mp.polylog(2,1-x)
    poly=mp.re(mp.polylog(3,x*x)-mp.polylog(2,x*x)*h)
    q=-mp.mpf(1) if x==1 else h/(1-x)
    def lower(word,z):
        if word==(-1,):return mp.log1p(z)
        if word==(-1,0):return mp.log(z)*mp.log1p(z)+mp.polylog(2,-z)
        if word==(0,-1):return -mp.polylog(2,-z)
        if word==(-1,-1):return mp.log1p(z)**2/2
        if word==(-2,-1):return mp.polylog(2,-z-1)+mp.log1p(z)*mp.log(z+2)+mp.pi**2/12
        if word==('r',-1):
            w=4/(mp.sqrt(z+4)+mp.sqrt(z))**2
            return mp.polylog(2,-w**3)/3-mp.polylog(2,-w)+mp.log(w)**2/2-mp.pi**2/18
        raise ValueError(word)
    g={}
    knots=sorted(set([mp.mpf(log_lower)]+[mp.mpf(v) for v in [-20,-5,0,3,6] if v<h]+[h]))
    for word in WORDS:
        def f(t):
            z=mp.exp(t);a=word[0]
            weight=mp.sqrt(z/(z+4)) if a=='r' else z/(z-a)
            return weight*lower(word[1:],z)
        g[word]=mp.quad(f,knots)
    terms = [
      -x*(273*x**3-982*x*x-2960*x-1744)/(24*(x+1)**2),
      (37*x**4-54*x**3-339*x*x-428*x-184)/(4*(x+1)**2)*l,
      x*x*(14*x**4+17*x**3-17*x*x-22*x-8)/(2*(1+x)**3)*q,
      -(4*x**6+35*x**5-31*x**4-755*x**3-1765*x*x-1506*x-440)/(2*(x+1)**2*(x+4))*l*l,
      -(x**6+7*x**5-28*x**4-239*x**3-449*x*x-338*x-88)/((x+1)**2*(x+4))*li2,
      -(x*x-16*x-23)*l**3/3,
      (x*x-x+2)*(li2_1-mp.pi**2/6),
      (x**4+7*x**3+x*x-3*x-2)/(x+1)**2*l*h,
      (x*x+2*x-6)*poly,
      -4*(x**5+26*x**4+146*x**3+316*x*x+288*x+96)/((x+1)**2*(x+4))*g[-2,-1],
      8*(x*x-4*x-6)*g[-1,-2,-1],
      4*(2*x*x-x-6)*g[-1,-1,0],
      2*(2*x*x-7*x-12)*g[-1,0,-1],
      -(5*x*x+32*x-8)*g[0,-1,-1],
      -3*(x-2)*(x+4)*y*(g[0,'r',-1]+2*g['r',-1,0]),
      -8*y*(x**4+3*x**3-18*x*x-68*x-24)/((x+4)*x)*g['r',0,-1],
      3*y*(5*x**4+14*x**3-96*x*x-352*x-128)/((x+4)*x)*g['r',-1,-1],
      -16*y*(x**4+2*x**3-24*x*x-80*x-48)/((x+4)*x)*g['r',-2,-1],
      (3*x**4+18*x**3+44*x*x-8*x-64)/x*y*g['r',-1],
      -6*y*(x**3-12*x-8)/x*g[-4,'r',-1],
      3*(3*x*x-8)*g['r','r',-1]]
    return mp.fsum(terms)/x**3
