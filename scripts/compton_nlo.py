"""Inclusive photonic Compton NLO, Lee--Schwartz--Zhang, arXiv:2102.06718.

Equation (6), evaluated through the defining iterated-integral ODEs (7),(8).
Units: GeV, GeV^-2; no e-e+ pair-production final states. x=2*omega/m.
SciPy is needed for a 13-component ODE, mpmath for classical polylogarithms.
The ODE implementation is intended for x >= 0.01, away from x=1 where the
removable logarithmic quotient is evaluated by its continuous limit.
"""
import math
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import spence

ALPHA = 1 / 137.035999084
MASS = 0.00051099895
GEV2_TO_MB = 0.3893793721
# The symbol 'r' is the variable letter with weight dx/sqrt(x*(x+4)).
WORDS = [(-2,-1), (-1,-2,-1), (-1,-1,0), (-1,0,-1), (0,-1,-1),
         (0,'r',-1), ('r',-1,0), ('r',0,-1), ('r',-1,-1),
         ('r',-2,-1), ('r',-1), (-4,'r',-1), ('r','r',-1)]


def iterated_integrals(xmax, rtol=2e-13, atol=2e-14, start=-36):
    def rhs(t, g):
        x = math.exp(t)
        log1 = math.log1p(x)
        # scipy.spence(z)=Li2(1-z). A local series avoids cancellation at x~0.
        li2 = -x + x*x/4 - x**3/9 if x < 1e-5 else spence(1+x)
        lower = {(-1,):log1, (-1,0):t*log1+li2,
                 (0,-1):-li2, (-1,-1):log1*log1/2,
                 (-2,-1):g[0], ('r',-1):g[10]}
        weights = {'r':math.sqrt(x/(x+4)), 0:1.,
                   -1:x/(x+1), -2:x/(x+2), -4:x/(x+4)}
        return [weights[w[0]]*lower[w[1:]] for w in WORDS]
    sol = solve_ivp(rhs, (start, math.log(xmax)), np.zeros(len(WORDS)),
                    method='DOP853', dense_output=True, rtol=rtol, atol=atol,
                    max_step=0.3)
    if not sol.success:
        raise RuntimeError(sol.message)
    return lambda x: dict(zip(WORDS, sol.sol(math.log(x))))


def born(omega):
    r = omega/MASS
    l = math.log1p(2*r)
    return 2*math.pi*ALPHA**2/MASS**2*((1+r)/r**2
          *(2*(1+r)/(1+2*r)-l/r)+l/(2*r)-(1+3*r)/(1+2*r)**2)


def correction(omega, G):
    x = 2*omega/MASS
    y = math.sqrt(x/(x+4))
    l, h = math.log1p(x), math.log(x)
    g = G(x)
    li2 = float(mp.polylog(2, -x))
    li2_1 = float(mp.polylog(2, 1-x))
    # Although each term is complex above x=1, this combination is real.
    with mp.workdps(35):
        xx = mp.mpf(x)
        poly = mp.polylog(3, xx**2)-mp.polylog(2, xx**2)*mp.log(xx)
        assert abs(mp.im(poly)) < mp.mpf('1e-25')
        poly = float(mp.re(poly))
    q = -1. if x == 1 else h/(1-x)
    terms = [
      -x*(273*x**3-982*x*x-2960*x-1744)/(24*(x+1)**2),
      (37*x**4-54*x**3-339*x*x-428*x-184)/(4*(x+1)**2)*l,
      x*x*(14*x**4+17*x**3-17*x*x-22*x-8)/(2*(1+x)**3)*q,
      -(4*x**6+35*x**5-31*x**4-755*x**3-1765*x*x-1506*x-440)/(2*(x+1)**2*(x+4))*l*l,
      -(x**6+7*x**5-28*x**4-239*x**3-449*x*x-338*x-88)/((x+1)**2*(x+4))*li2,
      -(x*x-16*x-23)*l**3/3,
      (x*x-x+2)*(li2_1-math.pi**2/6),
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
    return ALPHA**3/MASS**2/x**3*math.fsum(terms)


def high_energy_correction(omega):
    s = MASS*MASS+2*MASS*omega
    l = math.log(s/MASS**2)
    return ALPHA**3/s*(l**3/3-l*l/2+17*l/4-75/8-math.pi**2/2
                           +4*float(mp.zeta(3)))


if __name__ == '__main__':
    g = iterated_integrals(2*5.6/MASS)
    for w in [4.4,5.,5.47]:
        b, c = born(w), correction(w,g)
        print(w, b*GEV2_TO_MB, (b+c)*GEV2_TO_MB,
              'relative NLO',c/b, 'asymptotic difference',(high_energy_correction(w)-c)/c)
