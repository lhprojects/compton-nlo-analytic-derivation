"""Closed parameter moments for a real-photon on/off-shell QED vertex.

All integrals are evaluated in logarithms and dilogarithms. The moments
multiply the polynomial Dirac numerators in compton_physical_virtual.py.
Units: electron mass = 1. No published Compton correction is used.
"""
from functools import lru_cache
import mpmath as mp


def vertex_moments_mp(S):
    """Return J00,J10,J11,J20,J21,J22 in the current mp precision.

    Jnk = integral_0^1 dv v^k F_n(1+(S-1)v, -(S-1)v).
    The physical boundary is log(1-S-i0), Li2(S+i0).
    """
    S = mp.mpf(S)
    D = S-1
    if D == 0:
        raise ValueError('The on-shell endpoint requires the photon mass before integration')
    if S == 0:
        # Removable poles in the closed primitives, evaluated by their series.
        P = [mp.mpf(0), mp.mpf(1), mp.mpf(-2), mp.mpf(7)/4,
             -mp.mpf(15)/4, mp.mpf(31)/4]
    else:
        ell = mp.log(abs(1-S)) - (mp.j*mp.pi if S>1 else 0)
        li = mp.polylog(2,S)
        if S>1:
            li = mp.conj(li)
        P = [li,
             li+(1-1/S)*ell,
             2*S-2*li+(-S+1/S)*ell-1,
             li+(mp.mpf(3)/2-2/S+1/(2*S*S))*ell+1/(2*S),
             mp.mpf(5)*S/2-3*li+(-S-mp.mpf(3)/2+3/S-1/(2*S*S))*ell-1-1/(2*S),
             S*S-mp.mpf(15)*S/2+6*li+(-S*S/2+4*S-4/S+1/(2*S*S))*ell+4+1/(2*S)]
    bounds = [mp.pi**2/6,mp.pi**2/6,1-mp.pi**2/3,
              mp.mpf(1)/2+mp.pi**2/6,1-mp.pi**2/2,-2+mp.pi**2]
    return [(p-b)/D**power for p,b,power in zip(P,bounds,[1,1,2,1,2,3])]


@lru_cache(maxsize=8192)
def vertex_moments(S):
    # Extra precision controls cancellations at S=1 and the removable S=0 poles.
    scale = min(abs(S),abs(S-1))
    extra = max(0, -int(mp.floor(mp.log10(scale)))) if scale else 0
    with mp.workdps(40+3*extra):
        return tuple(complex(x) for x in vertex_moments_mp(S))
