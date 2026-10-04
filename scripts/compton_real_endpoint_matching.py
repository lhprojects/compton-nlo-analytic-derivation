"""Closed conversion of real z-cut finite parts to the common (xy)^eta scheme.

The input finite parts are coefficients before the overall 1/(8*s*(s-1)).
This module integrates no phase-space variable. It evaluates only the exact
logarithmic/dilogarithmic endpoint moments derived in the supplement.
"""
from compton_artifacts import read_artifact
import json
from functools import lru_cache
from pathlib import Path
import sympy as S
import mpmath as mp
from compton_box_analytic_matching_check import bose_matching
H=S.Symbol('H')


@lru_cache(maxsize=1)
def expressions():
    data=json.loads(read_artifact(Path(__file__).with_name('compton_real_soft_matching_integrated.json')))
    out={name:{key:S.sympify(data[name][key]) for key in ('pole','finite_shift')}
         for name in ('rational','logs','LW')}
    b=bose_matching(H)
    out['bose']={'pole':b['A'],'finite_shift':b['shift']}
    return out


def endpoint_matching_closed(h,dps=50):
    """Return poles and finite shifts. No finite energy integral is included."""
    with mp.workdps(dps+30):
        hn=mp.mpf(h)
        if hn<=0:raise ValueError('H=s-1 must be positive')
        return {name:{key:+S.lambdify(H,value,'mpmath')(hn) for key,value in values.items()}
                for name,values in expressions().items()}


def combine_real_finite(z_finite,uv_finite,uw_finite,h,dps=50):
    """Combine all independently supplied closed energy finite parts.

    z_finite must contain rational, logs, LW, bose in the z-cut convention.
    UV and UW are already in the symmetric eta convention. The Bose restoring
    term makes doubling UW valid. This function never fills a missing sector.
    """
    required={'rational','logs','LW','bose'}
    if set(z_finite)!=required:raise ValueError('All four z-cut finite parts are required')
    with mp.workdps(dps+30):
        shifts=endpoint_matching_closed(h,dps)
        finite=mp.fsum(mp.mpf(z_finite[k])+shifts[k]['finite_shift'] for k in required)
        finite+=mp.mpf(uv_finite)+2*mp.mpf(uw_finite)
        return finite/(8*mp.mpf(h)*(mp.mpf(h)+1))
