"""Check the derived low-energy truncations against the full analytic total."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S
import mpmath as mp
from compton_total_energy_closed import closed

folder=Path(__file__).resolve().parents[1]
source=json.loads(read_artifact(folder/'total_energy_low.json'))
ell=S.Symbol('ell')
coeff={int(n):S.lambdify(ell,S.sympify(c),'mpmath') for n,c in source['coefficients'].items()}
if __name__=='__main__':
    mp.mp.dps=95;rows=[]
    for xx in ['0.01','0.001','0.0001']:
        x=mp.mpf(xx);full=closed(x,60);errors={}
        for order in [2,3,4]:
            approx=mp.fsum(fn(mp.log(x))*x**n for n,fn in coeff.items() if n<=order)
            errors[order]=abs(full-approx)
        assert errors[4]<errors[3]<errors[2]
        scaled=errors[4]/(x**5*abs(mp.log(x)))
        row={'X':xx,'full_analytic':mp.nstr(full,60),
             'absolute_errors':{str(n):mp.nstr(e,15) for n,e in errors.items()},
             'error_divided_by_X5_abs_logX':mp.nstr(scaled,15)}
        rows.append(row);print(row,flush=True)
    write_artifact(folder/'total_energy_low_check.json' , json.dumps({'rows':rows},indent=2)+'\n')
