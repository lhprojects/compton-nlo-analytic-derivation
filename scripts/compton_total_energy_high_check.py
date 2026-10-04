"""Numerical checks of analytically derived high-energy constants and total."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S
import mpmath as mp
from compton_root_classical import (root_m1_m1,root_m1_0,root_0_m1,
    root_m2_m1,m4_root_m1,zero_root_m1)
from compton_total_energy_closed import closed

folder=Path(__file__).resolve().parents[1]
data=json.loads(read_artifact(folder/'total_energy_high.json'))
if __name__=='__main__':
    mp.mp.dps=110;x=mp.mpf('1e35');L=mp.log(x);C=mp.im(mp.polylog(2,mp.exp(mp.j*mp.pi/3)))
    values={'K_a':root_m1_m1(x)-L**3/6,
            'K_b':root_m1_0(x)-L**3/6+mp.pi**2*L/6,
            'K_c':root_0_m1(x)-L**3/6-mp.pi**2*L/6,
            'K_d':root_m2_m1(x)-L**3/6+mp.pi**2*L/12,
            'K_e':m4_root_m1(x)-L**3/6+mp.pi**2*L/18,
            'K_f':zero_root_m1(x)-L**3/6+mp.pi**2*L/18}
    checks={}
    for name,text in data['constants'].items():
        target=S.lambdify(S.Symbol('C'),S.sympify(text),'mpmath')(C)
        err=abs(values[name]-target)
        assert err<mp.mpf('1e-29'),(name,err)
        checks[name]={'exact_constant_value':mp.nstr(target,60),'finite_X_difference':mp.nstr(err,12)}
    rows=[];previous=mp.inf
    for xx in ['1000','100000','10000000']:
        x=mp.mpf(xx);L=mp.log(x);full=closed(x,65)
        approx=(L**3/3-L**2/2+17*L/4-mp.mpf(75)/8-mp.pi**2/2+4*mp.zeta(3))/x
        error=abs(full-approx);relative=error/abs(full)
        assert relative<previous;previous=relative
        row={'X':xx,'full_analytic':mp.nstr(full,65),'asymptotic':mp.nstr(approx,65),
             'relative_error':mp.nstr(relative,12),'error_divided_by_log3_over_X2':mp.nstr(error*x*x/L**3,12)}
        rows.append(row);print(row,flush=True)
    write_artifact(folder/'total_energy_high_check.json' , json.dumps({'constant_check_X':'1e35','constants':checks,'rows':rows},indent=2)+'\n')
