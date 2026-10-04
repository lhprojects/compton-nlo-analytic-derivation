"""Quadrature checks of the analytically integrated endpoint matching only."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S
import mpmath as mp
from compton_real_soft_matching import H,r
mp.mp.dps=70
raw=json.loads(read_artifact(Path(__file__).with_name('compton_real_soft_matching_integrated.json')))
report=[]
for name in ('rational','logs'):
    data=raw[name]
    f=S.lambdify((H,r),S.sympify(data['residue_sum']),'mpmath')
    A=S.lambdify(H,S.sympify(data['pole']),'mpmath')
    M=S.lambdify(H,S.sympify(data['log_moment']),'mpmath')
    for h in map(mp.mpf,['.6','9','99']):
        R=h/(h+1);fac=mp.log(h+1) if name=='logs' else 1
        aq=fac*mp.quad(lambda t:f(h,t),[0,R/2,R])
        mq=fac*mp.quad(lambda t:f(h,t)*mp.log(t),[0,R/2,R])
        residual=max(abs(A(h)-aq),abs(M(h)-mq))/max(1,abs(aq),abs(mq))
        assert residual<mp.mpf('1e-60'),residual
        item=dict(sector=name,H=str(h),pole=mp.nstr(A(h),55),
                  finite_shift=mp.nstr(A(h)*mp.log(h)-M(h),55),relative_residual=mp.nstr(residual,5))
        report.append(item);print(item,flush=True)
write_artifact(Path(__file__).resolve().parents[1]/'real_soft_matching_check.json' , json.dumps(report,indent=2)+'\n')
