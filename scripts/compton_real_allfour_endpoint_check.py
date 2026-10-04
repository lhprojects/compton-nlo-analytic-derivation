"""Check the combined primitive's endpoint poles against original real kernels."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S
import mpmath as mp
from compton_real_analytic_sector_eval import end_series
from compton_box_analytic_bose_eval import end_series as bose_series
folder=Path(__file__).resolve().parents[1]
d=json.loads(read_artifact(folder/'real_allfour_endpoint.json'))
s,ell,k=S.symbols('s ell log_2sqrt_s_over_splus1')
fA=S.lambdify((s,ell,k),S.sympify(d['simple_pole_integrand_coefficient']),'mpmath')
fD=S.lambdify((s,ell),S.sympify(d['double_pole_integrand_coefficient']),'mpmath')
mp.mp.dps=70;rows=[]
for value in ['1.01','1.6','10']:
    ss=mp.mpf(value);h=ss-1
    original=[end_series(h,sector,2) for sector in ['rational','W','logs']]+[bose_series(h,'bose',2)]
    A=sum(v.get(-1,(0,0))[0] for v in original)
    D=sum(v.get(-2,(0,0))[0] for v in original)
    aa=fA(ss,mp.log(ss),mp.log(2*mp.sqrt(ss)/(ss+1)));dd=fD(ss,mp.log(ss))
    errors={'A':abs(A-aa),'D':abs(D-dd)}
    for name,error in errors.items():assert error<mp.mpf('1e-55')*max(1,abs(A if name=='A' else D)),(value,name,error)
    row={'s':value,'simple_pole':mp.nstr(aa,60),'double_pole':mp.nstr(dd,60),'differences':{n:mp.nstr(e,8) for n,e in errors.items()}}
    rows.append(row);print(row,flush=True)
write_artifact(folder/'real_allfour_endpoint_check.json' , json.dumps({'dps':70,'source':'independent Laurent expansions of original four real kernels','rows':rows},indent=2)+'\n')
