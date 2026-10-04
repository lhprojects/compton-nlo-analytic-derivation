"""Assemble analytic real-sector poles and compare to the derived soft current."""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
import json,sys
import sympy as S
import mpmath as mp
from compton_real_soft_matching import H
from compton_soft_uv_check import uv_difference_closed
from compton_soft_uv_uw import uw_endpoint_closed
from compton_soft_analytic import ir_angular_closed
sys.set_int_max_str_digits(0)
mp.mp.dps=70
base=Path(__file__).parent
d=json.loads(read_artifact(base/'compton_real_soft_matching_integrated.json'))
b=json.loads(read_artifact(base/'compton_real_bose_subtraction.json'))
func={k:S.lambdify(H,S.sympify(d[k]['pole']),'mpmath') for k in ['rational','logs']}
func['bose']=S.lambdify(H,S.sympify(b['pole']),'mpmath')
records=[]
for h in ['.6','9','99']:
    hn=mp.mpf(h)
    parts={k:f(hn) for k,f in func.items()}
    parts['UV']=uv_difference_closed(h,'pole',dps=55)
    parts['2UW']=2*uw_endpoint_closed(h,'pole',dps=55)
    actual=sum(parts.values());expected=8*hn*ir_angular_closed(hn+1)
    residual=abs(actual-expected)/max(1,abs(expected))
    print('H',h,'sum',mp.nstr(actual,55),'expected',mp.nstr(expected,55),'relative',mp.nstr(residual,5),flush=True)
    assert residual<mp.mpf('1e-45')
    records.append(dict(H=h,parts={k:mp.nstr(v,60) for k,v in parts.items()},expected=mp.nstr(expected,60),relative_residual=mp.nstr(residual,8)))
    write_artifact(Path(__file__).resolve().parents[1]/'real_integrated_ir_check.json' , json.dumps(records,indent=2)+'\n')
