"""Compare the collected classical finite real sum with the original GPL primitives."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S
import mpmath as mp
from compton_real_analytic_sector_eval import real_sector_zcut
from compton_box_analytic_bose_eval import bose_zcut
folder=Path(__file__).resolve().parents[1]
d=json.loads(read_artifact(folder/'real_four_finite_collected.json'));r=S.symbols('r')
modules=[{'L2':lambda z:mp.re(mp.polylog(2,z)),'Li3':lambda z:mp.polylog(3,z),'logabs':lambda z:mp.log(abs(z))},'mpmath']
atoms=[S.lambdify(r,S.sympify(a),modules) for a in d['atoms']]
coeff=[S.lambdify(r,S.sympify(t['coefficient']),'mpmath') for t in d['terms']]
def closed(s):
    rr=mp.sqrt(s);values=[f(rr) for f in atoms]
    return mp.fsum(c(rr)*mp.fprod(v**p for v,p in zip(values,t['powers']) if p) for c,t in zip(coeff,d['terms']))
if __name__=='__main__':
    mp.mp.dps=75;rows=[]
    path=folder/'real_four_finite_check.json'
    for value in ['1.6','10']:
        ss=mp.mpf(value);h=ss-1;parts={}
        for sector in ['rational','W','logs']:
            parts[sector]=real_sector_zcut(h,sector,dps=30)['finite_zcut']
            print(value,sector,mp.nstr(parts[sector],20),flush=True)
        parts['D']=bose_zcut(h,dps=30)['finite_zcut']
        original=sum(parts.values());short=closed(ss);error=abs(short-original)
        assert error<mp.mpf('1e-24')*max(1,abs(original)),(value,error,short,original)
        row={'s':value,'collected':mp.nstr(short,60),'original_four_sum':mp.nstr(original,45),'absolute_difference':mp.nstr(error,10)}
        rows.append(row);write_artifact(path , json.dumps({'original_requested_dps':30,'collected_working_dps':75,'source':'fresh evaluations of original GPL primitives, not reference values','rows':rows},indent=2)+'\n')
        print(row,flush=True)
