"""Complete independently derived NLO coefficient in the collected energy basis."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
from collections import defaultdict
from functools import lru_cache
import sympy as S
import mpmath as mp
from compton_soft_analytic_gpl import evaluate_gpl_taylor
from compton_soft_uv_check import strip_trailing_zeros
from compton_root_classical import (groot,zero_root_m1,m4_root_m1,
    root_m1_m1,root_m2_m1,root_0_m1,root_m1_0)

folder=Path(__file__).resolve().parents[1]
X,Y=S.symbols('X Y')
@lru_cache(None)
def coefficients():
    data=json.loads(read_artifact(folder/'total_energy_collected.json'))
    return [(tuple(row['word']),S.lambdify((X,Y),S.sympify(row['coefficient']),'mpmath')) for row in data['terms']]

def closed(x,dps=50):
    """Return C1(X), where x=X=s_phys/m**2-1, with no Born term."""
    with mp.workdps(dps+40):
        x=mp.mpf(x)
        if x<=0:raise ValueError('X must be positive')
        y=mp.sqrt(x/(x+4));q=2/(x+2+mp.sqrt(x*(x+4)));l=mp.log(q)
        roots={('r','-1'):groot(x),('0','r','-1'):zero_root_m1(x),
               ('-4','r','-1'):m4_root_m1(x),('r','-1','-1'):root_m1_m1(x),
               ('r','-2','-1'):root_m2_m1(x),('r','0','-1'):root_0_m1(x),
               ('r','-1','0'):root_m1_0(x),
               ('r','r','-1'):-l**3/6-mp.polylog(3,-q**3)/9+mp.polylog(3,-q)+mp.pi**2*l/18+2*mp.zeta(3)/3}
        g10=mp.zeta(2) if x==1 else mp.log(x)*mp.log(abs(1-x))+mp.re(mp.polylog(2,x))
        special={('1','0'):g10,('0','1','0'):2*mp.re(mp.polylog(3,x))-mp.log(x)*mp.re(mp.polylog(2,x))}
        rational=defaultdict(mp.mpc);rest=mp.mpc(0)
        for word,fn in coefficients():
            if word==('0',) and x==1:
                # The pole multiplying log(X) is removable, with limit +1.
                rest+=1
                continue
            c=fn(x,y)
            if 'r' in word:rest+=c*roots[word]
            elif word in special:rest+=c*special[word]
            else:
                for w,v in strip_trailing_zeros(tuple(map(mp.mpf,word)),mp.log(x)).items():rational[w]+=c*v
        value,_,_=evaluate_gpl_taylor(rational,x)
        value+=rest
        assert abs(mp.im(value))<mp.mpf(10)**(-dps)*max(1,abs(value))
        return +mp.re(value)

if __name__=='__main__':
    mp.mp.dps=100;rows=[]
    for x,path in [('0.6','independent_closed_s1p6_70.json'),('9','independent_closed_s10_70.json'),('99','independent_closed_s100_75.json')]:
        old=json.loads(read_artifact(folder/path))
        ref=mp.mpf(old['components']['total'])
        val=closed(x,60);error=abs(val-ref)
        print(x,mp.nstr(val,65),'difference',mp.nstr(error,12),flush=True)
        assert error<mp.mpf('1e-55')*max(1,abs(ref)),(x,error)
        rows.append({'X':x,'new_analytic':mp.nstr(val,70),'original_report':path,'absolute_difference':mp.nstr(error,12)})
    delta=mp.mpf('1e-20');mid=closed(1,65)
    left,right=closed(1-delta,65),closed(1+delta,65)
    error=abs((left+right)/2-mid)
    assert error<mp.mpf('1e-38'),error
    report={'rows':rows,'X1_removable_limit':{
        'value':mp.nstr(mid,70),'offset':'1e-20',
        'symmetric_continuity_error':mp.nstr(error,12)}}
    write_artifact(folder/'total_energy_closed_check.json' , json.dumps(report,indent=2)+'\n')
    data=json.loads(read_artifact(folder/'total_energy_collected.json'))
    data['numerical_check']='passed against original independent analytic totals at X=0.6,9,99 and removable-limit continuity at X=1; see total_energy_closed_check.json'
    write_artifact(folder/'total_energy_collected.json' , json.dumps(data,indent=2)+'\n')
