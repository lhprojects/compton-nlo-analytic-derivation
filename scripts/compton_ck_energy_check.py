"""Cross-check new CK energy masters using quadrature only for verification."""
from compton_artifacts import write_artifact
import json
from pathlib import Path
import mpmath as mp
from compton_root_classical import (groot,zero_root_m1,m4_root_m1,
    root_m1_m1,root_m2_m1,root_0_m1)
from compton_soft_analytic_gpl import evaluate_gpl_taylor

def root_masters(x):
    x=mp.mpf(x);r=mp.sqrt(x*(x+4));g=groot(x)
    i7=3*(x+2)**2/(4*x*r)*g
    i2=(x+2)/r*(mp.mpf(3)/8*(zero_root_m1(x)-mp.log1p(x)*g
            +2*root_m1_m1(x)-m4_root_m1(x))
            +(root_m2_m1(x)+root_0_m1(x))/2)
    return i2,i7

if __name__=='__main__':
    mp.mp.dps=75;rows=[]
    for xx in ['0.001','0.6','9','100']:
        x=mp.mpf(xx);s=1+x;b=x/(x+2)
        a2,a7=root_masters(x)
        def kernel(u,kind):
            if u==1:return mp.mpf(0)
            k=mp.sqrt((1+x*u*u)/s);f=mp.atanh(u/k)
            return f/(k*(1-b*b*u*u))*(mp.atanh(b*u) if kind==2 else u)
        q2=mp.quad(lambda u:kernel(u,2),[0,mp.mpf('.5'),1])
        q7=mp.quad(lambda u:kernel(u,7),[0,mp.mpf('.5'),1])
        expr={tuple(map(mp.mpf,w)):mp.mpf(c)/4 for w,c in
              [((-1,0,-1),1),((-1,-1,-1),-1),((-1,-2,-1),1)]}
        a1,_,_=evaluate_gpl_taylor(expr,x)
        q1=mp.quad(lambda u:b*mp.atanh(b*u)*mp.atanh(u)/(1-b*b*u*u),[0,mp.mpf('.5'),1])
        assert max(abs(a2-q2),abs(a7-q7),abs(a1-q1))<mp.mpf('1e-65')*max(1,abs(q2),abs(q7))
        row={'X':xx,'I2':mp.nstr(a2,55),'I7':mp.nstr(a7,55),
             'I2_quadrature_difference':mp.nstr(abs(a2-q2),10),
             'I7_quadrature_difference':mp.nstr(abs(a7-q7),10),
             'beta_I1_quadrature_difference':mp.nstr(abs(a1-q1),10)}
        rows.append(row);print(row,flush=True)
    folder=Path(__file__).resolve().parents[1]
    write_artifact(folder/'ck_energy_masters_check.json' , json.dumps({'rows':rows},indent=2)+'\n')
