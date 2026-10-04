"""Numerical cross-check of the real q-pair endpoint reduction."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S
import mpmath as mp
folder=Path(__file__).resolve().parents[1]
d=json.loads(read_artifact(folder/'real_qpair_endpoint.json'));r=S.symbols('r')
expr=S.sympify(d['unit_interval_expression'])
f=S.lambdify(r,expr,modules=[{'L2':lambda z:mp.re(mp.polylog(2,z)),'logabs':lambda z:mp.log(abs(z))},'mpmath'])
mp.mp.dps=70;rows=[]
for value in ['1.000001','1.6','8.9','9','9.1','10','100']:
    s=mp.mpf(value);rr=mp.sqrt(s);length=s-rr
    def kernel(v):
        u=rr+length*v;delta=length*(1-v)
        Q=u*u-(s-3)*u+s;Qp=2*u-s+3
        C=(s-u*u)/(u*(u+1)**2*delta)
        return length*Qp*(-C if Q==0 else mp.log1p(-C*Q)/Q)
    points=[mp.mpf(0),mp.mpf(1)]
    if s>9:
        qp=(s-3+mp.sqrt((s-1)*(s-9)))/2
        if rr<qp<s:points.insert(1,(qp-rr)/length)
    numeric=mp.quad(kernel,points);closed=f(rr);error=abs(closed-numeric)
    assert error<mp.mpf('1e-55')*max(1,abs(closed)),(value,error)
    row={'s':value,'closed':mp.nstr(closed,60),'quadrature_difference':mp.nstr(error,8)}
    rows.append(row);print(row,flush=True)
write_artifact(folder/'real_qpair_endpoint_check.json' , json.dumps({'dps':70,'rows':rows,'source':'integral of Qprime/Q times log(A) with removable zeros treated analytically'},indent=2)+'\n')
