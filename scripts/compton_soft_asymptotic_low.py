"""Uniform threshold expansion of our soft integral, with exact coefficients."""
from compton_artifacts import write_artifact
import sympy as S,json
from pathlib import Path
x,c,u=S.symbols('x c u');N=6
ser=lambda f:S.series(f,x,0,N).removeO().expand()
s=1+x;t=-x*x*(1-c)/(2*s);E2=(x+2)**2/(4*s);beta2=x*x/(x+2)**2
M=1-t*u*(1-u);v2=ser(1-M/E2)
radial=sum(v2**j/S.Integer(2*j+1) for j in range(3))
K=ser((1-t/2)*ser(radial*ser(1/M)))
K=ser(S.integrate(K,(u,0,1))-sum(ser(beta2**j)/S.Integer(2*j+1) for j in range(3)))
y=x+t;B=ser(2*(x/y+y/x-4*(1/y-1/x)+4*(1/y-1/x)**2))
d=-t/2;bm1=ser(2*d/3-d*d/5)
IR=ser(S.integrate(ser(B*bm1),(c,-1,1)))
CK=ser(S.integrate(ser(B*K),(c,-1,1)))
A=ser(IR/s);Ks=ser(CK/s)
print('K=',S.collect(K,x),flush=True)
print('CIR=',IR,'CK=',CK,'A=',A,'CK/s=',Ks,flush=True)
report=dict(variable='x=s-1 ->0+',K=str(S.collect(K,x)),C_IR=str(IR),C_K=str(CK),A=str(A),C_K_over_s=str(Ks),remainder='O(x^6), uniformly in physical angle and soft parameter')
write_artifact(Path(__file__).resolve().parents[1]/'soft_asymptotic_low.json' , json.dumps(report,indent=2)+'\n')
