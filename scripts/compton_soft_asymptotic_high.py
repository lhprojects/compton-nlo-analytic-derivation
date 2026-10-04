"""Algebraic high-energy soft expansion from the original parameter integral."""
from compton_artifacts import write_artifact
import sympy as S,json
from pathlib import Path
r,w,L=S.symbols('r w L',positive=True);v=S.sqrt(1-r+r*w*w)
A=S.atanh(v)/v
# Exact identity underlying the finite shape derivative.
assert S.simplify(S.diff(w*A,w)-1/(1-w*w)-(1-r)*(A-1)/v**2)==0
F=(S.polylog(2,r)+S.log(r)*S.log(1-r))/2
assert S.simplify(S.expand_func(S.diff(F,r))+S.log(r)/(2*(1-r)))==0
# Elementary moments and the logarithmic moment zeta(3), obtained by the
# absolutely convergent expansion -log(1-r)=sum r^n/n.
z2=S.pi**2/6;z3=S.zeta(3)
C1=-z3
C2=z2/4-S.Rational(7,16)+(z2/4-S.Rational(5,16))+(S.Rational(1,2)-z2/4)
assert S.simplify(C2-(S.pi**2/24-S.Rational(1,4)))==0
CK=S.expand(4*((L*L/4-L/2)*(L+S.Rational(1,2))+S.pi**2*L/6+C1+C2))
IR=4*L*L-2*L-2*S.pi**2/3-5
report=dict(variable='x=s-1 -> infinity, L=ln x',A_leading_numerator=str(IR),C_K_over_s_leading_numerator=str(CK),A_remainder='O(L^2/x^2)',C_K_remainder='O(L^3/x^2)',finite_shape=str(F),shape_derivative='-ln(r)/(2*(1-r))',constant_moment=str(C1+C2),method='Exact parameter identity, integration by parts, convergent logarithmic series. No reference total or fit.')
write_artifact(Path(__file__).resolve().parents[1]/'soft_asymptotic_high.json' , json.dumps(report,indent=2)+'\n')
print(report)
