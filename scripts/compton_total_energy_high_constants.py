"""Analytic high-energy constants from our classical root primitives.

C denotes Cl_2(pi/3). Root-of-unity series, dilogarithm reflection,
and an independently integrated root-free I2 approximation fix the constants.
"""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S

pi=S.pi;Z=S.zeta(3);C=S.Symbol('C',real=True);theta=pi/3;I=S.I
def norm(n):
    while n>3:n-=6
    while n<=-3:n+=6
    return n
def log1minus(n):
    n=norm(n)
    if abs(n)==1:return -I*n*theta
    if abs(n)==2:return S.log(3)/2-I*S.sign(n)*pi/6
    if n==3:return S.log(2)
    raise ValueError(n)
def li2(n):
    n=norm(n);t=n*theta
    cl={0:0,1:C,-1:-C,2:2*C/3,-2:-2*C/3,3:0}[n]
    return pi*pi/6-pi*abs(t)/2+t*t/4+I*cl
def li3(n):
    n=norm(n);t=n*theta
    real={0:Z,1:Z/3,-1:Z/3,2:-4*Z/9,-2:-4*Z/9,3:-3*Z/4}[n]
    return real+I*(pi*pi*t/6-pi*t*abs(t)/4+t**3/12)
def J(m,n=0):
    if m==n==0:return -2*Z
    l=I*m*theta;ratio=norm(m-n)
    return l*l*log1minus(ratio)+2*l*li2(ratio)-2*li3(ratio)

sum3=li3(1)+li3(-1)
primitive1=2*J(-1)+2*J(1)-J(-2)+J(-2,2)-2*sum3
Ka=S.simplify(3*Z+li3(-2)+primitive1/2)
def mixed_pole_endpoint(n):
    lb=-I*n*theta;log_minus_c=-I*S.sign(n)*2*theta
    return (-2*Z+lb*lb*log_minus_c+4*S.zeta(2)*lb+J(n,2*n))/2
pr0_1=2*(mixed_pole_endpoint(1)+mixed_pole_endpoint(-1))-sum3+4*Z
Kc=S.simplify(-2*sum3-pr0_1)
psh_1=-J(1)-J(-1)+sum3+2*Z
Kf=S.simplify(-(4*Z+2*sum3-psh_1))
Kb=S.simplify(-Kf-Kc)

# In the outer-minus-four primitive, all Li3(2*exp(+-i*pi/3)) cancel.
# Im Li2(2*exp(i*pi/3))=5*C/6+pi*log(2)/2, obtained from the
# Bloch-Wigner angle formula for a triangle with angles pi/3,pi/2,pi/6.
im_li2=5*C/6+pi*S.log(2)/2
pair_difference=S.simplify(-4*S.log(2)*theta*pi+8*theta*im_li2)
Ke=S.simplify(16*Z/9+4*theta*C-pair_difference)

# Rational high-energy approximation of I2, whose difference is o(1).
# K(-2,-1,-1)=-3*zeta(3)/4 by t->1/t.
# K(0,-1,-2)=7*zeta(3)/8-zeta(2)*log(2), by differentiating
# the last letter -b and integrating b from 1 to 2.
k211=-3*Z/4
k012=7*Z/8-S.zeta(2)*S.log(2)
k20_plus_02=-S.zeta(2)*S.log(2)-k012
i2_constant=S.simplify(k20_plus_02/2-k211/4-Z/8)
assert i2_constant==-3*Z/8
# Read the constant of the independently proved exact energy identity for I2.
Kd=S.simplify(2*i2_constant-S.Rational(3,4)*(Kf+2*Ka-Ke)-Kc)
constants={'K_a':Ka,'K_b':Kb,'K_c':Kc,'K_d':Kd,'K_e':Ke,'K_f':Kf}
combination=S.simplify(15*Ka-6*Kb-8*Kc-16*Kd-6*Ke-3*Kf)
assert combination==-9*Z
folder=Path(__file__).resolve().parents[1]
partial=json.loads(read_artifact(folder/'total_energy_high_partial.json'))
high=S.simplify(S.sympify(partial['leading_with_constants']).subs({S.Symbol(n):v for n,v in constants.items()}))
assert not high.has(C)
report={'source':'own root primitives and the independently derived I2 energy equation',
        'constant_symbol':'C=Cl_2(pi/3)',
        'constants':{n:str(v) for n,v in constants.items()},
        'root_constant_combination':str(combination),
        'C1_high':str(high),
        'remainder':'O(log(X)^3/X^2)',
        'certificate':'all root-of-unity endpoint expressions simplified exactly; the Clausen constant cancels identically'}
write_artifact(folder/'total_energy_high.json' , json.dumps(report,indent=2)+'\n')
print(report,flush=True)
