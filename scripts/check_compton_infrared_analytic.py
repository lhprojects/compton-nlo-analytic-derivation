"""Check the infrared coefficient of the complete closed virtual amplitude."""
from compton_artifacts import read_artifact
import math
import numpy as np
import sympy as sy
import json
from pathlib import Path
from compton_amplitude_analytic import AnalyticAmplitude
from compton_physical_virtual import cm_event
from check_compton_amplitudes import slash,I,pol,born,adj,dot


def coefficients(s,c,lam):
    loop=AnalyticAmplitude(lam);p,k,kp,pp=cm_event(s,c);T=0j;B=0.
    for e in pol(k):
        for ep in pol(kp):
            a=loop.full(p,k,kp,e,ep);b=born(p,k,kp,e,ep)
            T+=np.trace((slash(pp)+I)@a@(slash(p)+I)@adj(b))/4
            B+=np.trace((slash(pp)+I)@b@(slash(p)+I)@adj(b)).real/4
    ze=dot(p,pp);be=ze*math.acosh(ze)/math.sqrt(ze*ze-1)
    return T.real/4,B*(be-1)


if __name__=='__main__':
    x,y=sy.symbols('x y')
    data=json.loads(read_artifact(Path(__file__).with_name('compton_box_analytic_trace.json')))
    coefficient=lambda ch,key:sy.sympify(data[ch][key])
    born_square=2*(x/y+y/x-4*(1/y-1/x)+4*(1/y-1/x)**2)
    residual=sy.cancel(coefficient('s','Ct')+coefficient('s','Db')/x
                       +coefficient('u','Ct')-coefficient('u','Db')/y
                       -4*born_square*(1+(x-y)/2))
    assert residual==0
    print('Exact box infrared identity residual:',residual)
    for s,c in [(1.6,.2),(10.,-.4),(100.,.7)]:
        v1,A=coefficients(s,c,1e-5);v2,_=coefficients(s,c,1e-8)
        residue=(v1-v2)/math.log(1e3)-A
        print('s=',s,'c=',c,'coefficient=',A,'residual=',residue)
        assert abs(residue)<1e-10*max(1,abs(A))
