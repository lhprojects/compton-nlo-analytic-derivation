"""Check the exact rational trace coefficients against explicit Dirac matrices."""
from compton_artifacts import read_artifact
import json
from pathlib import Path
import numpy as np
import sympy as sy
from compton_physical_virtual import PhysicalLoop,cm_event,moments
from compton_vertex_analytic import vertex_moments
from check_compton_amplitudes import slash,I,born,adj,prop,metric

records=json.loads(read_artifact(Path(__file__).with_name('compton_virtual_trace_analytic.json')))
x,y=sy.symbols('x y')
functions={name:{key:sy.lambdify((x,y),sy.sympify(value),'numpy') for key,value in data.items()}
           for name,data in records.items()}


def evaluate(s,c,lam):
    x=s-1;y=x-x*x/(2*s)*(1-c)
    params={'cz':-4-4*np.log(lam),'constant':1}
    for channel,S in [('s',s),('u',1-y)]:
        for suffix,j in zip(['00','10','11','20','21','22'],vertex_moments(S)):
            params[channel+'J'+suffix]=j
        F,_=moments(np.array([S]),np.array([1-S]))
        params[channel+'a']=-.5-2*S*F[1][0]+S*F[2][0]
        params[channel+'b']=2+4*S*F[1][0]
    return sum(f(x,y)*params[key] for data in functions.values() for key,f in data.items())


def matrix(s,c,lam):
    loop=PhysicalLoop(lam=lam);p,k,kp,pp=cm_event(s,c);total=0j
    # This is a diagram subset. Use the same covariant polarization
    # contraction as the exact trace; physical sums agree only for the
    # full Ward-invariant amplitude including the spanning graphs.
    for i,e in enumerate(np.eye(4)):
        for j,ep in enumerate(np.eye(4)):
            A=slash(e);B=slash(ep);out=0*I
            for Q,a,b in [(p+k,A,B),(p-kp,B,A)]:
                S=prop(Q)
                out+=loop.vertex(pp,Q,b)@S@a+b@S@loop.vertex(Q,p,a)+b@S@loop.sigma(Q)@S@a
            total+=metric[i]*metric[j]*np.trace((slash(pp)+I)@out@(slash(p)+I)@adj(born(p,k,kp,e,ep)))/4
    return total


if __name__=='__main__':
    for s,c in [(1.6,.2),(10,-.4),(100,.7)]:
        a=evaluate(s,c,1e-6);b=matrix(s,c,1e-6)
        error=abs(a-b)/max(1,abs(a))
        print(s,c,a,b,error)
        assert error<1e-10
