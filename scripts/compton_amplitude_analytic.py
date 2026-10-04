"""Complete one-loop Compton amplitude with no numerical loop integration.

Uses the independently derived self-energy, vertex and tensor-box closed
expressions. Numerical phase-space integration is not performed here.
"""
import numpy as np
from compton_physical_virtual import PhysicalLoop,cm_event
from compton_box_analytic import AnalyticBox
from check_compton_amplitudes import pol,sandwich


class AnalyticAmplitude(PhysicalLoop):
    def __init__(self,lam=1e-5):
        self.lam=lam
        self.cz=-4-4*np.log(lam)
        self.boxes={}

    def spanning(self,p,pp,Q,A,B):
        key=tuple(np.concatenate([p,pp,Q]))
        if key not in self.boxes:
            self.boxes[key]=AnalyticBox(p,pp,Q,self.lam)
        return self.boxes[key].matrix(A,B)


def check():
    for s,c in [(1.6,.2),(10.,-.4),(100.,.7)]:
        for lam in [1e-5,1e-8]:
            loop=AnalyticAmplitude(lam);p,k,kp,pp=cm_event(s,c)
            norm=max(1,sandwich(loop.full(p,k,kp,pol(k)[0],pol(kp)[0]),p,pp))
            residual=max(sandwich(loop.full(p,k,kp,k,pol(kp)[0]),p,pp),
                         sandwich(loop.full(p,k,kp,pol(k)[0],kp),p,pp))/norm
            print('s=',s,'c=',c,'lambda=',lam,'Ward residual=',residual,flush=True)
            assert residual<1e-9

if __name__=='__main__':check()
