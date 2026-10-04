"""Direct quadrature of the independently derived one-loop parameter kernels.

m=1, mu=1, and the common alpha/(4*pi) is factored out. A finite photon
mass above the physical cut makes this a smooth real quadrature test.
This checks assembly and Ward identities, not the inclusive physical total.
"""
from compton_artifacts import write_artifact
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from check_compton_amplitudes import I, gamma, metric, slash, prop, pol, dot, sandwich, adj


def sv(p):return np.einsum('...a,a,aij->...ij',p,metric,np.array(gamma))
def vv(p):return np.sum(p*p*metric,axis=-1)
def contract(A):return sum(metric[j]*gamma[j]@A@gamma[j] for j in range(4))

class Loop:
    def __init__(self,n,lam):
        self.lam=lam
        a,w=leggauss(n);self.x=(a+1)/2;self.w=w/2
        v=self.x;D=v*v+(1-v)*lam*lam;L=-np.log(D)
        self.cm=np.sum(self.w*(2*(1+v)*L-2*v))
        self.cz=np.sum(self.w*(-2*(1-v)*L+2*(1-v)+4*v*(1-v)*(1+v)/D))
    def sigma(self,P):
        v=self.x;D=v+(1-v)*self.lam**2-v*(1-v)*dot(P,P);L=-np.log(D)
        a=np.sum(self.w*(-2*(1-v)*L+2*(1-v)))
        b=np.sum(self.w*(4*L-2))
        return a*slash(P)+b*I-self.cm*I-(slash(P)-I)*self.cz
    def vertex(self,Pp,P,E):
        u,v=np.meshgrid(self.x,self.x,indexing='ij');wu,wv=np.meshgrid(self.w,self.w,indexing='ij')
        x=u.ravel();y=((1-u)*v).ravel();z=1-x-y;weight=(wu*wv*(1-u)).ravel()
        K=x[:,None]*Pp+y[:,None]*P;q=Pp-P
        D=(x+y)+z*self.lam**2-x*z*dot(Pp,Pp)-y*z*dot(P,P)-x*y*dot(q,q)
        N=contract((sv(Pp-K)+I)@E@(sv(P-K)+I))
        return np.einsum('a,aij->ij',weight,-N/D[:,None,None]+(-2*np.log(D)-4)[:,None,None]*E)+self.cz*E
    def spanning(self,p,pp,Q,A,B):
        u,v,w=np.meshgrid(self.x,self.x,self.x,indexing='ij')
        wu,wv,ww=np.meshgrid(self.w,self.w,self.w,indexing='ij')
        a1=u.ravel();a2=((1-u)*v).ravel();a3=((1-u)*(1-v)*w).ravel();a0=1-a1-a2-a3
        weight=(wu*wv*ww*(1-u)**2*(1-v)).ravel()
        K=a1[:,None]*p+a2[:,None]*Q+a3[:,None]*pp
        D=vv(K)+a0*self.lam**2+a2*(1-dot(Q,Q))
        if np.min(D)<=0:raise ValueError('This real quadrature test must stay below the photon-mass cut')
        A1=sv(p-K)+I;A2=sv(Q-K)+I;A3=sv(pp-K)+I
        N=contract(A3@B@A2@A@A1)
        box=0*N
        for j in range(4):
            g=gamma[j]
            box+=2*metric[j]*(g@B@g@A@A1+g@B@A2@A@g+A3@B@g@A@g)
        box=contract(box)
        return np.einsum('a,aij->ij',weight,N/D[:,None,None]**2-box/(4*D[:,None,None]))
    def ordered(self,p,pp,Q,A,B):
        S=prop(Q)
        return (self.vertex(pp,Q,B)@S@A+B@S@self.vertex(Q,p,A)
                +B@S@self.sigma(Q)@S@A+self.spanning(p,pp,Q,A,B))
    def full(self,p,k,kp,e,ep):
        pp=p+k-kp;A=slash(e);B=slash(ep)
        return self.ordered(p,pp,p+k,A,B)+self.ordered(p,pp,p-kp,B,A)

def main():
    omega=.3;c=.2;wp=omega/(1+omega*(1-c))
    p=np.array([1.,0,0,0]);k=np.array([omega,0,0,omega])
    kp=wp*np.array([1.,np.sqrt(1-c*c),0,c]);pp=p+k-kp
    e=pol(k)[0];ep=pol(kp)[0]
    rows=[];previous=None
    for n in [12,20,28]:
        loop=Loop(n,.7)
        R=loop.full(p,k,kp,e,ep)
        scale=max(1,sandwich(R,p,pp))
        ward=max(sandwich(loop.full(p,k,kp,k,ep),p,pp),
                 sandwich(loop.full(p,k,kp,e,kp),p,pp))/scale
        relative=None if previous is None else np.linalg.norm(R-previous)/np.linalg.norm(R)
        rows.append(dict(gauss_order=n,ward_residual=ward,change_from_previous=relative))
        previous=R
    assert rows[-1]['ward_residual']<1e-10
    assert rows[-1]['change_from_previous']<1e-10
    report=dict(omega_over_m=omega,cos_theta=c,lambda_over_m=.7,quadrature=rows)
    print(json.dumps(report,indent=2))
    write_artifact(Path(__file__).resolve().parents[1]/'virtual_checks.json' , json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
