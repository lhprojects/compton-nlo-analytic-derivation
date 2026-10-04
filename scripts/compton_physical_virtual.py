"""Independent physical Compton loop amplitude. Common alpha/(4*pi) removed.

Analytic radial parameter integration with the Feynman boundary value.
The scalar box is evaluated analytically in dilogarithms. QCDLoop
is optional and used only to cross-check that scalar reduction. All
numerator-dependent terms are derived and integrated here. The photon-mass
limit is taken in the finite terms.
"""
import ctypes
import math
from functools import lru_cache
from scipy.special import spence
import os
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from check_compton_amplitudes import I,gamma,metric,slash,prop,pol,dot,adj,born,sandwich
from check_compton_virtual import sv,contract,vv
from compton_vertex_analytic import vertex_moments

_lib=None
def scalar_box_qcdloop(S,t,lam):
    global _lib
    if _lib is None:
        _lib=ctypes.CDLL(os.environ.get('COMPTON_SCALAR_LIBRARY',str(Path(os.environ.get('XDG_CACHE_HOME',Path.home()/'.cache'))/'qftbook/compton-integration/compton_scalar_check.so')))
        _lib.compton_box.argtypes=[ctypes.c_double]*3+[ctypes.POINTER(ctypes.c_double)]
    r=(ctypes.c_double*6)();_lib.compton_box(S,t,-1.,r)
    if abs(r[4])+abs(r[5])>1e-20:
        raise ArithmeticError('Unexpected double infrared pole in massive-electron box')
    return complex(r[0],r[1])+2*np.log(lam)*complex(r[2],r[3])

@lru_cache(maxsize=4096)
def box_t_functions(t):
    T=-t
    if T<0:raise ValueError('Physical Compton momentum transfer must be nonpositive')
    g0=1. if T==0 else 4*math.asinh(math.sqrt(T)/2)/math.sqrt(T*(T+4))
    a=math.sqrt(T)/2
    h=math.sqrt(1+a*a)
    k=a/h
    if k<.2:
        # Li2(k)-Li2(-k)=2k sum_j k^(2j)/(2j+1)^2.
        ga=sum(k**(2*j)/(2*j+1)**2 for j in range(16))/(h*h)
    else:
        ga=(spence(1-k)-spence(1+k))/(2*a*h)
    return g0,ga

@lru_cache(maxsize=8192)
def scalar_box(S,t,lam):
    # Derived by isolating b=0, integrating rho and then b analytically.
    # D0 = [g0(t) log(lambda/(-(S-1)-i0)) + ga(t)]/(S-1).
    D=S-1;g0,ga=box_t_functions(t)
    log=math.log(lam/abs(D))+(1j*math.pi if D>0 else 0)
    return (g0*log+ga)/D


def moments(A,B,C=None):
    """F_n=integral rho^n/(A*rho+B-i0), R_n with squared denominator."""
    A,B=np.broadcast_arrays(np.asarray(A,dtype=complex),np.asarray(B,dtype=complex))
    C=A+B if C is None else np.asarray(C,dtype=complex)
    logB=np.log(B.astype(complex));logB=np.where(B.real<0,np.log(-B)-1j*np.pi,logB)
    F=[(np.log(C)-logB)/A]
    for n in (1,2):F.append((1/n-B*F[-1])/A)
    R=[1/(B*C)]
    for n in (1,2,3):R.append((F[n-1]-B*R[-1])/A)
    small=np.abs(A/B)<.15
    if np.any(small):
        z=-A[small]/B[small]
        for n in range(3):
            F[n]=np.array(F[n]);F[n][small]=sum(z**j/(n+j+1) for j in range(30))/B[small]
        for n in range(4):
            R[n]=np.array(R[n]);R[n][small]=sum((j+1)*z**j/(n+j+1) for j in range(30))/B[small]**2
    return F,R

class PhysicalLoop:
    def __init__(self,n=48,lam=1e-5):
        a,w=leggauss(n);self.x=(a+1)/2;self.w=w/2;self.lam=lam
        self.cz=-4-4*np.log(lam)
    def sigma(self,P):
        S=dot(P,P);F,_=moments(np.array([S]),np.array([1-S]));f1,f2=F[1][0],F[2][0]
        a=-.5-2*S*f1+S*f2;b=2+4*S*f1
        return a*slash(P)+(b-4)*I-(slash(P)-I)*self.cz
    def vertex(self,Pp,P,E):
        s1=dot(P,P);s2=dot(Pp,Pp)
        if abs(s1-1)<1e-7:
            on,off,S=P,Pp,s2
        elif abs(s2-1)<1e-7:
            on,off,S=Pp,P,s1
        else:
            raise ValueError('Closed vertex requires one on-shell electron')
        j00,j10,j11,j20,j21,j22=vertex_moments(S)
        H=slash(on);Z=slash(off-on)
        P1=slash(P)+I;P2=slash(Pp)+I
        N0=contract(P2@E@P1)
        U0=-contract(H@E@P1+P2@E@H)
        U1=-contract(Z@E@P1+P2@E@Z)
        V0=contract(H@E@H)
        V1=contract(H@E@Z+Z@E@H)
        V2=contract(Z@E@Z)
        return (-N0*j00-U0*j10-U1*j11-V0*j20-V1*j21-V2*j22
                +E*(-1.5+j20+(S-1)*j21+self.cz))
    def vertex_quadrature(self,Pp,P,E):
        """Independent parameter quadrature, for cross-checks only."""
        u=self.x
        den=u**3+(1-u)**3;v=u**3/den
        weight=self.w*3*u**2*(1-u)**2/den**2
        K=v[:,None]*Pp+(1-v[:,None])*P
        s1=dot(P,P);s2=dot(Pp,Pp)
        if abs(s1-1)<1e-7:s1=1.
        if abs(s2-1)<1e-7:s2=1.
        A=v*s2+(1-v)*s1;B=v*(1-s2)+(1-v)*(1-s1)
        F,_=moments(A,B,np.ones_like(A));K=sv(K);P1=slash(P)+I;P2=slash(Pp)+I
        N0=contract(P2@E@P1)
        N1=-contract(K@E@P1+P2@E@K)
        N2=contract(K@E@K)
        val=-F[0][:,None,None]*N0-F[1][:,None,None]*N1-F[2][:,None,None]*N2
        val+=(-1.5+A*F[2])[:,None,None]*E
        return np.einsum('a,aij->ij',weight,val)+self.cz*E
    def spanning(self,p,pp,Q,A,B):
        # b is the fraction of the internal electron, z partitions external lines.
        b,z=np.meshgrid(self.x,self.x,indexing='ij');wb,wz=np.meshgrid(self.w,self.w,indexing='ij')
        db=b**3+(1-b)**3;dz=z**3+(1-z)**3
        jb=3*b**2*(1-b)**2/db**2;jz=3*z**2*(1-z)**2/dz**2
        b=b**3/db;z=z**3/dz
        weight=(wb*wz*jb*jz*(1-b)).ravel();b=b.ravel();z=z.ravel()
        K=b[:,None]*Q+(1-b)[:,None]*(z[:,None]*p+(1-z)[:,None]*pp)
        cc=1-dot(p-pp,p-pp)*(1-b)**2*z*(1-z)
        bb=b*(1-dot(Q,Q));aa=cc-bb;F,R=moments(aa,bb,cc);K=sv(K)
        P1=slash(p)+I;P2=slash(Q)+I;P3=slash(pp)+I
        N0=contract(P3@B@P2@A@P1)
        N1=-contract(K@B@P2@A@P1+P3@B@K@A@P1+P3@B@P2@A@K)
        N2=contract(K@B@K@A@P1+K@B@P2@A@K+P3@B@K@A@K)
        N3=-contract(K@B@K@A@K)
        box0=0*I;box1=0*K
        for j in range(4):
            g=gamma[j]
            box0+=2*metric[j]*(g@B@g@A@P1+g@B@P2@A@g+P3@B@g@A@g)
            box1-=2*metric[j]*(g@B@g@A@K+g@B@K@A@g+K@B@g@A@g)
        box0=contract(box0);box1=contract(box1)
        val=N1*R[1][:,None,None]+N2*R[2][:,None,None]+N3*R[3][:,None,None]
        val-=.25*(box0*F[1][:,None,None]+box1*F[2][:,None,None])
        return N0*scalar_box(dot(Q,Q),dot(p-pp,p-pp),self.lam)+np.einsum('a,aij->ij',weight,val)
    def ordered(self,p,pp,Q,A,B):
        S=prop(Q)
        return (self.vertex(pp,Q,B)@S@A+B@S@self.vertex(Q,p,A)
          +B@S@self.sigma(Q)@S@A+self.spanning(p,pp,Q,A,B))
    def full(self,p,k,kp,e,ep):
        pp=p+k-kp;A=slash(e);B=slash(ep)
        return self.ordered(p,pp,p+k,A,B)+self.ordered(p,pp,p-kp,B,A)

def cm_event(s,c):
    E=(s+1)/(2*np.sqrt(s));q=(s-1)/(2*np.sqrt(s))
    p=np.array([E,0,0,-q]);k=np.array([q,0,0,q]);kp=q*np.array([1,np.sqrt(1-c*c),0,c])
    return p,k,kp,p+k-kp

def interference(s,c,n=48,lam=1e-5,ward=False):
    loop=PhysicalLoop(n,lam);p,k,kp,pp=cm_event(s,c)
    tr=0j;bt=0.
    for e in pol(k):
        for ep in pol(kp):
            R=loop.full(p,k,kp,e,ep);B=born(p,k,kp,e,ep)
            tr+=np.trace((slash(pp)+I)@R@(slash(p)+I)@adj(B))/4
            bt+=np.trace((slash(pp)+I)@B@(slash(p)+I)@adj(B)).real/4
    if ward:
        norm=max(1,sandwich(loop.full(p,k,kp,pol(k)[0],pol(kp)[0]),p,pp))
        print('ward',max(sandwich(loop.full(p,k,kp,k,pol(kp)[0]),p,pp),sandwich(loop.full(p,k,kp,pol(k)[0],kp),p,pp))/norm,flush=True)
    return tr,bt

if __name__=='__main__':
    for lam in [1e-3,1e-5,1e-7]:
        for n in [24,48,96]:
            tr,b=interference(1.6,.2,n,lam,True)
            print(lam,n,tr,b,flush=True)
