"""Independent analytic rank-three Compton box reduction.
Every scalar integral is elementary or a dilogarithm. No quadrature.
Dimensional UV residues are propagated to preserve rational finite terms.
"""
import itertools
import numpy as np
from scipy.special import spence
from check_compton_amplitudes import I,gamma,metric,slash,dot
from check_compton_virtual import contract

def li2(z): return spence(1-complex(z))
def pair(p,f): return np.array([p,f],dtype=complex)
def dimension(x,c): return np.array([x[0]/c,x[1]/c+2*x[0]/c**2])

class AnalyticBox:
    def __init__(self,p,pp,Q,lam):
        self.q=np.array([np.zeros(4),p,Q,pp]);self.mass2=np.array([0.,1.,1.,1.])
        self.S=dot(Q,Q);self.t=dot(p-pp,p-pp);self.lam=lam
        self.T=max(0.,-self.t);self.a=np.sqrt(self.T)/2;self.h=np.sqrt(1+self.a**2)
        self.r=np.arcsinh(self.a);self.g0=1. if self.a==0 else self.r/(self.a*self.h)
        self.cache={}
    def scalar(self,ids):
        n=len(ids)
        if n==1:
            m=self.mass2[ids[0]]
            return pair(m,0 if m==0 else m*(1-np.log(m)))
        if n==2:
            a,b=ids;s=dot(self.q[a]-self.q[b],self.q[a]-self.q[b])
            if a==0:
                if abs(s)<1e-10: val=1+s/2+s*s/6
                elif abs(s-1)<1e-10: val=2
                else: val=2+(1-s)/s*np.log(complex(1-s,-0.0))
            elif abs(s)<1e-10: val=s/6+s*s/60
            else:
                aa=np.sqrt(-s)/2
                val=2-2*np.sqrt(1+aa*aa)/aa*np.arcsinh(aa)
            return pair(1,val)
        if n==3:
            if ids==(1,2,3): val=-.5 if self.T==0 else -2*self.r**2/self.T
            elif ids==(0,1,3):
                if self.a<1e-5: hlog=self.T/6-self.T**2/20
                else: hlog=(2*self.r*np.log(2*self.h)-self.r**2-li2(-np.exp(-2*self.r))+li2(-1))/(self.a*self.h)
                val=self.g0*np.log(self.lam)-hlog/2
            else:
                val=-(li2(complex(self.S,1e-100 if self.S>1 else 0))-np.pi**2/6)/(self.S-1)
            return pair(0,val)
        if n==4:
            k=self.a/self.h;ga=1. if self.a==0 else (li2(k)-li2(-k))/(2*self.a*self.h)
            D=self.S-1;log=np.log(self.lam/abs(D))+(1j*np.pi if D>0 else 0)
            return pair(0,(self.g0*log+ga)/D)
        raise ValueError(ids)
    def shifted(self,ids,rank,v):
        out=np.zeros((2,)+(4,)*rank,dtype=complex)
        for bits in itertools.product((0,1),repeat=rank):
            slots=[j for j,b in enumerate(bits) if b];mom=self.tensor(ids,len(slots))
            for abc in itertools.product(range(4),repeat=rank):
                coefficient=np.prod([v[abc[j]] for j,b in enumerate(bits) if not b])
                out[(slice(None),)+abc]+=coefficient*mom[(slice(None),)+tuple(abc[j] for j in slots)]
        return out
    def tensor(self,ids,rank):
        key=(ids,rank)
        if key in self.cache:return self.cache[key]
        if rank==0:
            result=self.scalar(ids);self.cache[key]=result;return result
        n=len(ids);origin=self.q[ids[0]]
        if n==2 and rank==1:
            v=self.q[ids[1]]-origin;B=self.scalar(ids)
            if self.mass2[ids[0]]==self.mass2[ids[1]]: coefficient=B/2
            else:
                s=dot(v,v)
                if abs(s)<1e-10: coefficient=pair(.5,.25+s/6+s*s/24)
                else: coefficient=(self.scalar((ids[1],))-self.scalar((ids[0],))+(s+self.mass2[ids[0]]-self.mass2[ids[1]])*B)/(2*s)
            result=np.einsum('i,z->zi',origin,B)+np.einsum('i,z->zi',v,coefficient)
            self.cache[key]=result;return result
        basis=self.q[list(ids[1:])]-origin;cov=basis*metric
        inv=np.linalg.inv(cov@basis.T);dual=inv@basis
        projector=np.diag(metric)-np.einsum('ai,aj->ij',dual,basis)
        f=np.einsum('ai,ai->a',cov,basis)+self.mass2[ids[0]]-self.mass2[list(ids[1:])]
        previous=self.shifted(ids,rank-1,-origin);basepinch=self.shifted(ids[1:],rank-1,-origin)
        V=np.array([(basepinch-self.shifted(ids[:j]+ids[j+1:],rank-1,-origin)+f[j-1]*previous)/2 for j in range(1,n)])
        if rank==1: canonical=np.einsum('ai,az->zi',dual,V)
        elif rank==2:
            H=np.einsum('bj,azj->azb',cov,V)
            canonical=np.einsum('ai,bj,azb->zij',dual,dual,H)
            trace=self.shifted(ids[1:],0,-origin)+self.mass2[ids[0]]*self.scalar(ids)
            Z=dimension(trace-np.einsum('ab,azb->z',inv,H),4-len(basis))
            canonical+=np.einsum('ij,z->zij',projector,Z)
        elif rank==3:
            H=np.einsum('bj,ck,azjk->azbc',cov,cov,V)
            canonical=np.einsum('ai,bj,ck,azbc->zijk',dual,dual,dual,H)
            trace=self.shifted(ids[1:],1,-origin)+self.mass2[ids[0]]*self.shifted(ids,1,-origin)
            leftover=np.einsum('ai,zi->az',cov,trace)-np.einsum('bc,azbc->az',inv,H)
            Z=np.array([dimension(row,4-len(basis)) for row in leftover]);transverse=np.einsum('ai,az->zi',dual,Z)
            canonical+=np.einsum('ij,zk->zijk',projector,transverse)+np.einsum('ik,zj->zijk',projector,transverse)+np.einsum('jk,zi->zijk',projector,transverse)
        else:raise ValueError(rank)
        result=canonical.copy()
        for bits in itertools.product((0,1),repeat=rank):
            if all(bits):continue
            slots=[j for j,b in enumerate(bits) if b];mom=self.shifted(ids,len(slots),-origin)
            for abc in itertools.product(range(4),repeat=rank):
                coefficient=np.prod([origin[abc[j]] for j,b in enumerate(bits) if not b])
                result[(slice(None),)+abc]+=coefficient*mom[(slice(None),)+tuple(abc[j] for j in slots)]
        self.cache[key]=result;return result
    def matrix(self,A,B):
        p,Q,pp=self.q[1:];factors=[slash(pp)+I,slash(Q)+I,slash(p)+I];out=np.zeros((4,4),dtype=complex)
        for bits in itertools.product((0,1),repeat=3):
            rank=sum(bits);mom=self.tensor((0,1,2,3),rank)[1]
            for abc in itertools.product(range(4),repeat=rank):
                seq=iter(abc);mats=[]
                for bit,fixed in zip(bits,factors):
                    if bit:
                        j=next(seq);mats.append(-metric[j]*gamma[j])
                    else:mats.append(fixed)
                coefficient=mom[abc] if rank else mom
                out+=coefficient*contract(mats[0]@B@mats[1]@A@mats[2])
        return out
