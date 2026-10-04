"""Independent Dirac-algebra checks of the Compton kernels, with m=e=1.

Tests physical spin/polarization sums, all real-photon Ward identities,
Bose symmetry, the soft limit, and the off-shell loop-vertex Ward identity.
Run with a Python environment containing NumPy.
"""
from compton_artifacts import write_artifact
import itertools
import json
from pathlib import Path
import numpy as np

I = np.eye(4, dtype=complex)
z = np.zeros((2,2), dtype=complex)
sigma = [np.array([[0,1],[1,0]]), np.array([[0,-1j],[1j,0]]), np.diag([1,-1])]
gamma = [np.diag([1,1,-1,-1])] + [np.block([[z,s],[-s,z]]) for s in sigma]
metric = np.array([1,-1,-1,-1])

def dot(a,b): return np.dot(a*metric,b)
def slash(a): return sum(a[j]*metric[j]*gamma[j] for j in range(4))
def prop(a): return (slash(a)+I)/(dot(a,a)-1)
def adj(a): return gamma[0]@a.conj().T@gamma[0]
def pol(k):
    n=k[1:]/np.linalg.norm(k[1:])
    ref=np.array([1.,0,0]) if abs(n[0])<0.8 else np.array([0.,1,0])
    a=np.cross(n,ref);a/=np.linalg.norm(a)
    return [np.r_[0.,a],np.r_[0.,np.cross(n,a)]]
def born(p,k,kp,e,ep):
    return slash(ep)@prop(p+k)@slash(e)+slash(e)@prop(p-kp)@slash(ep)
def real(p,k,k1,k2,e,e1,e2):
    ks=[k,-k1,-k2];es=list(map(slash,[e,e1,e2]));r=0*I
    for a,b,c in itertools.permutations(range(3)):
        r+=es[c]@prop(p+ks[a]+ks[b])@es[b]@prop(p+ks[a])@es[a]
    return r
def spin(K,p,pp): return float(np.trace((slash(pp)+I)@K@(slash(p)+I)@adj(K)).real)/4
def sandwich(K,p,pp): return np.linalg.norm((slash(pp)+I)@K@(slash(p)+I))
def final(omega,soft):
    p=np.array([1.,0,0,0]);k=np.array([omega,0,0,omega])
    n1=np.array([1.,0.6,0.,0.8]);n2=np.array([1.,-0.8,0.48,0.36])
    k1=soft*n1;A=p+k-k1
    w2=(dot(A,A)-1)/(2*dot(A,n2));k2=w2*n2
    return p,k,k1,k2,A-k2

def main():
    worst_born=worst_ward=worst_vertex=worst_bose=0.
    for w,c in [(0.1,-.4),(1.,.3),(20.,-.7)]:
        p=np.array([1.,0,0,0]);k=np.array([w,0,0,w]);wp=w/(1+w*(1-c))
        kp=wp*np.array([1.,np.sqrt(1-c*c),0,c]);pp=p+k-kp
        sq=sum(spin(born(p,k,kp,e,ep),p,pp) for e in pol(k) for ep in pol(kp))
        target=2*(w/wp+wp/w-(1-c*c))
        worst_born=max(worst_born,abs(sq/target-1))
        for e,ep in itertools.product(pol(k),pol(kp)):
            scale=sandwich(born(p,k,kp,e,ep),p,pp)
            for K in [born(p,k,kp,k,ep),born(p,k,kp,e,kp)]:
                worst_ward=max(worst_ward,sandwich(K,p,pp)/max(scale,1))
    for w in [.3,1.,7.]:
        p,k,k1,k2,pp=final(w,.04*w)
        assert abs(dot(pp,pp)-1)<1e-12
        for es in itertools.product(pol(k),pol(k1),pol(k2)):
            R=real(p,k,k1,k2,*es);scale=max(1,sandwich(R,p,pp))
            worst_bose=max(worst_bose,np.linalg.norm(R-real(p,k,k2,k1,es[0],es[2],es[1]))/max(1,np.linalg.norm(R)))
            for j,mom in enumerate([k,k1,k2]):
                replaced=list(es);replaced[j]=mom
                worst_ward=max(worst_ward,sandwich(real(p,k,k1,k2,*replaced),p,pp)/scale)
    P=np.array([2.,.3,-.4,.1]);Pp=np.array([1.7,-.2,.1,.6]);ell=np.array([.2,.7,.3,-.4])
    A=prop(Pp-ell);B=prop(P-ell)
    lhs=sum(metric[j]*gamma[j]@A@slash(Pp-P)@B@gamma[j] for j in range(4))
    rhs=sum(metric[j]*gamma[j]@(B-A)@gamma[j] for j in range(4))
    worst_vertex=np.linalg.norm(lhs-rhs)/np.linalg.norm(rhs)
    soft=[]
    for delta in [1e-3,1e-4,1e-5,1e-6]:
        p,k,k1,k2,pp=final(1.,delta)
        sq=sum(spin(real(p,k,k1,k2,*es),p,pp) for es in itertools.product(pol(k),pol(k1),pol(k2)))
        # Project the hard system to exact two-body kinematics at fixed direction.
        n=k2/k2[0];w2=1/(1+(1-n[3]));kb=w2*n;pb=p+k-kb
        b=sum(spin(born(p,k,kb,e,ep),p,pb) for e in pol(k) for ep in pol(kb))
        J=pb/dot(pb,k1)-p/dot(p,k1)
        ratio=sq/(-dot(J,J)*b)
        soft.append({'energy_over_mass':delta,'real_over_eikonal_born':ratio})
    assert worst_born<1e-11 and worst_ward<1e-10 and worst_vertex<1e-12 and worst_bose<1e-12
    assert abs(soft[-1]['real_over_eikonal_born']-1)<1e-5
    report=dict(born_relative_error=worst_born,ward_relative_residual=worst_ward,
                offshell_vertex_ward_residual=worst_vertex,bose_relative_residual=worst_bose,soft=soft)
    print(json.dumps(report,indent=2))
    out=Path(__file__).resolve().parents[1]/'amplitude_checks.json'
    out.parent.mkdir(parents=True,exist_ok=True);write_artifact(out , json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
