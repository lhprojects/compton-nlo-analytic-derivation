"""Cross-check closed vertex moments and matrices against parameter integration."""
import mpmath as mp
import numpy as np
from compton_vertex_analytic import vertex_moments_mp
from compton_physical_virtual import PhysicalLoop,cm_event
from check_compton_amplitudes import slash,pol


def check():
    mp.mp.dps=60
    largest=mp.mpf(0)
    for S in map(mp.mpf,['-4','0','0.49','1.6','10']):
        closed=vertex_moments_mp(S)
        for value,(n,k) in zip(closed,[(0,0),(1,0),(1,1),(2,0),(2,1),(2,2)]):
            def integrand(v):
                A=1+(S-1)*v;B=1-A
                if abs(A)<mp.mpf('0.01'):
                    F=mp.fsum((-A/B)**j/(n+j+1) for j in range(40))/B
                else:
                    F=-(mp.log(abs(B))-(mp.j*mp.pi if B<0 else 0))/A
                    for j in range(1,n+1):
                        F=(mp.mpf(1)/j-B*F)/A
                return v**k*F
            q=mp.quad(integrand,[0,mp.mpf('.5'),1])
            err=abs(q-value)
            largest=max(largest,err)
            assert err<mp.mpf('1e-50'),(S,n,k,err)
    print('Maximum scalar moment residual:',mp.nstr(largest,8))
    for s,c in [(1.6,.2),(10.,-.4),(100.,.7)]:
        p,k,kp,pp=cm_event(s,c)
        for order in [32,64,128]:
            loop=PhysicalLoop(order)
            errors=[]
            for Pp,P,eps in [(pp,p+k,pol(kp)[0]),(p+k,p,pol(k)[0]),
                              (pp,p-kp,pol(k)[0]),(p-kp,p,pol(kp)[0])]:
                a=loop.vertex(Pp,P,slash(eps))
                b=loop.vertex_quadrature(Pp,P,slash(eps))
                errors.append(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a)))
            print('s=',s,'order=',order,'maximum matrix residual=',max(errors))
        assert max(errors)<1e-9


if __name__=='__main__':
    check()
