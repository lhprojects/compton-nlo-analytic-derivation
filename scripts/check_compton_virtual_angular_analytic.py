"""Numerical cross-checks of the independently integrated virtual subset."""
import mpmath as mp
import sympy as S
from compton_virtual_angular_analytic import crossed_integrand,x,y,cz,crossed_primitive_mp,coefficient_closed
from check_compton_virtual_trace_analytic import evaluate


def check():
    mp.mp.dps=50
    r=crossed_integrand();fun=S.lambdify((x,y,cz),r,'mpmath')
    for energy in ['1.6','10','100']:
        s=mp.mpf(energy);xx=s-1;zz=-4-4*mp.log(mp.mpf('1e-6'))
        def integrand(yy):
            if yy==1:
                return (6*zz*xx**2+6*zz*xx-18*zz-24*xx**2+4*mp.pi**2*xx**2-4*mp.pi**2*xx+15*xx+4*mp.pi**2+48)/(3*xx)
            a,b,c=fun(xx,yy,zz)
            return a+b*mp.log(yy)+c*(mp.pi**2/6-mp.polylog(2,1-yy))
        lo=xx/(1+xx);hi=xx
        points=[lo,hi] if hi<1 else [lo,mp.mpf('.999'),mp.mpf('1.001'),hi]
        q=mp.quad(integrand,points)
        closed=crossed_primitive_mp(xx,hi,zz)-crossed_primitive_mp(xx,lo,zz)
        print('s=',s,'crossed primitive residual=',mp.nstr(q-closed,5),flush=True)
        assert abs(q-closed)<mp.mpf('1e-40')
        # A separate double-precision trace quadrature verifies the overall
        # phase-space Jacobian and the direct-channel endpoint evaluation.
        from scipy.integrate import quad
        sf=float(s)
        check=quad(lambda c:evaluate(sf,c,1e-6).real/(4*sf),-1,1,
                   epsabs=1e-8,epsrel=1e-9)[0]
        result=coefficient_closed(s,'1e-6')
        print('subset coefficient=',mp.nstr(result,30),'trace quadrature=',check,flush=True)
        assert abs(float(result)-check)<1e-7*max(1,abs(float(result)))

if __name__=='__main__':check()
