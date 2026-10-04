"""Independent one-dimensional quadrature check of the closed Bose primitive.

Quadrature appears only in this cross-check, never in the production interface.
"""
import sympy as S
import mpmath as mp
import compton_box_analytic_bose_eval as B


def check(h='9',dps=75):
    with mp.workdps(dps):
        h=mp.mpf(h);s=h+1;r=mp.sqrt(s)
        kernels={k:S.lambdify((B.H,B.u),e,'mpmath') for k,e in B.kernel_expressions('bose')}
        def f(u):
            return (kernels['rational'](h,u)+kernels['L1'](h,u)*mp.log(u*u/s)
                    +kernels['L2'](h,u)*mp.log(u*(u-1)/(s-u))
                    +kernels['L3'](h,u)*mp.log((u-1)*(u+s)**2/(u*(u+1)**2*(s-u))))
        delta=(s-r)/400
        v,_=B.primitive_value(h,'bose',s-delta)
        v-=B.base_constant(h,'bose')
        points=[r+(s-r)*q for q in (0,mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.75'),mp.mpf('.95'))]+[s-delta]
        independent=mp.quad(f,points)
        return v-independent


if __name__=='__main__':
    difference=check()
    print('closed primitive minus quadrature:',mp.nstr(difference,30),flush=True)
    assert abs(difference)<mp.mpf('1e-55'), difference
    lo=B.bose_zcut('9',35);hi=B.bose_zcut('9',60)
    with mp.workdps(90):
        assert abs(lo['finite_zcut']-hi['finite_zcut'])<mp.mpf('1e-30')*max(1,abs(hi['finite_zcut']))
        print('35 vs 60 digits:',mp.nstr(lo['finite_zcut']-hi['finite_zcut'],30),flush=True)
        print('finite z-cut:',mp.nstr(hi['finite_zcut'],65),flush=True)
