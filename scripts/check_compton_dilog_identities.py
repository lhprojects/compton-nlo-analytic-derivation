"""Numerical implementation checks of identities already certified analytically."""
import mpmath as mp
import sympy as S
from compton_dilog_identities import certify, ordinary, quadratic, s, Ls, Lh, Lp, c, Li


def check():
    certify()
    mp.mp.dps = 70
    for expression, target in [(ordinary, lambda v: 1-v),
                                (quadratic, lambda v: -(v-1)**2/(4*v))]:
        evaluate = S.lambdify((s, Ls, Lh, Lp, c), expression,
                             [{'L2': lambda z: mp.polylog(2, z)}, 'mpmath'])
        for value in ['1.0001', '1.02', '1.6', '10', '1000']:
            v = mp.mpf(value)
            residual = evaluate(v, mp.log(v), mp.log(v-1), mp.log(v+1), mp.log(2))-mp.polylog(2, target(v))
            assert abs(residual) < mp.mpf('1e-60')
    print('PASS: real dilogarithm implementations at five points; constants fixed analytically')


if __name__ == '__main__':
    check()
