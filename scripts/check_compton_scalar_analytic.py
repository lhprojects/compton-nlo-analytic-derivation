"""Cross-check the independently derived scalar-box dilogarithm.

The production formula contains no numerical parameter integration.
Quadrature is used here only to check that formula against its defining
parameter integral. No published Compton total is imported.
"""
import mpmath as mp


def check(dps=50):
    mp.mp.dps = dps
    for transfer in ('0.001', '0.09', '5', '7800'):
        magnitude = mp.mpf(transfer)
        a = mp.sqrt(magnitude) / 2
        h = mp.sqrt(1 + a*a)
        kappa = a / h
        analytic = (mp.polylog(2, kappa) - mp.polylog(2, -kappa))/(2*a*h)

        def integrand(z):
            v = magnitude*z*(1-z)
            if not v:
                return mp.mpf(1)
            return mp.atan(mp.sqrt(v))/(mp.sqrt(v)*(1+v))

        quadrature = mp.quad(integrand, [0, mp.mpf('.5'), 1])
        residual = abs(analytic-quadrature)/max(1, abs(analytic))
        assert residual < mp.power(10, -dps+5)
        print(f'-t = {transfer}: ga = {mp.nstr(analytic, 35)}, '
              f'absolute residual = {mp.nstr(abs(analytic-quadrature), 5)}')


if __name__ == '__main__':
    check()
