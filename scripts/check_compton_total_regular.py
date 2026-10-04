"""Numerical checks of the separately certified regular and real-axis forms."""
import mpmath as mp
from compton_total_regular import coefficients, certify_real_pair, closed_regular, paired_real
from compton_total_energy_closed import closed
from compton_root_classical import (all_roots, groot, zero_root_m1, m4_root_m1,
    root_m1_m1, root_m2_m1, root_0_m1, root_m1_0)


def main():
    coefficients()
    certify_real_pair()
    with mp.workdps(80):
        for value in ['0.0002', '0.02', '0.6', '0.99999999999999999999', '1',
                      '1.00000000000000000001', '9', '99', '2000', '200000']:
            x = mp.mpf(value)
            roots = all_roots(x)
            for word, fn in [(('r','-1'),groot), (('0','r','-1'),zero_root_m1),
                             (('-4','r','-1'),m4_root_m1),
                             (('r','-1','-1'),root_m1_m1),
                             (('r','-2','-1'),root_m2_m1),
                             (('r','0','-1'),root_0_m1),
                             (('r','-1','0'),root_m1_0)]:
                original = fn(x)
                assert abs(roots[word]-original) < mp.mpf('1e-70')*max(1,abs(original)), (value,word)
            # The eighth word is checked independently through its derivative.
            derivative = mp.diff(lambda t: all_roots(t)[('r','r','-1')], x)
            expected = groot(x)/mp.sqrt(x*(x+4))
            assert abs(derivative-expected) < mp.mpf('1e-65')*max(1,abs(expected)), value
            old, new = closed(x, dps=35), closed_regular(x, dps=35)
            error = abs(new-old)
            assert error < mp.mpf('1e-30')*max(1, abs(old)), (value, error)
            pair = (mp.re(mp.polylog(3,x*x))-mp.log(x)*mp.re(mp.polylog(2,x*x)))/2
            assert abs(pair-paired_real(x)) < mp.mpf('1e-65'), value
            print(value, 'total difference', mp.nstr(error, 8), flush=True)
        # These coefficients follow from recursive analytic kernel integration
        # in compton_total_energy_low.py, not from numerical fitting. The
        # remainder O(X**5 abs(log(X))) is negligible at these two points.
        for value in ['1e-18', '1e-24']:
            x = mp.mpf(value)
            m = mp.log(x)
            low = -(80*m-21)*x**2/45+(120*m-113)*x**3/90-(3125*m-1827)*x**4/1125
            result = closed_regular(value, dps=35)
            error = abs(result-low)/abs(low)
            assert mp.isfinite(result) and error < mp.mpf('1e-32'), (value, result, error)
            print(value, 'relative low-energy difference', mp.nstr(error, 8), flush=True)
    print('PASS: shared root words, regular total and real pair at ten energies, including X=1')
    print('PASS: low-energy cancellation at X=1e-18 and 1e-24 against the analytically derived series')


if __name__ == '__main__':
    main()
