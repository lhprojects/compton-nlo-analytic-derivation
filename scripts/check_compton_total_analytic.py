"""Recompute the analytic sectors and compare their sum with the collected form.

Run directly with Python; no preliminary generation command is required.
This checks analytic evaluation, not numerical phase-space integration.
"""

import argparse
import mpmath as mp

from compton_total_independent_closed import independent_total
from compton_total_energy_closed import closed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--s', default='1.6', help='s_phys / m^2 (greater than 1)')
    parser.add_argument('--dps', type=int, default=15)
    args = parser.parse_args()
    if args.dps < 15:
        parser.error('--dps must be at least 15')
    with mp.workdps(args.dps + 40):
        s = mp.mpf(args.s)
        if not mp.isfinite(s) or s <= 1:
            parser.error('--s must be finite and greater than 1')

        def progress(name, value, parts):
            print(f'{name}: {mp.nstr(value, args.dps)}', flush=True)

        parts = independent_total(s, dps=args.dps, progress=progress)
        collected = closed(s - 1, dps=args.dps + 10)
        difference = abs(parts['total'] - collected)
        tolerance = mp.mpf(10) ** (-args.dps) * max(1, abs(collected))
        print(f'Collected form: {mp.nstr(collected, args.dps)}', flush=True)
        print(f'Absolute difference: {mp.nstr(difference, 8)}', flush=True)
        print(f'Tolerance: {mp.nstr(tolerance, 8)}', flush=True)
        if not mp.isfinite(difference) or difference > tolerance:
            raise ArithmeticError('Sector sum and collected analytic total disagree')
        print('PASS: analytic sector assembly and infrared pole cancellation', flush=True)


if __name__ == '__main__':
    main()
