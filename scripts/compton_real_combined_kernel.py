"""Collect the four real kernels before outer integration, with exact pole tests.

Only rational kernels obtained from the inner integrals are read. Neither
GPL primitives nor a published total enter this calculation.
"""
import json
from pathlib import Path
import sympy as S
from compton_artifacts import read_artifact, write_artifact

BASE = Path(__file__).resolve().parent
H, u, ell = S.symbols('H u ell')
s = H + 1
Q = u**2 + 2*u + s


def load(name):
    return json.loads(read_artifact(BASE / name))


def pole_order(expr):
    """Multiplicity of the irreducible polynomial Q in the reduced denominator."""
    denominator = S.Poly(S.cancel(expr).as_numer_denom()[1], u, domain=S.QQ.frac_field(H, ell))
    divisor = S.Poly(Q, u, domain=denominator.domain)
    order = 0
    while denominator.degree() >= divisor.degree():
        quotient, remainder = denominator.div(divisor)
        if not remainder.is_zero:
            break
        order += 1
        denominator = quotient
    return order


def build():
    rational = load('compton_real_analytic_rational_final_kernels.json')
    bose = load('compton_box_analytic_bose_kernels.json')
    logs = load('compton_real_analytic_logs_inner.json')
    single = load('compton_real_analytic_single_sector.json')
    kernels, orders = {}, {}
    for key in ('rational', 'L1', 'L2', 'L3'):
        a, b = S.sympify(rational[key]), S.sympify(bose[key])
        # Reduce this sum before introducing any quadratic root or GPL word.
        combined = S.cancel(a+b)
        assert pole_order(combined) == 0, key
        extra = ell*S.sympify(logs[key]) if key != 'L3' else S.S.Zero
        if key == 'L1':
            extra += S.sympify(single['u_rational_single_log'])
        assert pole_order(extra) == 0, key
        kernels[key] = S.factor(S.cancel(combined+extra))
        assert S.cancel(kernels[key]-a-b-extra) == 0, key
        assert pole_order(kernels[key]) == 0, key
        orders[key] = [pole_order(a), pole_order(b), pole_order(combined)]
    kernels['L1L2'] = S.factor(S.sympify(single['u_rational_two_logs']))
    assert pole_order(kernels['L1L2']) == 0
    assert any(a > 0 and b > 0 for a, b, _ in orders.values())
    return kernels, orders


def main():
    kernels, orders = build()
    report = {
        'source': 'independent rational inner-integral kernels R0, D, logs and W',
        'variables': {'H': 's-1', 'ell': 'ln(s)', 'domain': 's>1, sqrt(s)<u<s'},
        'kernel': 'rational + L1*L1_coefficient + L2*L2_coefficient + L3*L3_coefficient + L1*L2*L1L2_coefficient',
        'logarithms': {'L1': 'ln(u**2/s)', 'L2': 'ln(u*(u-1)/(s-u))',
                       'L3': 'ln((u-1)*(u+s)**2/(u*(u+1)**2*(s-u)))'},
        'coefficients': {key: str(value) for key, value in kernels.items()},
        'eliminated_factor': str(Q),
        'pole_orders_R0_D_sum': orders,
        'certificate': 'Exact rational reconstruction and polynomial division over Q(H,ell)[u]',
        'scope': 'Common finite-cutoff integrand. The original endpoint conversion remains necessary.'}
    write_artifact(BASE/'real_combined_kernel.json', json.dumps(report, indent=2)+'\n')
    print('Exact kernel reconstruction and quadratic-factor cancellation verified.')
    print('Pole orders in R0, D and their sum:', orders)


if __name__ == '__main__':
    main()
