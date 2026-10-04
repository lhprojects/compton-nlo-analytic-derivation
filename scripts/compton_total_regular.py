"""Exact regular presentation of the independently collected Compton total.

Coefficients are reconstructed algebraically. Real-axis polylogarithm
identities are certified by differential equations and boundary constants.
No numerical matching or published total is used to obtain this form.
"""
from functools import lru_cache
import json
from pathlib import Path
import sympy as S
from compton_artifacts import read_artifact

X, Y, L, M = S.symbols('X Y L M')
ROOT = Path(__file__).resolve().parents[1]


@lru_cache(None)
def coefficients():
    rows = json.loads(read_artifact(ROOT/'total_energy_compact.json'))['terms']
    result = {}
    for row in rows:
        word = tuple(row['word'])
        original = X**3*S.sympify(row['coefficient'])
        regular = original - (M/(X-1) if not word else 0)
        expression = S.S.Zero
        for (l, m, y), value in S.Poly(S.expand(regular), L, M, Y).terms():
            choices = [S.factor(value), S.apart(value, X)]
            shortest = min(choices, key=lambda a: (S.count_ops(a), len(str(a))))
            expression += shortest*L**l*M**m*Y**y
        assert S.cancel(expression-regular) == 0, word
        assert S.cancel(expression).as_numer_denom()[1].subs(X, 1) != 0, word
        result[word] = expression
    assert S.cancel(result[('0', '-1', '0')]-result[('0', '1', '0')]) == 0
    return result


def certify_real_pair():
    # For t=ln X>0, D=Li2(exp(-2t)), H=Li3(exp(-2t)),
    # J=ln(1-exp(-2t)). Their derivatives are D'=2J and H'=-2D.
    t, D, H, J = S.symbols('t D H J')
    I2 = S.pi**2/3-2*t**2-D
    I3 = H-4*t**3/3+2*S.pi**2*t/3
    def derivative(expr):
        return S.diff(expr, t)+2*J*S.diff(expr, D)-2*D*S.diff(expr, H)
    # ln(exp(2t)-1)=2t+J follows from positive factors for t>0.
    assert S.expand(derivative(I2)+2*(2*t+J)) == 0
    assert S.expand(derivative(I3)-2*I2) == 0
    boundary = {t: 0, D: S.pi**2/6, H: S.zeta(3)}
    assert S.simplify(I2.subs(boundary)-S.pi**2/6) == 0
    assert S.simplify(I3.subs(boundary)-S.zeta(3)) == 0
    pair = (H+t*D)/2+t**3/3+S.pi**2*t/6
    assert S.expand((I3-t*I2)/2-pair) == 0
    assert S.simplify(pair.subs(boundary)-S.zeta(3)/2) == 0
    h = S.symbols('h')
    assert S.series(S.log(1+h)/h, h, 0, 3).removeO() == 1-h/2+h**2/3
    return True


def paired_real(x):
    """Evaluate the certified pair with real polylogarithm arguments in (0,1]."""
    import mpmath as mp
    x = mp.mpf(x)
    if x <= 0:
        raise ValueError('X must be positive')
    m = mp.log(x)
    if x <= 1:
        return (mp.polylog(3, x*x)-m*mp.polylog(2, x*x))/2
    return (mp.polylog(3, 1/(x*x))+m*mp.polylog(2, 1/(x*x)))/2+m**3/3+mp.pi**2*m/6


@lru_cache(None)
def numerical_coefficients():
    result = dict(coefficients())
    result.pop(('0', '1', '0'))
    return [(word, S.lambdify((X, Y, L, M), value, 'mpmath'))
            for word, value in result.items()]


def _closed_regular_at_precision(x):
    """Evaluate the unchanged exact formula at the active working precision."""
    import mpmath as mp
    from collections import defaultdict
    from compton_soft_analytic_gpl import evaluate_gpl_taylor
    from compton_soft_uv_check import strip_trailing_zeros
    from compton_root_classical import all_roots
    x = mp.mpf(x)
    y = mp.sqrt(x/(x+4))
    m = mp.log(x)
    roots = all_roots(x)
    special = {('1','0'):mp.zeta(2)-mp.polylog(2,1-x),
               ('0','-1','0'):paired_real(x)}
    rest = mp.mpc(1 if x == 1 else m/(x-1))
    rational = defaultdict(mp.mpc)
    for word, fn in numerical_coefficients():
        value = fn(x, y, mp.log1p(x), m)
        if word in roots:
            rest += value*roots[word]
        elif word in special:
            rest += value*special[word]
        else:
            for w, c in strip_trailing_zeros(tuple(map(mp.mpf,word)), m).items():
                rational[w] += value*c
    value, _, _ = evaluate_gpl_taylor(rational, x)
    return (value+rest)/x**3


def closed_regular(x, dps=50):
    """Evaluate the exact form, checking relative precision by repetition.

    Endpoint subtraction followed by division by X**3 loses digits as X
    approaches zero. The root coefficient can scale as X**(-7/2), while
    the total starts at X**2 log(X). Six guard digits per small-X decade
    cover this algebraic amplification. Repetition checks numerical
    convergence, not the analytic coefficients, constants or branches.
    """
    import mpmath as mp
    if not isinstance(dps, int) or dps < 1:
        raise ValueError('dps must be a positive integer')
    with mp.workdps(dps+50):
        argument = mp.mpf(x)
        if not mp.isfinite(argument) or argument <= 0:
            raise ValueError('X must be finite and positive')
        guards = 50+int(mp.ceil(6*max(0, -mp.log10(argument))))
    previous = None
    for extra in (0, 30, 60):
        with mp.workdps(dps+guards+extra):
            total = _closed_regular_at_precision(x)
            if mp.isfinite(total) and previous is not None and mp.isfinite(previous):
                scale = max(abs(total), abs(previous))
                tolerance = mp.mpf(10)**(-dps)*scale
                if abs(total-previous) <= tolerance and abs(mp.im(total)) <= tolerance:
                    return +mp.re(total)
            previous = total
    raise ArithmeticError(f'regular total did not converge to {dps} digits at X={x}')


if __name__ == '__main__':
    table = coefficients()
    certify_real_pair()
    print(f'PASS: {len(table)} exact coefficient reconstructions, all regular at X=1')
    print('PASS: real-axis inversion derivatives and analytic boundary constants')
