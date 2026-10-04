"""Exact certificates for the real dilogarithm substitutions, for s > 1.

No sampled values or published cross section enter these identities.
Positive logarithm factors determine the real branch throughout the domain.
"""
import sympy as S

s, Ls, Lh, Lp, c = S.symbols('s Ls Lh Lp c')
Li = S.Function('L2')
ordinary = Li(1/s) - Ls*Lh + Ls**2/2 - S.pi**2/6
elementary = S.expand(
    2*(Lh-Lp)*(c-Lp) - 2*(Lh-c-Ls)*(Lp-c-Ls)
    + (c+Ls-Lp)**2 - (2*Lp-2*c-Ls)**2/2
    - 2*c*(Lp-Ls) + c**2 - S.pi**2/6)
quadratic = 4*Li(2/(s+1)) - 2*Li(1/s) - 2*Li(1/(s+1)) + elementary


def certify():
    # These are exact factorizations of 1-z with positive factors for s>1.
    log_one_minus = {
        1/s: Lh-Ls,
        1/(s+1): Ls-Lp,
        2/(s+1): Lh-Lp,
    }

    def derivative(expr):
        atoms = sorted(expr.atoms(Li), key=str)
        symbols = S.symbols(f'dummy0:{len(atoms)}')
        polynomial = expr.xreplace(dict(zip(atoms, symbols)))
        result = (S.diff(polynomial, s) + S.diff(polynomial, Ls)/s
                  + S.diff(polynomial, Lh)/(s-1)
                  + S.diff(polynomial, Lp)/(s+1))
        for atom, symbol in zip(atoms, symbols):
            z = atom.args[0]
            result -= S.diff(polynomial, symbol)*log_one_minus[z]*S.diff(z, s)/z
        return result.xreplace(dict(zip(symbols, atoms)))

    targets = [1-s, -(s-1)**2/(4*s)]
    target_logs = [Ls, 2*Lp-2*c-Ls]
    certificates = []
    for name, expression, z, logarithm in zip(
            ['ordinary', 'quadratic'], [ordinary, quadratic], targets, target_logs):
        residual = S.factor(S.cancel(derivative(expression) + logarithm*S.diff(z, s)/z))
        if residual != 0:
            raise ArithmeticError(f'{name}: nonzero derivative residual {residual}')
        # Coefficients of log(s-1) vanish at least linearly. Hence their
        # product tends to zero as s->1+, rather than supplying a constant.
        log_coefficient = S.expand(expression).coeff(Lh)
        if S.simplify(log_coefficient.subs({Ls: 0, Lp: c})) != 0:
            raise ArithmeticError(f'{name}: unresolved endpoint logarithm')
        boundary = expression.subs({Ls: 0, Lp: c, Lh: 0}).subs(s, 1)
        boundary = boundary.subs({Li(1): S.pi**2/6,
                                 Li(S.Rational(1, 2)): S.pi**2/12-c**2/2})
        boundary = S.simplify(boundary)
        if boundary != 0:
            raise ArithmeticError(f'{name}: nonzero endpoint constant {boundary}')
        certificates.append({'identity': name, 'derivative_residual': '0',
                             'boundary_difference_at_s_1': '0'})
    return certificates


if __name__ == '__main__':
    for certificate in certify():
        print(certificate)
