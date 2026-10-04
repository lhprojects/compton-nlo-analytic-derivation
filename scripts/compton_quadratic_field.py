"""Exact arithmetic in Q(r,y)[j,k]/(j^2-1+r^2,k^2-r^4+10r^2-9).

Reduce after every field operation to avoid expanding huge unreduced common
denominators. This is symbolic rational arithmetic, with no sampling.
"""
from functools import lru_cache
import sympy as S
from sympy.polys.fields import field

r, y, j, k = S.symbols('r y j k')
F, rf, yf = field('r,y', S.QQ)
ZERO = (F.zero,)*4
ONE = (F.one, F.zero, F.zero, F.zero)
fj, fk = 1-rf**2, rf**4-10*rf**2+9


def add(a, b):
    return tuple(x+z for x, z in zip(a, b))


def multiply(a, b):
    out = list(ZERO)
    for m, x in enumerate(a):
        if not x:
            continue
        for n, z in enumerate(b):
            if not z:
                continue
            value = x*z
            if (m & n) & 1:
                value *= fj
            if (m & n) & 2:
                value *= fk
            out[m ^ n] += value
    return tuple(out)


def inverse(a):
    if not any(a[1:]):
        return (1/a[0], F.zero, F.zero, F.zero)
    product = ONE
    for mask in (1, 2, 3):
        conjugate = tuple(-c if ((m & mask).bit_count() % 2) else c
                          for m, c in enumerate(a))
        product = multiply(product, conjugate)
    norm = multiply(a, product)
    assert not any(norm[1:]) and norm[0]
    return tuple(c/norm[0] for c in product)


@lru_cache(None)
def convert(expr):
    if not expr.has(j, k):
        return (F.from_expr(expr), F.zero, F.zero, F.zero)
    if expr == j:
        return (F.zero, F.one, F.zero, F.zero)
    if expr == k:
        return (F.zero, F.zero, F.one, F.zero)
    if expr.is_Add:
        out = ZERO
        for term in expr.args:
            out = add(out, convert(term))
        return out
    if expr.is_Mul:
        out = ONE
        for term in expr.args:
            out = multiply(out, convert(term))
        return out
    if expr.is_Pow and expr.exp.is_Integer:
        power = int(expr.exp)
        value = convert(expr.base)
        if power < 0:
            value, power = inverse(value), -power
        out = ONE
        while power:
            if power & 1:
                out = multiply(out, value)
            value = multiply(value, value)
            power >>= 1
        return out
    raise ValueError(expr)


def reduce_expression(expr):
    substitutions = {}
    for power in expr.atoms(S.Pow):
        if power.exp.is_Rational and power.exp.q == 2:
            base = S.expand(power.base)
            if base == 1-r*r:
                substitutions[power] = j**(2*power.exp)
            elif base == r**4-10*r*r+9:
                substitutions[power] = k**(2*power.exp)
            else:
                raise ValueError(('unexpected quadratic radical', power))
    return convert(expr.xreplace(substitutions))
