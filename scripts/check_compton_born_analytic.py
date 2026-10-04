"""Exact Born normalization and angular endpoint checks, without sample values."""
import sympy as S
X, y = S.symbols('X y', positive=True)
Bamp = 2*(X/y+y/X-4*(1/y-1/X)+4*(1/y-1/X)**2)
P = y**2/(2*X)+(4/X+4/X**2)*y+(X-4-8/X)*S.log(y)-4/y
assert S.cancel(S.diff(P, y)-Bamp/2) == 0
eta = X/2
B = 2*((1+eta)/eta**2*(2*(1+eta)/(1+2*eta)-S.log(1+2*eta)/eta)
       + S.log(1+2*eta)/(2*eta)-(1+3*eta)/(1+2*eta)**2)
endpoint = 2*(P.subs(y, X)-P.subs(y, X/(1+X)))/X**2
assert S.cancel(S.expand_log(endpoint-B, force=False)) == 0
assert S.limit(B, X, 0) == S.Rational(8, 3)
print('PASS: Born primitive, physical endpoints, total normalization and Thomson limit')
