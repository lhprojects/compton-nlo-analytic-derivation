"""GPL shuffle algebra and causal energy-kernel pullbacks.

This module contains no published cross section or reference data.
"""
from pathlib import Path
from collections import defaultdict
from functools import lru_cache
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
q=S.Symbol('q')

@lru_cache(None)
def shuffle(a,b):
    if not a:return {b:1}
    if not b:return {a:1}
    out=defaultdict(int)
    for w,c in shuffle(a[1:],b).items():out[(a[0],)+w]+=c
    for w,c in shuffle(a,b[1:]).items():out[(b[0],)+w]+=c
    return dict(out)


def multiply(a,b):
    out=defaultdict(lambda:S.S.Zero)
    for u,c in a.items():
        for v,d in b.items():
            for w,n in shuffle(u,v).items():out[w]+=c*d*n
    return dict(out)


def kernel_map():
    omega=(1+S.sqrt(3)*S.I)/2
    golden=(3+S.sqrt(5))/2
    return {0:{S.Integer(1):2,S.Integer(0):-1},-1:{omega:1,S.conjugate(omega):1,S.Integer(0):-1},-2:{S.I:1,-S.I:1,S.Integer(0):-1},-4:{S.Integer(-1):2,S.Integer(0):-1},1:{golden:1,1/golden:1,S.Integer(0):-1},'r':{S.Integer(0):-1}}


def verify_kernels(maps):
    X=(1-q)**2/q
    for a,terms in maps.items():
        rhs=sum(c/(q-b) for b,c in terms.items())
        lhs=-1/q if a=='r' else S.diff(X,q)/(X-a)
        assert S.simplify(lhs-rhs)==0,a
    # Explicit real boundary logarithm, needed for trailing zero words.
    assert S.simplify(S.diff(2*S.log(1-q)-S.log(q),q)-S.diff(X,q)/X)==0
