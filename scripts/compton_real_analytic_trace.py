"""Exact six-diagram double-Compton trace, independently of total formulae.

The photon polarizations are summed using -g after summing all six orders.
The Clifford reduction uses only anticommutation and trace cyclicity.
Run with the project numerical Python (SymPy is required).
"""
from compton_artifacts import write_artifact
from functools import lru_cache
from itertools import permutations, product
from collections import defaultdict
from pathlib import Path
import json
import sympy as S


def canon(word):
    if not word: return ()
    w=tuple(word)
    return min(w[i:]+w[:i] for i in range(len(w)))


def add(out, terms, fac=1):
    for mon,c in terms: out[mon]+=fac*c


@lru_cache(None)
def trace(w):
    if not w: return (((),4),)
    if len(w)%2: return ()
    indices=set(i for i in w if i<0)
    if indices:
        choices=[]
        for j in indices:
            where=[i for i,x in enumerate(w) if x==j]
            for start,end in (where,where[::-1]):
                choices.append(((end-start-1)%len(w),start))
        n,start=min(choices)
        w=w[start:]+w[:start]
        middle=w[1:n+1];tail=w[n+2:]
        out=defaultdict(int)
        add(out,trace(canon(middle+tail)),4*(-1)**n)
        for i in range(n):
            moved=middle[:i]+middle[i+1:]+(middle[i],)+tail
            add(out,trace(canon(moved)),2*(-1)**i)
    else:
        out=defaultdict(int)
        for i in range(1,len(w)):
            pair=tuple(sorted((w[0],w[i])))
            for mon,c in trace(canon(w[1:i]+w[i+1:])):
                out[tuple(sorted(mon+(pair,)))]+=(-1)**(i+1)*c
    return tuple((mon,c) for mon,c in sorted(out.items()) if c)


# Momentum basis is (p,k,k1,k2), m=1.
vec=[(1,0,0,0),(1,1,-1,-1),(1,1,0,0),(1,0,-1,0),
     (1,0,0,-1),(1,0,-1,-1),(1,1,0,-1),(1,1,-1,0)]
initial=[2,3,4]
final=[5,6,7]
orders=list(permutations(range(3)))


def interference(i,j):
    a,b,c=orders[i];aa,bb,cc=orders[j]
    # Trace [(p'+1) E_c (q_c+1) E_b (q_a+1) E_a
    #        (p+1) E_aa (q_aa+1) E_bb (q_cc+1) E_cc].
    slots=[1,-c-1,final[c],-b-1,initial[a],-a-1,
           0,-aa-1,initial[aa],-bb-1,final[cc],-cc-1]
    positions=[k for k,v in enumerate(slots) if v>=0]
    out=defaultdict(int)
    for bits in product((0,1),repeat=6):
        if sum(bits)%2: continue
        keep=dict(zip(positions,bits))
        w=tuple(v for k,v in enumerate(slots) if v<0 or keep[k])
        add(out,trace(canon(w)),-1) # Three polarization sums -g.
    return [(mon,S.Rational(c,4)) for mon,c in out.items() if c]


def main():
    H,x,y,U,V=S.symbols('H x y U V')
    z=x+y-H
    gram=S.Matrix([[1,H/2,U/2,V/2],
       [H/2,0,(x-U)/2,(y-V)/2],
       [U/2,(x-U)/2,0,z/2],
       [V/2,(y-V)/2,z/2,0]])
    dots={(i,j):S.expand((S.Matrix(vec[i]).T*gram*S.Matrix(vec[j]))[0])
          for i in range(8) for j in range(i,8)}
    denominator=[]
    for a,b,c in orders:
        denominator.append([H,-U,-V][a]*[-(U+V-z),H-y,H-x][c])
    records=[]
    for i in range(6):
        for j in range(i,6):
            numerator=S.expand(sum(c*S.prod(dots[p] for p in mon)
                                   for mon,c in interference(i,j)))
            numerator=S.factor(numerator)
            records.append(dict(i=i,j=j,multiplicity=1 if i==j else 2,
                                numerator=str(numerator),
                                denominator=str(denominator[i]*denominator[j])))
            print(i,j,len(str(numerator)),flush=True)
    out=Path(__file__).with_name('compton_real_analytic_trace.json')
    write_artifact(out , json.dumps(dict(symbols='H=s-1,x=2P.k1,y=2P.k2,U=2p.k1,V=2p.k2',
             orders=orders,terms=records),indent=2)+'\n')
    print(trace.cache_info(),flush=True)

if __name__=='__main__': main()
