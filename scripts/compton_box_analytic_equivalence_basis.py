"""Exact reference-basis conversion for an analytic equality audit.

This file uses LSZ Eq. (6) only as the comparison target. It never changes
or supplies coefficients to the independently integrated cross section.
"""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
from collections import defaultdict
from functools import lru_cache
from types import SimpleNamespace
import ast,itertools,json
import sympy as S

x,q,Y=S.symbols('x q Y')
from compton_gpl_algebra import ROOT, shuffle, multiply, kernel_map, verify_kernels

def reference_x_words():
    # Read the already transcribed reference terms without executing quadrature.
    source=read_artifact(ROOT/'scripts/compton_reference_high_precision.py')
    tree=ast.parse(source)
    terms=next(n.value for n in ast.walk(tree) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='terms' for t in n.targets))
    l,h,d2,d21,p=S.symbols('l h d2 d21 p')
    words=[(-2,-1),(-1,-2,-1),(-1,-1,0),(-1,0,-1),(0,-1,-1),(0,'r',-1),('r',-1,0),('r',0,-1),('r',-1,-1),('r',-2,-1),('r',-1),(-4,'r',-1),('r','r',-1)]
    gs={w:S.Symbol('g'+str(k)) for k,w in enumerate(words)}
    env=dict(x=x,y=Y,l=l,h=h,li2=d2,li2_1=d21,poly=p,q=h/(1-x),g=gs,mp=SimpleNamespace(pi=S.pi))
    expr=sum(eval(compile(ast.Expression(terms),'<LSZ Eq6>','eval'),{'__builtins__':{}},env))/x**3
    basis={l:{(-1,):1},h:{(0,):1},d2:{(0,-1):-1},d21:{():S.pi**2/6,(1,0):-1},p:{(0,1,0):2,(0,-1,0):2}}
    basis.update({v:{w:1} for w,v in gs.items()})
    generators=list(basis);out=defaultdict(lambda:S.S.Zero)
    for powers,c in S.Poly(S.expand(expr),*generators).terms():
        piece={():c}
        for symbol,power in zip(generators,powers):
            for _ in range(power):piece=multiply(piece,basis[symbol])
        for w,v in piece.items():out[w]+=v
    return {w:S.factor(c) for w,c in out.items() if S.factor(c)!=0}


def convert():
    maps=kernel_map();verify_kernels(maps)
    original=reference_x_words();out=defaultdict(lambda:S.S.Zero)
    for w,c in original.items():
        coefficient=c.subs({x:(1-q)**2/q,Y:(1-q)/(1+q)})
        factors=[list(maps[a].items()) for a in w]
        for selection in itertools.product(*factors):
            letters=tuple(S.simplify(a) for a,_ in selection)
            out[letters]+=coefficient*S.prod(c for _,c in selection)
    out={w:S.factor(c) for w,c in out.items() if S.cancel(c)!=0}
    derivative=defaultdict(lambda:S.S.Zero)
    for w,c in out.items():
        derivative[w]+=S.diff(c,q)
        if w:derivative[w[1:]]+=c/(q-w[0])
    derivative={w:S.factor(c) for w,c in derivative.items() if S.cancel(c)!=0}
    serialize=lambda d:[dict(word=list(map(str,w)),coefficient=str(c)) for w,c in sorted(d.items(),key=lambda kv:(len(kv[0]),str(kv[0])))]
    report=dict(source='Lee, Schwartz and Zhang, arXiv:2102.06718, Eq. (6)',variable='x=(1-q)^2/q, 0<q<1',basepoint='q=1, real left tangential logarithm log(1-q)',normalization='sigma_1=alpha_T^3/m^2*C(s), x=s-1',kernel_maps={str(a):{str(b):str(c) for b,c in v.items()} for a,v in maps.items()},x_basis=serialize(original),q_basis=serialize(out),q_derivative=serialize(derivative),status='Reference conversion proved. Independent total has not yet been reduced to this basis.')
    target=ROOT/'compton_box_analytic_equivalence_basis.json'
    write_artifact(target , json.dumps(report,indent=2)+'\n')
    print('Reference x words:',len(original),'q words:',len(out),'derivative words:',len(derivative),flush=True)
    print('All six kernel pullbacks proved exactly.',flush=True)


if __name__=='__main__':convert()
