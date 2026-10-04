"""Analytic rational inner primitives for the real-radiation UV/UW sectors.

Inputs are independently generated angular coefficients. No literature total
or numerical integral is used. Logarithms are real log-absolute-values on the
physical interval; separate sectors must share the same infrared boundary.
"""
from compton_artifacts import artifact_path, read_artifact, write_artifact
from pathlib import Path
import json
import sympy as S
H,v,x=S.symbols('H v x')


def compute(name):
    source=Path(__file__).with_name('compton_real_analytic_'+name+'.txt')
    expr=S.sympify(read_artifact(source))
    cached=source.with_name(source.stem+'_apart.txt')
    parts=S.sympify(read_artifact(cached)) if artifact_path(cached).exists() else S.apart(expr,x)
    records=[]
    for term in S.Add.make_args(parts):
        num,den=S.fraction(term)
        poly=S.Poly(den,x)
        if poly.degree()==0:
            for (power,),c in S.Poly(term,x).terms():
                records.append(dict(kind='polynomial',coefficient=str(S.factor(c)),power=int(power)))
        else:
            # Each apart denominator is a power of a single linear factor.
            # Read its repeated root from its leading two coefficients.
            power=poly.degree()
            root=S.cancel(-poly.nth(power-1)/(power*poly.LC()))
            assert S.cancel(poly.as_expr()-poly.LC()*(x-root)**power)==0
            c=S.cancel(term*(x-root)**power)
            records.append(dict(kind='pole',coefficient=str(S.factor(c)),root=str(S.factor(root)),power=int(power)))
    reconstructed=sum(S.sympify(r['coefficient'])*(x**r['power'] if r['kind']=='polynomial' else (x-S.sympify(r['root']))**(-r['power'])) for r in records)
    assert S.cancel(reconstructed-expr)==0
    write_artifact(Path(__file__).with_name('compton_real_inner_'+name+'.json') , json.dumps(records,indent=2)+'\n')
    print(name,len(records),'exact rational terms',flush=True)
    # Each pole term integrates to c log|x-a| (power 1), or
    # c (x-a)^(1-n)/(1-n). Polynomial terms use c x^(n+1)/(n+1).
    return records

if __name__=='__main__':
    compute('J01')
    compute('J02')
