"""Exact inverse word transformation of the independent virtual answer.

The tensor inverse of the six kernel pullbacks supplies every candidate word.
No reference coefficients, reference word list or reference file is read.
"""
from compton_artifacts import read_artifact, write_artifact
from collections import defaultdict
import itertools,json,time
import sympy as S
from compton_gpl_algebra import ROOT,kernel_map
q,X,Y,ll=S.symbols('q X Y ll')


def build():
    folder=ROOT
    vs=json.loads(read_artifact(folder/'compton_box_analytic_equivalence_virtual_total.json'))
    # JSON input letters and kernel letters must use the same algebraic form.
    # In particular 1/((3+sqrt(5))/2) simplifies to (3-sqrt(5))/2.
    maps={a:{S.simplify(b):v for b,v in terms.items()}
          for a,terms in kernel_map().items()}
    alphabet=list(maps)
    letters=sorted(set().union(*(set(c) for c in maps.values())),key=str)
    M=S.Matrix([[maps[a].get(b,0) for a in alphabet] for b in letters])
    pivots=M.T.rref()[1]
    assert len(pivots)==len(alphabet)
    inverse=M[list(pivots),:].inv()
    dual={b:{} for b in letters}
    for j,i in enumerate(pivots):
        dual[letters[i]]={a:inverse[k,j] for k,a in enumerate(alphabet) if inverse[k,j]}
    source={tuple(S.simplify(S.sympify(a)) for a in row['word']):S.sympify(row['coefficient']) for row in vs['independent_virtual']}
    result=defaultdict(lambda:S.S.Zero)
    for word,c in source.items():
        for choice in itertools.product(*(dual[b].items() for b in word)):
            result[tuple(a for a,_ in choice)]+=c*S.prod(v for _,v in choice)
    result={w:S.factor(S.cancel(c)) for w,c in result.items() if S.cancel(c)!=0}
    # Reconstruct ALL source coefficients, including letters outside the pivot
    # rows and words whose reconstructed coefficient vanishes.
    rebuilt=defaultdict(lambda:S.S.Zero)
    for word,c in result.items():
        for choice in itertools.product(*(maps[a].items() for a in word)):
            rebuilt[tuple(S.simplify(b) for b,_ in choice)]+=c*S.prod(v for _,v in choice)
    for word in set(source)|set(rebuilt):
        assert S.cancel(rebuilt[word]-source.get(word,0))==0,word
    original=sorted(result,key=lambda w:(len(w),str(w)))
    coefficients=[result[w] for w in original]
    print(f'{len(source)} source words reconstructed from {len(original)} independently selected X words.',flush=True)
    mod=q*q-(X+2)*q+1;out=[]
    for w,c in zip(original,coefficients):
        num,den=S.fraction(c)
        remainder=S.rem(S.rem(num,mod,q)*S.invert(den,mod,q),mod,q)
        aa=S.factor(remainder.coeff(q,0));bb=S.factor(remainder.coeff(q,1))
        coefficient=S.factor(aa+bb*(X+2)/2-bb*(X+4)*Y/2)
        out.append(dict(word=list(map(str,w)),coefficient=str(coefficient)))
    report=dict(source='Tensor inverse of six kernel pullbacks applied to independent virtual coefficients',candidate_source='All words through weight three in the six physical kernels; no literature word list',definition='X=s-1, Y=sqrt(X/(X+4)), ll=log(lambda/m)',coefficients=out,exact_116_equations_verified=True)
    write_artifact(folder/'compton_box_analytic_equivalence_virtual_x.json' , json.dumps(report,indent=2)+'\n')
    print('Converted all coefficient prefactors to Q(X)+sqrt(X/(X+4))*Q(X).',flush=True)


if __name__=='__main__':build()
