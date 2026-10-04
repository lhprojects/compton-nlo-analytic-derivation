"""Low-energy expansion from the complete compact independent formula."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
from collections import defaultdict
import sympy as S
from compton_box_analytic_equivalence_asymptotic import gseries

folder=Path(__file__).resolve().parents[1]
X,Y,L,M,t,ell=S.symbols('X Y L M t ell')
rows=json.loads(read_artifact(folder/'total_energy_compact.json'))['terms']
total=defaultdict(lambda:S.S.Zero)
for row in rows:
    word=tuple(a if a=='r' else int(a) for a in row['word'])
    co=S.sympify(row['coefficient']).subs({X:t*t,Y:t/S.sqrt(t*t+4),L:S.log(1+t*t),M:ell})
    pref=S.series(co,t,0,10).removeO().expand()
    for term in S.Add.make_args(pref):
        n=int(term.as_powers_dict().get(t,0));a=term/t**n
        for (p,k),b in gseries(word).items():
            if n+p<=8:total[n+p]+=a*b*ell**k
    print('expanded',word,flush=True)
total={n:S.factor(a) for n,a in total.items() if S.factor(a)!=0}
assert all(n%2==0 and n>=4 for n in total),total
assert all(S.degree(a,ell)<=1 for a in total.values()),total
expr=sum(a*X**(S.Rational(n,2)) for n,a in total.items())
report={'source':'complete compact independent cross section, expanded by recursive kernel integration',
        'log_symbol':'ell=ln(X)', 'C1_low':str(expr),
        'coefficients':{str(n//2):str(a) for n,a in sorted(total.items())},
        'remainder':'O(X^5*abs(log(X))) as X->0+',
        'certificate':'all negative powers, constant, linear, half-integer powers and log powers above one cancel exactly'}
write_artifact(folder/'total_energy_low.json' , json.dumps(report,indent=2)+'\n')
print('C1 =',expr,flush=True)
