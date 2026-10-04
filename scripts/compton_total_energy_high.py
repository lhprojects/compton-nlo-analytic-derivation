"""High-energy logarithms of the independent total; expose unfixed constants."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S
from compton_box_analytic_equivalence_asymptotic import high_partial

folder=Path(__file__).resolve().parents[1]
X,Y,L,v=S.symbols('X Y L v')
rows=json.loads(read_artifact(folder/'total_energy_collected.json'))['terms']
values=high_partial([],return_values=True)
# Only these weight-two functions multiply coefficients of order X^0.
corrections={(-1,-1):L, (0,-1):-S.S.One, ('r',-1):2*L+1}
total=0
for row in rows:
    word=tuple(a if a=='r' else int(a) for a in row['word'])
    co=S.sympify(row['coefficient']).subs({X:1/v,Y:1/S.sqrt(1+4*v)})
    pref=S.series(co,v,0,2).removeO().expand()
    if pref.coeff(v,0)!=0:assert word in corrections or word==()
    total+=pref*(values[word]+v*corrections.get(word,0))
total=S.series(S.expand(total),v,0,2).removeO().expand()
assert S.simplify(total.coeff(v,0))==0
result=S.collect(total.coeff(v,1),L)/X
out={'source':'independent collected total and GPL kernel asymptotics',
     'leading_with_constants':str(result),
     'constants':'K_a...K_f are the finite large-X constants of the six root weight-three functions in the original independent asymptotic basis',
     'status':'Intermediate expression before endpoint constants are inserted; the exact evaluation is in total_energy_high.json.',
     'remainder':'O(log(X)^3/X^2) after the leading constants are evaluated'}
write_artifact(folder/'total_energy_high_partial.json' , json.dumps(out,indent=2)+'\n')
print(result,flush=True)
