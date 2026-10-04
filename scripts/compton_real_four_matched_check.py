"""Check the shortened matched sum against our original endpoint and shift formulas."""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
import json
import sympy as S
import mpmath as mp
folder=Path(__file__).resolve().parents[1]
s,x,Ls,Lh,Lp=S.symbols('s x Ls Lh Lp');mp.mp.dps=75
modules=[{'L2':lambda z:mp.polylog(2,z),'Li3':lambda z:mp.polylog(3,z)},'mpmath']
d=json.loads(read_artifact(folder/'real_four_matched_short.json'))
a=[S.lambdify((s,Ls,Lh,Lp),S.sympify(v),modules) for v in d['atoms']]
c=[S.lambdify(s,S.sympify(t['coefficient']),'mpmath') for t in d['terms']]
m=json.loads(read_artifact(folder/'matching_simplified.json'))
mc=[S.lambdify(x,S.sympify(t['coefficient']),'mpmath') for t in m['terms']]
original=json.loads(read_artifact(folder/'real_four_finite_check.json'));rows=[]
for old in original['rows']:
    ss=mp.mpf(old['s']);h=ss-1
    av=[f(ss,mp.log(ss),mp.log(h),mp.log(ss+1)) for f in a]
    short=mp.fsum(f(ss)*mp.fprod(v**p for v,p in zip(av,t['powers']) if p) for f,t in zip(c,d['terms']))
    mv=[mp.log(ss),mp.log(h),mp.log((ss+1)/2),mp.polylog(2,-h),mp.polylog(2,-h*h/(4*ss))]
    shift=mp.fsum(f(h)*mp.fprod(v**p for v,p in zip(mv,t['powers']) if p) for f,t in zip(mc,m['terms']))
    # The stored independent original value has 45 significant digits.
    old_text=old['original_four_sum'].strip('()').split(' + ')[0].split(' - ')[0]
    reference=mp.mpf(old_text)+shift;error=abs(short-reference)
    assert error<mp.mpf('1e-40')*max(1,abs(reference)),(ss,error)
    rows.append({'s':str(ss),'matched':mp.nstr(short,60),'difference_from_original_plus_shifts':mp.nstr(error,9),
                 'reference_precision_note':'original GPL sum printed to 45 significant digits'})
write_artifact(folder/'real_four_matched_check.json' , json.dumps({'rows':rows,'comparison':'own original GPL endpoints plus own finite matching shifts'},indent=2)+'\n')
print(json.dumps(rows,indent=2))
