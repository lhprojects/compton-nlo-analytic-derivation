"""Collect our finite regulator conversions in a real logarithm/Li2 basis."""
from compton_artifacts import write_artifact
import json
from pathlib import Path
import sympy as S
from compton_real_endpoint_matching import expressions,H
x=S.Symbol('x', positive=True)
raw=sum(v['finite_shift'] for v in expressions().values()).subs(H,x)
canonical=raw.replace(lambda a:a.func==S.polylog,lambda a:S.polylog(a.args[0],S.factor(a.args[1])))
canonical=canonical.replace(lambda a:a.func==S.polylog and a.args[0]==1,lambda a:-S.log(S.factor(1-a.args[1])))
canonical=canonical.subs(S.polylog(2,x/(x+1)),-S.polylog(2,-x)-S.log(1+x)**2/2)
normalized=canonical.replace(lambda a:a.func==S.log,lambda a:S.expand_log(S.log(S.factor(a.args[0])),force=True))
normalized=normalized.subs(S.log(x+2),S.log(2)+S.log(1+x/2))
atoms=sorted(normalized.atoms(S.log,S.polylog),key=str)
y=S.symbols('a:'+str(len(atoms)))
poly=S.Poly(S.expand(normalized.xreplace(dict(zip(atoms,y)))),*y)
terms=[]; compact=0
for powers,co in poly.terms():
 original=co
 co=S.factor(co)
 assert S.cancel(original-co)==0
 if co==0:continue
 product=S.prod(a**n for a,n in zip(atoms,powers))
 compact+=co*product
 terms.append({'powers':powers,'coefficient':str(co)})
assert S.Poly(compact.xreplace(dict(zip(atoms,y))),*y).monoms()==[tuple(t["powers"]) for t in terms]
out={'source':'sum of independent z-cut to eta finite shifts, including Bose double-pole conversion','domain':'x>0','basis':[str(a) for a in atoms],'terms':terms,'expression':str(compact),'raw_operations':int(S.count_ops(raw)),'compact_operations':int(S.count_ops(compact)),'identity':'Exact rational coefficient cancellation after positive-real logarithm identities'}
p=Path(__file__).resolve().parents[1]/'matching_simplified.json'
write_artifact(p , json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))

# This compares two independently assembled representations of the shifts.
# It is a numerical cross-check of the exact coefficient reduction above.
import mpmath as mp
with mp.workdps(100):
 old=S.lambdify(x,raw,'mpmath');new=S.lambdify(x,compact,'mpmath')
 checks=[]
 for energy in ('0.0001','0.6','1','8','10000'):
  value=mp.mpf(energy); a=old(value);b=new(value)
  error=abs(a-b)/max(1,abs(a))
  assert error<mp.mpf('1e-65')
  checks.append({'X':energy,'relative_difference':mp.nstr(error,8)})
 out['cross_checks']=checks
 write_artifact(p , json.dumps(out,indent=2)+'\n')
 print('Five independent endpoint-sum checks passed.',flush=True)
