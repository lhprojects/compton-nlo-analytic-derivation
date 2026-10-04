"""Check the removable S=0 vertex endpoint against exact Taylor limits."""
from compton_artifacts import write_artifact
import json
from pathlib import Path
import sympy as S
import mpmath as mp
from compton_vertex_analytic import vertex_moments_mp

def check():
 a=S.Symbol('a');ell=S.log(1-a);d=S.polylog(2,a)
 P=[d,d+(1-1/a)*ell,2*a-2*d+(-a+1/a)*ell-1,
    d+(S.Rational(3,2)-2/a+1/(2*a*a))*ell+1/(2*a),
    5*a/2-3*d+(-a-S.Rational(3,2)+3/a-1/(2*a*a))*ell-1-1/(2*a),
    a*a-S.Rational(15,2)*a+6*d+(-a*a/2+4*a-4/a+1/(2*a*a))*ell+4+1/(2*a)]
 constants=[S.series(p,a,0,1).removeO() for p in P]
 assert constants==[0,1,-2,S.Rational(7,4),-S.Rational(15,4),S.Rational(31,4)]
 bounds=[S.pi**2/6,S.pi**2/6,1-S.pi**2/3,S.Rational(1,2)+S.pi**2/6,1-S.pi**2/2,-2+S.pi**2]
 powers=[1,1,2,1,2,3];expected=[(p-b)/(-1)**k for p,b,k in zip(constants,bounds,powers)]
 with mp.workdps(110):
  exact=vertex_moments_mp(0);errors=[abs(v-mp.mpf(str(e.evalf(105)))) for v,e in zip(exact,expected)]
  assert max(errors)<mp.mpf('1e-104')
  trials=[]
  for h in ['1e-12','1e-24','1e-36']:
   delta=mp.mpf(h)
   for sign in [-1,1]:
    values=vertex_moments_mp(sign*delta);err=max(abs(v-e) for v,e in zip(values,exact))
    # Linear approach follows from analyticity at S=0; both sides are real.
    assert err<10*delta
    trials.append(dict(S=mp.nstr(sign*delta,8),max_difference=mp.nstr(err,20)))
 report=dict(primitive_limits=list(map(str,constants)),moment_limits=list(map(str,expected)),J10_zero='pi**2/6 - 1',exact_limit_error=mp.nstr(max(errors),20),two_sided_checks=trials,impact='The total-cross-section production calls vertex_moments_mp only at S=s>1. Its crossed contribution uses the separate integrated primitive. Existing total curves are unchanged.')
 path=Path(__file__).resolve().parents[1]/'vertex_zero_check.json';write_artifact(path , json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2))
if __name__=='__main__':check()
