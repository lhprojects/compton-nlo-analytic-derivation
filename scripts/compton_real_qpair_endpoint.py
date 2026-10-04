"""Finite endpoint of the summed auxiliary q-pair, from its own primitive identity."""
from compton_artifacts import write_artifact
import sympy as S
import mpmath as mp
from pathlib import Path
import json
r=S.symbols('r',positive=True);ss=r*r
Li=S.Function('L2');Log=S.Function('logabs')
def loga(z):
    z=S.factor(z)
    if z==1 or z==-1:return S.S.Zero
    n,d=S.fraction(z)
    def expand(p):
        c,f=S.factor_list(p,r)
        return S.log(abs(c))+sum(k*Log(S.factor(a)) for a,k in f)
    return S.expand(expand(n)-expand(d))
def l2(z):
    z=S.factor(z)
    if z==0:return S.S.Zero
    if z==1:return S.pi**2/6
    if z==-1:return -S.pi**2/12
    return Li(z)
def g(c,b):
    # Constant coefficient at u=s after expansion in delta=s-u.
    if c==ss:
        if b==ss:return loga(ss-r)**2/2
        return l2((ss-r)/(ss-b))-loga((ss-b)/(r-b))*loga(ss-r)
    if c==r:return -l2((ss-r)/(b-r))
    gc=loga((ss-c)/(r-c))
    if c==b:return gc**2/2
    return l2((r-c)/(b-c))-l2((ss-c)/(b-c))-loga((b-r)/(b-c))*gc
weights={S.S.Zero:-1,S.S.One:1,ss:-1,-ss:2,-S.S.One:-2}
outer={r:-1,-r:-1,S.S.Zero:1,ss:1,-S.S.One:2}
K=4*ss*(ss-1)/(ss+1)**2
expr=S.expand(S.pi**2/6+loga(K)**2/2+sum(c*k*g(a,b) for a,c in outer.items() for b,k in weights.items()))
terms={}
for atom in expr.atoms(Li):
    coef=S.expand(expr).coeff(atom)
    if coef:terms[str(atom.args[0])]=str(coef)
report={'variable':'r=sqrt(s)>1','derivation':'finite endpoint of eq:comp-real-qpair-reduction','Li2_terms':terms,
        'elementary':str(S.expand(expr.xreplace({a:0 for a in expr.atoms(Li)}))),
        'expression':str(expr),'branch':'L2(z)=Re Li2(z), logabs(z)=log|z|'}
folder=Path(__file__).resolve().parents[1]
write_artifact(folder/'real_qpair_endpoint.json' , json.dumps(report,indent=2)+'\n')
if __name__=='__main__':
    print(json.dumps(report,indent=2))
# Real reflection and inversion identities, applied only after proving the
# sign of each argument over the entire interval r>1.
h=S.symbols('h',positive=True)
def sign_domain(e):
    e=S.factor(e.subs(r,h+1));n,d=S.fraction(e)
    def sign_poly(p):
        cs=S.Poly(S.expand(p),h).all_coeffs()
        if all(v>=0 for v in cs) and any(v>0 for v in cs):return 1
        if all(v<=0 for v in cs) and any(v<0 for v in cs):return -1
        raise ValueError(('sign not established',p))
    return sign_poly(n)*sign_poly(d)
def unit(z):
    z=S.factor(z)
    assert sign_domain(z)>0 and sign_domain(1-z)>0
    if z==S.Rational(1,2):return S.pi**2/12-S.log(2)**2/2
    other=S.factor(1-z)
    # A fixed algebraic choice, not a sample-dependent branch choice.
    if (S.count_ops(other),str(other)) < (S.count_ops(z),str(z)):
        return S.pi**2/6-loga(z)*loga(1-z)-Li(other)
    return Li(z)
def canonical(z):
    z=S.factor(z)
    if z in [0,1,-1]:return l2(z)
    if sign_domain(z)<0:return -unit(z/(z-1))-loga(1-z)**2/2
    if sign_domain(1-z)<0:return S.pi**2/3-loga(z)**2/2-unit(1/z)
    return unit(z)
converted=S.expand(expr.xreplace({v:canonical(v.args[0]) for v in expr.atoms(Li)}))
converted=S.expand(converted.xreplace({v:Li(S.factor(v.args[0])) for v in converted.atoms(Li)}))
report['unit_interval_expression']=str(converted)
report['unit_interval_Li2_terms']={str(v.args[0]):str(converted.coeff(v)) for v in converted.atoms(Li) if converted.coeff(v)!=0}
report['unit_interval_elementary']=str(S.expand(converted.xreplace({v:0 for v in converted.atoms(Li)})))
report['argument_sign_certificate']='Signs proved from rational polynomials in h=r-1>0 with sign-definite coefficients'
write_artifact(folder/'real_qpair_endpoint.json' , json.dumps(report,indent=2)+'\n')
print('Reduced dilogarithm arguments:',len(report['unit_interval_Li2_terms']))
print(json.dumps(report['unit_interval_Li2_terms'],indent=2))
