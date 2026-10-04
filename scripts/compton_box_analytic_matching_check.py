"""Analytic moving-boundary soft matching, independently checked by exact toys.

No published cross section and no numerical quadrature enter the result.
"""
from compton_artifacts import read_artifact
import json
from pathlib import Path
import sympy as S
import mpmath as mp
H, r = S.symbols('H r')


def moments(f, upper):
    """Exact rational integral and logarithmic moment using pole derivatives."""
    plain = logmoment = S.S.Zero
    u = S.Dummy('pole')
    for term in S.Add.make_args(S.apart(f, r)):
        num, den = S.fraction(term)
        if not den.has(r):
            for (n,), c in S.Poly(term, r).terms():
                plain += c*upper**(n+1)/(n+1)
                logmoment += c*upper**(n+1)*(S.log(upper)/(n+1)-S.Rational(1,(n+1)**2))
            continue
        factors = [(p,n) for p,n in S.factor_list(den,r)[1] if p.has(r)]
        assert len(factors)==1
        p,n=factors[0]
        assert S.degree(p,r)==1
        pole=S.solve(p,r)[0]
        c=S.cancel(term*(r-pole)**n)
        assert not c.has(r)
        primitive=S.log(1-upper/u)
        logarithmic=S.log(upper)*primitive+S.polylog(2,upper/u)
        plain += c*S.diff(primitive,u,n-1).subs(u,pole)/S.factorial(n-1)
        logmoment += c*S.diff(logarithmic,u,n-1).subs(u,pole)/S.factorial(n-1)
    return plain,logmoment


def bose_matching(h=H):
    """Return closed P, A, log moment B, regulator correction C, finite shift.

    All coefficients are unnormalized integral coefficients.
    """
    h=S.sympify(h)
    if isinstance(h,S.Float):h=S.Rational(str(h))
    data=json.loads(read_artifact(Path(__file__).with_name('compton_real_bose_subtraction.json')))
    a=S.sympify(data['Y_double']).subs(H,h)
    b=S.sympify(data['Y_simple']).subs(H,h)
    q=h/(h+1);q1=-1/(h+1)**2
    A0,M0=moments(b,q)
    P=moments(r*a,q)[0]
    C=moments((1-r)*a,q)[0]/h
    A=A0+q1*a.subs(r,q)
    B=M0+q1*a.subs(r,q)*S.log(q)
    return dict(P=P,A=A,B=B,C=C,shift=A*S.log(h)-B-C)


def toy_checks(h=9,dps=90):
    mp.mp.dps=dps
    h=mp.mpf(h);s=h+1;q=h/s
    eps=mp.mpf('1e-22');cut=mp.mpf('1e-30')
    xp=(h+cut+mp.sqrt((h-cut)**2-4*cut))/2
    xm=s*cut/xp
    def mellin(t,n):
        return h**(2*t+2-n)*mp.gamma(t-n+1)*mp.gamma(t+1)/mp.gamma(2*t-n+3)*(mp.hyp2f1(t+1,t-n+1,2*t-n+3,q)-1)
    errors={}
    for n in (2,3):
        if n==2:
            A=q;P=0
            F=2*q*mp.log(h)-mp.log(s)/s
            FP=2*q*mp.log(h)-mp.log(s)-q
            zvalue=q*mp.log(xp/xm)+mp.log((s-xp)/(s-xm))/s+cut*(1/xp-1/xm)
            shift=q*mp.log(h)-q*(mp.log(q)-1)
        else:
            A=-1/s**2;P=q*q/2
            C=(q-q*q/2)/h
            F=-1/s-(2*mp.log(h)+mp.log(s))/s**2-C
            FP=-1/s-2*mp.log(h)/s**2
            def primitive(x):return -q/x-mp.log(x)/s**2+mp.log(s-x)/s**2+cut/(2*x*x)
            zvalue=primitive(xp)-primitive(xm)
            shift=-mp.log(s)/s**2-C
        errors[f'n{n}_mellin']=(mellin(eps,n)+mellin(-eps,n))/2-F
        errors[f'n{n}_zcut']=zvalue-P/cut+A*mp.log(cut)-FP
        errors[f'n{n}_matching']=FP+shift-F
    return errors


if __name__=='__main__':
    for name,value in bose_matching(S.Integer(9)).items():
        print(name,S.N(value,65),flush=True)
    for name,value in toy_checks().items():
        print(name,mp.nstr(value,8),flush=True)
