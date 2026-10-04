"""Analytic outer-t primitives for the UW endpoint coefficient and finite match."""
from compton_artifacts import read_artifact
import json
from pathlib import Path
from collections import defaultdict
import sympy as S
import mpmath as mp
from compton_soft_uv_check import rational_parts,primitive,combine,checked_physical_real
from compton_soft_analytic_gpl import evaluate_gpl_taylor


def endpoint_coefficients(H,component):
    d=json.loads(read_artifact(Path(__file__).with_name('compton_real_analytic_UW_endpoints.json')))
    hs,t=S.symbols('H t');h=S.Rational(str(H))
    LH1=S.sympify(d['LH1']).subs(hs,h);LH2=S.sympify(d['LH2']).subs(hs,h)
    TH1=S.sympify(d['TH1']).subs(hs,h);TH2=S.sympify(d['TH2']).subs(hs,h)
    Rlog=LH1;Rrat=-LH2*t/(1+h*t)-TH1
    if component=='finite_endpoint':
        C=-1/h+t*t/(1+h*t)
        Rlog=Rlog*S.log(h)+LH2*C
        Rrat=Rrat*S.log(h)-TH2*C
    return t,[S.cancel(Rrat),S.cancel(Rlog)]


def build(H,component):
    t,coefficients=endpoint_coefficients(H,component);h=mp.mpf(H);letter=-1/h
    out=defaultdict(mp.mpc)
    for R,word in zip(coefficients,[(),(letter,)]):
        for kind,data,co in rational_parts(R,t):
            if kind=='f':
                a,n=data
                if abs(a-letter)<mp.mpf(10)**(-mp.mp.dps+8):a=letter
                data=(a,n)
            combine(out,primitive(kind,data,word,mp.mpf(1)),co)
    c0=out.get((mp.mpf(0),),mp.mpf(0));scale=max(1,max(abs(c) for c in out.values()))
    if abs(c0)<scale*mp.mpf(10)**(-mp.mp.dps+15):out.pop((mp.mpf(0),),None)
    else:raise ValueError('uncancelled origin logarithm')
    out={word:co for word,co in out.items() if abs(co)>mp.mpf(10)**(-mp.mp.dps+15)}
    integrand=coefficients[0]+coefficients[1]*S.log(1+S.Rational(str(H))*t)
    return out,integrand,t


def uw_endpoint_closed(H,component='pole',dps=40,diagnostics=False):
    guards=25+int(mp.ceil(6*abs(mp.log10(mp.mpf(H)))))
    with mp.workdps(dps+guards):
        expr,_,_=build(H,component)
        value,_,_=evaluate_gpl_taylor(expr,mp.mpf(1))
        result,metadata=checked_physical_real(value,dps,'UW '+component)
        return (result,metadata) if diagnostics else result

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--H',default='9');p.add_argument('--digits',type=int,default=65);p.add_argument('--component',choices=['pole','finite_endpoint'],default='pole');p.add_argument('--check',action='store_true');a=p.parse_args()
    mp.mp.dps=a.digits;ex,f,t=build(a.H,a.component);value,nwords,nsteps=evaluate_gpl_taylor(ex,mp.mpf(1))
    print(a.H,a.component,'terms',len(ex),'words',nwords,'series steps',nsteps,flush=True)
    print('analytic',mp.nstr(value,a.digits-10),flush=True)
    if a.check:
        ser=S.series(f,t,0,14).removeO();low=S.series(f,t,0,2)
        print('origin',low,flush=True)
        fun=S.lambdify(t,f,'mpmath');local=S.lambdify(t,ser,'mpmath');cut=mp.mpf('1e-8')
        with mp.workdps(a.digits+50):
            check=mp.quad(lambda x:local(x) if x<cut else fun(x),[0,cut,mp.mpf('.5'),1])
        print('check',mp.nstr(check,a.digits-10),'difference',mp.nstr(value-check,8),flush=True)
