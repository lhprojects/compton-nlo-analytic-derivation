"""Laurent coefficients of the UV sector's analytically regulated x integral.

The seven rational coefficients come from the independent real trace.
The exact finite-regulator integral is Euler beta times Gauss 2F1.
No numerical quadrature is used in this reduction or its default check.
"""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
import json
import sympy as S
import mpmath as mp
H,v=S.symbols('H v');k=(v-1)**2/(H*v)


def coefficients(name='J01'):
    records=json.loads(read_artifact(Path(__file__).with_name('compton_real_inner_'+name+'.json')))
    pole=0;finite=0
    for r in records:
        c=S.sympify(r['coefficient']);n=r['power']
        if r['kind']=='polynomial':
            finite+=c*H**(n+1)/(n+1)
            continue
        a=S.sympify(r['root'])
        if a==0:
            if n==1:pole+=c;finite+=2*c*S.log(H)
            else:finite+=c*H**(1-n)*(k**(n-1)-2)/(n-1)
        elif a==H:
            if n==1:pole-=c;finite+=c*(S.log(1-k)-2*S.log(H))
            else:finite+=c*(-1)**n*H**(1-n)*(-2+(-k/(1-k))**(n-1))/(n-1)
        else:
            assert S.cancel(a-H/k)==0
            if n==1:finite+=c*S.log(1-k)
            else:finite+=c*((H-a)**(1-n)-(-a)**(1-n))/(1-n)
    return S.factor(pole),S.expand(finite)


def exact_monomial(kind,n,h,kk,eta):
    if kind=='constant':return h**(1+2*eta)*mp.beta(1+eta,1+eta)*mp.hyp2f1(eta,1+eta,2+2*eta,kk)
    if kind=='zero':return h**(1-n+2*eta)*mp.beta(1+eta-n,1+eta)*mp.hyp2f1(eta,1+eta-n,2+2*eta-n,kk)
    if kind=='H':return (-1)**n*h**(1-n+2*eta)*mp.beta(1+eta,1+eta-n)*mp.hyp2f1(eta,1+eta,2+2*eta-n,kk)
    return (-h/kk)**(-n)*h**(1+2*eta)*mp.beta(1+eta,1+eta)*mp.hyp2f1(n+eta,1+eta,2+2*eta,kk)


def check():
    mp.mp.dps=80;eta=mp.mpf('1e-25')
    for h,kk in [(mp.mpf('.6'),mp.mpf('.1')),(mp.mpf(9),mp.mpf('.6'))]:
        cases=[('constant',0,0,h),('zero',1,1,2*mp.log(h)),
               ('zero',2,0,(kk-2)/h),('H',1,-1,mp.log(1-kk)-2*mp.log(h)),
               ('H',2,0,-(2-kk)/(h*(1-kk))),('other',1,0,mp.log(1-kk)),
               ('other',2,0,kk*kk/(h*(1-kk)))]
        for kind,n,pole,finite in cases:
            residual=exact_monomial(kind,n,h,kk,eta)-pole/eta-finite
            print(h,kk,kind,n,mp.nstr(residual,6))
            assert abs(residual)<mp.mpf('1e-20')
    p,f=coefficients()
    write_artifact(Path(__file__).with_suffix('.json') , json.dumps({'pole':str(p),'finite':str(f)},indent=2)+'\n')

if __name__=='__main__':check()
