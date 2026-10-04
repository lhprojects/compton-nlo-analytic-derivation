from compton_artifacts import read_artifact
from compton_soft_uv_check import checked_physical_real
"""Outer UW finite integral as a generated finite GPL endpoint expression."""
import json
from pathlib import Path
from collections import defaultdict
import sympy as S
import mpmath as mp
from compton_soft_uv_check import rational_parts,primitive,combine,shuffle,strip_trailing_zeros
from compton_soft_analytic_gpl import evaluate_gpl_taylor

INPUT=Path(__file__).with_name('compton_real_analytic_UW_finite.json')


def multiply(left,right):
    out=defaultdict(mp.mpc)
    for a,ca in left.items():
        for b,cb in right.items():
            for word in shuffle(a,b):out[word]+=ca*cb
    return dict(out)


def atoms_gpl(H):
    h=mp.mpf(H);p=(1+mp.sqrt(1+4/h))/2;m=1-p;a=-1/h;b=-2/h
    d3=defaultdict(mp.mpc)
    for r in [p,m]:
        d3[(mp.mpf(1),r)]-=1;d3[(a,r)]-=1;d3[(mp.mpf(0),r)]+=1
    d3[(a,mp.mpf(0))]+=1
    return [{(a,):mp.mpf(1)}, {():mp.log(h)}, {():mp.log(2),(b,):mp.mpf(1)},
            {(p,):mp.mpf(1),(m,):mp.mpf(1)}, {(mp.mpf(0),):mp.mpf(1)},
            {(mp.mpf(0),a):-mp.mpf(1)}, {(a,):-mp.log(2),(a,b):-mp.mpf(1)},dict(d3)]


def atoms_direct(H,t):
    h=mp.mpf(H);t=mp.mpf(t)
    return [mp.log(1+h*t),mp.log(h),mp.log(2+h*t),mp.log(1+h*t-h*t*t),mp.log(t),
            mp.polylog(2,-h*t),mp.polylog(2,-1-h*t)+mp.pi**2/12,
            mp.polylog(2,-(1-t)*(1+h*t)/t)-mp.polylog(2,-(1-t)/t)]


def regularize(expr,Z):
    out=defaultdict(mp.mpc)
    for word,co in expr.items():combine(out,strip_trailing_zeros(word,mp.log(Z)),co)
    cutoff=mp.mpf(10)**(-mp.mp.dps+12)
    return {word:co for word,co in out.items() if abs(co)>cutoff}


def check_atoms(H,t):
    result=[]
    for expression,target in zip(atoms_gpl(H),atoms_direct(H,t)):
        ex=regularize(expression,mp.mpf(t))
        if set(ex).issubset({()}):value=ex.get((),0)
        else:value,_,_=evaluate_gpl_taylor(ex,mp.mpf(t))
        result.append(abs(value-target))
    return result


def coefficients(H):
    d=json.loads(read_artifact(INPUT));hs,t=S.symbols('H t');h=S.Rational(str(H))
    return t,[(record['powers'],S.cancel(S.sympify(record['coefficient']).subs(hs,h))) for record in d['finite']]


def expression(H,Z):
    h=mp.mpf(H);t,records=coefficients(H);atoms=atoms_gpl(H);out=defaultdict(mp.mpc)
    p=(1+mp.sqrt(1+4/h))/2
    canonical=[mp.mpf(0),mp.mpf(1),-1/h,-2/h,p,1-p,mp.j/mp.sqrt(h),-mp.j/mp.sqrt(h)]
    checks=[]
    for powers,R in records:
        words={():mp.mpf(1)}
        for power,atom in zip(powers,atoms):
            for _ in range(power):words=multiply(words,atom)
        parts=rational_parts(R,t)
        for kind,data,c in parts:
            if kind=='f':
                a,n=data
                for b in canonical:
                    if abs(a-b)<mp.mpf(10)**(-mp.mp.dps+8)*(1+abs(b)):
                        a=b;break
                data=(a,n)
            for word,weight in words.items():combine(out,primitive(kind,data,word,Z),c*weight)
    return regularize(out,Z)


def conv(a,b,n):
    return [sum(a[j]*b[k-j] for j in range(max(0,k-len(b)+1),min(k,len(a)-1)+1)) for k in range(n+1)]


def div(a,b,n):
    a=[mp.mpf(value) for value in a];b=[mp.mpf(value) for value in b]
    out=[]
    for k in range(n+1):out.append(((a[k] if k<len(a) else 0)-sum(b[j]*out[k-j] for j in range(1,min(k,len(b)-1)+1)))/b[0])
    return out


def log_series(a,n):
    derivative=[(j+1)*a[j+1] for j in range(len(a)-1)]
    quotient=div(derivative,a,max(0,n-1))
    return [mp.log(a[0])]+[quotient[j-1]/j for j in range(1,n+1)]


def li2_zero_series(a,n):
    result=[mp.mpf(0)]*(n+1);power=[mp.mpf(1)]+[mp.mpf(0)]*n
    for k in range(1,n+1):
        power=conv(power,a,n)
        for j in range(n+1):result[j]+=power[j]/(k*k)
    return result


def endpoint_atoms(H,n):
    h=mp.mpf(H)
    ell=log_series([1+h,-h],n);lh=[mp.log(h)]+[mp.mpf(0)]*n
    l2=log_series([2+h,-h],n);lp=log_series([1,h,-h],n);lt=log_series([1,-1],n)
    q0=div(ell,[1,-1],n-1);d0=[mp.polylog(2,-h)]+[q0[j-1]/j for j in range(1,n+1)]
    q2=div([h*a for a in l2],[1+h,-h],n-1);d2=[mp.polylog(2,-1-h)+mp.pi**2/12]+[q2[j-1]/j for j in range(1,n+1)]
    z1=div([0,-(1+h),h],[1,-1],n);z2=div([0,-1],[1,-1],n)
    d31=li2_zero_series(z1,n);d32=li2_zero_series(z2,n);d3=[a-b for a,b in zip(d31,d32)]
    return [ell,lh,l2,lp,lt,d0,d2,d3]


def endpoint_series(H,degree):
    t,records=coefficients(H);d=S.Symbol('d');trans=[];maximum=0
    def cv(x):return mp.mpf(str(S.N(x,mp.mp.dps+10)))
    for powers,R in records:
        num,den=map(lambda q:S.Poly(q,d),S.fraction(S.cancel(R.subs(t,1-d))))
        order=min(power[0] for power,c in den.terms() if c)
        trans.append((powers,num,den,order));maximum=max(maximum,order)
    N=degree+maximum;atoms=endpoint_atoms(H,N);result=defaultdict(mp.mpc)
    for powers,num,den,order in trans:
        nc=[cv(num.nth(j)) for j in range(N+1)]
        dc=[cv(den.nth(j+order)) for j in range(max(0,den.degree()-order)+1)]
        rational=div(nc,dc,N)
        prod=[mp.mpf(1)]+[mp.mpf(0)]*N
        for power,atom in zip(powers,atoms):
            for _ in range(power):prod=conv(prod,atom,N)
        full=conv(rational,prod,N)
        for j,value in enumerate(full):
            if j-order<=degree:result[j-order]+=value
    bad=max([abs(c) for j,c in result.items() if j<0]+[mp.mpf(0)])
    assert bad<mp.mpf(10)**(-mp.mp.dps+15),('uncancelled endpoint poles',bad)
    return [result[j] for j in range(degree+1)]


def origin_atoms(H,n):
    h=mp.mpf(H)
    ell=log_series([1,h],n);lh=[mp.log(h)]+[mp.mpf(0)]*n
    l2=log_series([2,h],n);lp=log_series([1,h,-h],n)
    d0=li2_zero_series([0,-h]+[mp.mpf(0)]*max(0,n-1),n)
    q2=div([-h*a for a in l2],[1,h],n-1)
    d2=[mp.mpf(0)]+[q2[j-1]/j for j in range(1,n+1)]
    A1=div(lp,[1,-1],n-1);A2=div([-h*a for a in lp],[1,h],n-1)
    analytic=[mp.mpf(0)];logarithm=[mp.mpf(0)]
    for j in range(1,n+1):
        aa=A1[j-1]+A2[j-1]+lp[j]
        bb=h*(-h)**(j-1)
        analytic.append(aa/j-bb/(j*j));logarithm.append(bb/j)
    return [{0:ell},{0:lh},{0:l2},{0:lp},{1:[mp.mpf(1)]+[mp.mpf(0)]*n},
            {0:d0},{0:d2},{0:analytic,1:logarithm}]


def origin_series(H,degree):
    t,records=coefficients(H);trans=[];maximum=0
    def cv(x):return mp.mpf(str(S.N(x,mp.mp.dps+10)))
    for powers,R in records:
        num,den=map(lambda q:S.Poly(q,t),S.fraction(R))
        order=min(power[0] for power,c in den.terms() if c)
        trans.append((powers,num,den,order));maximum=max(maximum,order)
    N=degree+maximum;atoms=origin_atoms(H,N);result=defaultdict(mp.mpc)
    for powers,num,den,order in trans:
        nc=[cv(num.nth(j)) for j in range(N+1)]
        dc=[cv(den.nth(j+order)) for j in range(max(0,den.degree()-order)+1)]
        rational=div(nc,dc,N);prod={0:[mp.mpf(1)]+[mp.mpf(0)]*N}
        for power,atom in zip(powers,atoms):
            for _ in range(power):
                following={}
                for l1,co1 in prod.items():
                    for l2,co2 in atom.items():
                        cc=conv(co1,co2,N)
                        if l1+l2 not in following:following[l1+l2]=cc
                        else:following[l1+l2]=[a+b for a,b in zip(following[l1+l2],cc)]
                prod=following
        for logs,co in prod.items():
            full=conv(rational,co,N)
            for j,value in enumerate(full):
                if j-order<=degree:result[(j-order,logs)]+=value
    bad=max([abs(c) for (j,l),c in result.items() if j<0]+[mp.mpf(0)])
    assert bad<mp.mpf(10)**(-mp.mp.dps+15),('uncancelled lower endpoint poles',bad)
    return {(j,l):c for (j,l),c in result.items() if j>=0}


def quadrature_check(H,dps=60,upper=1):
    # Verification only: the production closed routine does not call this.
    with mp.workdps(dps+45):
        t,records=coefficients(H);functions=[(powers,S.lambdify(t,S.factor(R),'mpmath')) for powers,R in records]
        lower=origin_series(H,16);upper_coefficients=endpoint_series(H,16);cut=mp.mpf('1e-6')/(1+mp.mpf(H))
        def f(q):
            if q<cut:return sum(c*q**n*mp.log(q)**l for (n,l),c in lower.items())
            if 1-q<cut:return sum(c*(1-q)**n for n,c in enumerate(upper_coefficients))
            values=atoms_direct(H,q);total=mp.mpf(0)
            for powers,R in functions:
                term=R(q)
                for power,value in zip(powers,values):term*=value**power
                total+=term
            return total
        upper=mp.mpf(upper)
        intervals=[0,cut,mp.mpf('.1'),mp.mpf('.5')]
        if upper==1:intervals += [1-cut,mp.mpf(1)]
        else:intervals += [upper]
        result=mp.quad(f,intervals)
        return +mp.re(result)


def real_uw_finite_closed(H,dps=40,diagnostics=False):
    h=mp.mpf(H)
    guards=30+int(mp.ceil(5*(mp.mpf(dps+18)/9+abs(mp.log10(1+h)))))
    with mp.workdps(dps+guards):
        degree=8;delta=mp.power(10,-mp.mpf(dps+18)/(degree+1))/(1+h)
        Z=1-delta;ex=expression(H,Z);partial,words,steps=evaluate_gpl_taylor(ex,Z)
        co=endpoint_series(H,degree+2)
        tail=sum(co[n]*delta**(n+1)/(n+1) for n in range(degree+1))
        error=abs(sum(co[n]*delta**(n+1)/(n+1) for n in range(degree+1,degree+3)))
        value,imaginary=checked_physical_real(partial+tail,dps,'UW finite')
        if error>mp.power(10,-dps)*max(1,abs(value))/100:
            raise ArithmeticError('insufficient analytic endpoint tail precision')
        return (value,{'tail_estimate':error,'words':words,'steps':steps,'delta':delta,**imaginary}) if diagnostics else value

def real_uw_closed(H,dps=40):
    """Return UW {'pole','finite'} before the overall 1/(8*s*H).

The crossed VW contribution is not included. Both entries use analytic GPL
expressions and convergent endpoint series, without numerical quadrature.
"""
    from compton_soft_uv_uw import uw_endpoint_closed
    return {'pole':uw_endpoint_closed(H,'pole',dps),'finite':real_uw_finite_closed(H,dps)}


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--H',default='9');p.add_argument('--digits',type=int,default=40);p.add_argument('--atoms',action='store_true');p.add_argument('--check',action='store_true');args=p.parse_args()
    mp.mp.dps=args.digits
    if args.atoms:
        for t in ['.1','.5','.9']:print(t,[mp.nstr(v,5) for v in check_atoms(args.H,t)],flush=True)
    else:
        value,data=real_uw_finite_closed(args.H,args.digits,True)
        print(mp.nstr(value,args.digits),data,flush=True)
        if args.check:
            check=quadrature_check(args.H,args.digits)
            print('quadrature check',mp.nstr(check,args.digits),'difference',mp.nstr(value-check,10),flush=True)
