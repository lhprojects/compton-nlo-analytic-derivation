"""Analytic GPL integration of the UV difference's final v coordinate.

Consumes the independently generated inner-integral Laurent coefficients.
No source-paper coefficients are used. Quadrature is optional verification.
"""
from compton_artifacts import read_artifact
import json
from pathlib import Path
from collections import defaultdict
from functools import lru_cache
import sympy as S
import mpmath as mp
from compton_soft_analytic_gpl import evaluate_gpl_taylor,combine


@lru_cache(None)
def shuffle(a,b):
    if not a:return (b,)
    if not b:return (a,)
    return tuple((a[0],)+w for w in shuffle(a[1:],b))+tuple((b[0],)+w for w in shuffle(a,b[1:]))


def strip_trailing_zeros(word,logz):
    """Shuffle reduction into regular-origin words times powers of log(z)."""
    if not word:return {():mp.mpf(1)}
    if word[-1]!=0:return {word:mp.mpf(1)}
    if all(a==0 for a in word):return {():logz**len(word)/mp.factorial(len(word))}
    prefix=word[:-1];allwords=shuffle(prefix,(mp.mpf(0),));multiplicity=allwords.count(word)
    out=defaultdict(mp.mpc)
    combine(out,strip_trailing_zeros(prefix,logz),logz/multiplicity)
    for other in allwords:
        if other!=word:combine(out,strip_trailing_zeros(other,logz),-1/mp.mpf(multiplicity))
    return dict(out)


def word_series(word,degree):
    if not word:return [mp.mpf(1)]+[mp.mpf(0)]*degree
    if word[-1]==0:
        # Only the log-free Taylor coefficients enter the finite-part constant.
        out=[mp.mpf(0)]*(degree+1)
        for regular,weight in strip_trailing_zeros(word,mp.mpf(0)).items():
            if weight:
                values=word_series(regular,degree)
                for n in range(degree+1):out[n]+=weight*values[n]
        return out
    tail=word_series(word[1:],degree);a=word[0];out=[mp.mpf(0)]
    if a==0:
        assert len(word)>1
        out.extend(tail[n]/n for n in range(1,degree+1))
    else:
        for n in range(degree):out.append((tail[n]-n*out[n])/(-a*(n+1)))
    return out


def primitive(kind,data,word,Z):
    """Finite-part normalized rational/GPL primitive, including z=0 poles."""
    out=defaultdict(mp.mpc)
    if not word:
        if kind=='p':out[()]=Z**(data+1)/(data+1)
        else:
            a,n=data
            if n==1:out[(a,)]=1
            elif a==0:out[()]=Z**(1-n)/(1-n)
            else:out[()]=((Z-a)**(1-n)-(-a)**(1-n))/(1-n)
        return out
    b,tail=word[0],word[1:]
    if kind=='p':
        n=data;out[word]+=Z**(n+1)/(n+1)
        for j in range(n+1):combine(out,primitive('p',n-j,tail,Z),-b**j/(n+1))
        combine(out,primitive('f',(b,1),tail,Z),-b**(n+1)/(n+1))
    else:
        a,n=data
        if n==1:out[(a,)+word]=1;return out
        m=n-1;out[word]-=1/(m*(Z-a)**m)
        if a==0:out[()]+=word_series(word,m)[m]/m
        if a==b:combine(out,primitive('f',(a,m+1),tail,Z),1/mp.mpf(m))
        else:
            for j in range(1,m+1):
                combine(out,primitive('f',(a,j),tail,Z),(-1)**(m-j)/(a-b)**(m-j+1)/m)
            combine(out,primitive('f',(b,1),tail,Z),(-1)**m/(a-b)**m/m)
    return out


def rational_parts(expr,z):
    num,den=map(lambda p:S.Poly(p,z),S.fraction(S.cancel(expr)))
    P,_=S.div(num,den)
    def cv(value):
        re,im=S.N(value,mp.mp.dps+10).as_real_imag()
        return mp.mpc(str(re),str(im))
    terms=[('p',int(power[0]),cv(co)) for power,co in P.terms() if co]
    # These GPL kernels have linear and quadratic denominator factors.
    # Calling roots() on the full numeric rational polynomial first invokes
    # integer-basis heuristics, which can factor enormous integer coefficients.
    # Factor the polynomial, then use the exact elementary root formulas.
    roots=[]
    for factor,multiplicity in S.factor_list(den.as_expr(),z)[1]:
        polynomial=S.Poly(factor,z)
        if polynomial.degree()==1:
            a,b=polynomial.all_coeffs()
            factor_roots={-b/a:1}
        elif polynomial.degree()==2:
            a,b,c=polynomial.all_coeffs()
            discriminant=S.cancel(b*b-4*a*c)
            if discriminant==0:
                factor_roots={-b/(2*a):2}
            else:
                radical=S.sqrt(discriminant)
                factor_roots={(-b+radical)/(2*a):1,(-b-radical)/(2*a):1}
        else:
            factor_roots=S.roots(polynomial.as_expr(),z)
        roots.extend((cv(a),order*multiplicity) for a,order in factor_roots.items())
    assert sum(order for _,order in roots)==den.degree()
    nc=[cv(v) for v in num.all_coeffs()];lead=cv(den.LC())
    for a,order in roots:
        def reduced(v):
            result=mp.polyval(nc,v)/lead
            for b,n in roots:
                if b!=a:result/=(v-b)**n
            return result
        for n in range(1,order+1):
            co=mp.diff(reduced,a,order-n)/mp.factorial(order-n)
            terms.append(('f',(a,n),co))
    return terms


def logarithm_words(powers,H):
    H=mp.mpf(H);rp=(H+2+mp.sqrt(H*(H+4)))/2;rm=1/rp
    l1={(-mp.mpf(1),):mp.mpf(1)}
    l2={(rp-1,):mp.mpf(1),(rm-1,):mp.mpf(1),(-mp.mpf(1),):-mp.mpf(1)}
    result={():mp.mpf(1)}
    for times,terms in zip(powers,[l1,l2]):
        for _ in range(times):
            following=defaultdict(mp.mpc)
            for left,lc in result.items():
                for right,rc in terms.items():
                    for word in shuffle(left,right):following[word]+=lc*rc
            result=dict(following)
    return result


def inner_records(component,H):
    hs,v=S.symbols('H v');h=S.Rational(str(H));kk=(v-1)**2/(h*v)
    groups=defaultdict(lambda:S.Integer(0))
    for name,outer,sign in [('J01',1,1),('UV_taylor',0,-1)]:
        path=Path(__file__).with_name('compton_real_inner_'+name+'.json')
        for record in json.loads(read_artifact(path)):
            c=S.sympify(record['coefficient']).subs(hs,h)*sign;n=record['power']
            if record['kind']=='polynomial':
                if component=='finite':groups[(outer,0)]+=c*h**(n+1)/(n+1)
                continue
            root=S.sympify(record['root']).subs(hs,h)
            if root==0:
                if n==1:
                    if component=='pole':groups[(outer,0)]+=c
                    else:groups[(outer,0)]+=2*c*S.log(h)
                elif component=='finite':groups[(outer,0)]+=c*h**(1-n)*(kk**(n-1)-2)/(n-1)
            elif root==h:
                if n==1:
                    if component=='pole':groups[(outer,0)]-=c
                    else:
                        groups[(outer,1)]+=c
                        groups[(outer,0)]-=2*c*S.log(h)
                elif component=='finite':groups[(outer,0)]+=c*(-1)**n*h**(1-n)*(-2+(-kk/(1-kk))**(n-1))/(n-1)
            else:
                assert S.cancel(root-h/kk)==0
                if component=='finite':
                    if n==1:groups[(outer,1)]+=c
                    else:groups[(outer,0)]+=c*((h-root)**(1-n)-(-root)**(1-n))/(1-n)
    return [dict(powers=list(powers),coefficient=str(S.factor(S.cancel(c)))) for powers,c in groups.items() if c!=0]


def build(path,component,H,Z=None,verbose=True):
    records=inner_records(component,H) if path=='inner' else json.loads(read_artifact(Path(path)))[component]
    h=S.Rational(str(H));Hnum=mp.mpf(H);Z=Hnum if Z is None else mp.mpf(Z)
    hs,v,z=S.symbols('H v z');out=defaultdict(mp.mpc)
    integrand=0;checks=[]
    for record in records:
        coefficient=S.sympify(record['coefficient']).subs(hs,h).subs(v,z+1)
        parts=rational_parts(coefficient,z)
        words=logarithm_words(record['powers'],Hnum)
        canonical=[mp.mpf(0),-mp.mpf(1),-mp.mpf(2)]
        canonical+=list(logarithm_words([0,1],Hnum).keys())
        canonical=[a[0] if isinstance(a,tuple) else a for a in canonical]
        normalized=[]
        for kind,data,c in parts:
            if kind=='f':
                root,n=data
                for a in canonical:
                    if abs(root-a)<mp.mpf(10)**(-mp.mp.dps+8)*(1+abs(a)):
                        root=a;break
                data=(root,n)
            normalized.append((kind,data,c))
        parts=normalized
        for word,weight in words.items():
            for kind,data,c in parts:combine(out,primitive(kind,data,word,Z),weight*c)
        a,b=record['powers'];k=z*z/(h*(z+1))
        integrand+=coefficient*S.log(1+z)**a*S.log(1-k)**b
        target=S.lambdify(z,coefficient,'mpmath')
        zz=Z*mp.mpf('.371')
        rec=mp.mpc(0)
        for kind,data,c in parts:
            if kind=='p':rec+=c*zz**data
            else:
                root,n=data;rec+=c/(zz-root)**n
        checks.append(abs(rec-target(zz)))
    # An uncancelled G(0;z)=ln z would indicate a nonintegrable primitive.
    logarithm=out.get((mp.mpf(0),),mp.mpf(0))
    if verbose:print('partial-fraction reconstruction',mp.nstr(max(checks),7),'logarithm coefficient',mp.nstr(logarithm,7),flush=True)
    scale=max([abs(c) for c in out.values()]+[mp.mpf(1)])
    if abs(logarithm)<mp.mpf(10)**(-mp.mp.dps+15)*scale:out.pop((mp.mpf(0),),None)
    else:raise ValueError('nonzero logarithmic endpoint primitive')
    out={word:co for word,co in out.items() if abs(co)>mp.mpf(10)**(-mp.mp.dps+15)}
    return out,integrand,z


def checked_physical_real(value,dps,label):
    """Reject a wrong branch or lost cancellation before taking the real part."""
    residual=abs(mp.im(value))
    tolerance=mp.power(10,-dps)*max(1,abs(value))/100
    if residual>tolerance:
        raise ArithmeticError(label+' has unexpected imaginary residual '+str(residual))
    return +mp.re(value), {'imaginary_residual':+residual,'imaginary_tolerance':+tolerance}


def uv_difference_closed(H,component='finite',dps=40,diagnostics=False):
    """Analytic last-v integral of the complete UV difference sector.

Returns the coefficient before any overall sector prefactor supplied by the
parent real-emission calculation. No quadrature or source-paper result is used.
"""
    if component not in ('pole','finite'):raise ValueError('component must be pole or finite')
    if mp.mpf(H)<=0:raise ValueError('H must be positive')
    guards=25+int(mp.ceil(6*abs(mp.log10(mp.mpf(H)))))
    with mp.workdps(dps+guards):
        expr,_,_=build(str(Path(__file__).with_name('compton_real_uv_analytic.json')),component,H,verbose=False)
        value,words,steps=evaluate_gpl_taylor(expr,mp.mpf(H))
        result,imaginary=checked_physical_real(value,dps,'UV '+component)
        data={'terms':len(expr),'gpl_words':words,'continuation_steps':steps,'working_dps':mp.mp.dps,**imaginary}
    return (result,data) if diagnostics else result


def real_uv_closed(H,dps=40):
    """Return {'pole': P, 'finite': F} before the overall 1/(8*s*H).

Both entries are evaluated from finite GPL expressions. No quadrature is used.
"""
    return {component:uv_difference_closed(H,component,dps) for component in ('pole','finite')}


def main():
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--input',default=str(Path(__file__).with_name('compton_real_uv_analytic.json')));p.add_argument('--component',choices=['pole','finite'],default='pole');p.add_argument('--H',default='9');p.add_argument('--digits',type=int,default=70);p.add_argument('--check',action='store_true');args=p.parse_args()
    mp.mp.dps=args.digits
    expr,integrand,z=build(args.input,args.component,args.H)
    value,nwords,nsteps=evaluate_gpl_taylor(expr,mp.mpf(args.H))
    print('component',args.component,'H',args.H,'terms',len(expr),'words',nwords,'series steps',nsteps,flush=True)
    print('analytic GPL endpoint',mp.nstr(value,args.digits-10),flush=True)
    if args.check:
        # The input combination is proved regular by its local series first.
        ser=S.series(integrand,z,0,2).removeO().expand()
        print('combined origin series',ser,flush=True)
        assert all(power>=0 for power in [term.as_powers_dict().get(z,0) for term in S.Add.make_args(ser)])
        f=S.lambdify(z,integrand,'mpmath');h=mp.mpf(args.H)
        local=S.series(integrand,z,0,14).removeO()
        flocal=S.lambdify(z,local,'mpmath');cut=min(h/100,mp.mpf('1e-8'))
        def stable(q):return flocal(q) if abs(q)<cut else f(q)
        # The analytic Taylor polynomial removes cancellation at the origin.
        with mp.workdps(args.digits+50):
            check=mp.quad(stable,[0,cut,h/4,h/2,h])
        print('quadrature cross-check',mp.nstr(check,args.digits-10),'difference',mp.nstr(value-check,8),flush=True)

if __name__=='__main__':main()
