"""Generate the finite GPL endpoint expression, then independently check it.

The expression generator uses rational partial fractions and finite integration
by parts only. The ODE evaluates the resulting special functions for checking;
it does not integrate the Compton soft kernel used to derive the expression.
"""
import mpmath as mp
import sympy as sp
from collections import defaultdict
from compton_soft_analytic import soft_log_letters,soft_rational_kernels,born_angular,soft_angular_reduced


def kernels(s):
    w=sp.Symbol('w');x=sp.Rational(str(s))-1;ss=x+1
    beta=x/(x+2);U=(w-1/w)/(2*sp.sqrt(x));kk=(w+1/w)**2/(4*ss)
    D=1-kk
    a1=-(x**3-6*x*x-16*x-8)/x**2;a2=D*(x*x+2*x+4)/x;a3=-x*D**2
    b1=(x**3-2*x*x-16*x-8)/x**2;b2=4*D*(x+2)/x**2
    h0=a2*U/(2*(1-U*U))+a3*(U/(4*(1-U*U)**2)+3*U/(8*(1-U*U)))+b2*U/(2*kk*(kk-U*U))
    h1=a1+a2/2+3*a3/8
    hk=b1*2*sp.sqrt(ss)*w/(w*w+1)+b2*4*ss**sp.Rational(3,2)*w**3/(w*w+1)**3
    pref=4/beta*sp.diff(U,w)/((x+2)**2/(4*ss)*(1-beta**2*U**2))
    return w,[sp.cancel(pref*h) for h in [h0,h1,hk]]


def partials(s):
    w,rs=kernels(s)
    s=mp.mpf(s);x=s-1;beta=x/(s+1);a=mp.sqrt(x)/beta
    r=1/(mp.sqrt(1+a*a)+a)
    common=[(r,1),(-r,1),(1/r,1),(-1/r,1)]
    poles=[[(mp.mpf(0),3),(mp.j,1),(-mp.j,1)]+common,
           [(mp.mpf(0),4)]+common,
           [(mp.j,2),(-mp.j,2)]+common]
    out=[]
    def cv(value):return mp.mpf(str(sp.N(value,mp.mp.dps+5)))
    for R,roots in zip(rs,poles):
        num,den=map(lambda t:sp.Poly(t,w),sp.fraction(R))
        P,_=sp.div(num,den)
        # Store polynomial part in z=w-1, matching GPL origin.
        Z=sp.Symbol('z');pol=sp.Poly(P.as_expr().subs(w,Z+1),Z)
        terms=[('p',int(power[0]),cv(coeff)) for power,coeff in pol.terms() if coeff]
        nc=[cv(v) for v in num.all_coeffs()];lead=cv(den.LC())
        for root,order in roots:
            def reduced(v):
                value=mp.polyval(nc,v)/lead
                for other,degree in roots:
                    if other!=root:value/=(v-other)**degree
                return value
            for n in range(1,order+1):
                coefficient=mp.diff(reduced,root,order-n)/mp.factorial(order-n)
                terms.append(('f',(root-1,n),coefficient))
        out.append(terms)
    return out


def combine(result,terms,factor=1):
    for word,value in terms.items():result[word]+=factor*value


def primitive(kind,data,word,Z):
    """Algebraic GPL reduction, returned as {word: endpoint coefficient}."""
    result=defaultdict(mp.mpc)
    if not word:
        if kind=='p':result[()]=Z**(data+1)/(data+1)
        else:
            a,n=data
            if n==1:result[(a,)]=1
            else:result[()]=((Z-a)**(1-n)-(-a)**(1-n))/(1-n)
        return result
    b,tail=word[0],word[1:]
    if kind=='p':
        n=data;result[word]+=Z**(n+1)/(n+1)
        for j in range(n+1):combine(result,primitive('p',n-j,tail,Z),-b**j/(n+1))
        combine(result,primitive('f',(b,1),tail,Z),-b**(n+1)/(n+1))
    else:
        a,n=data
        if n==1:result[(a,)+word]=1;return result
        m=n-1;result[word]-=1/(m*(Z-a)**m)
        if a==b:combine(result,primitive('f',(a,m+1),tail,Z),1/mp.mpf(m))
        else:
            for k in range(1,m+1):
                coefficient=(-1)**(m-k)/(a-b)**(m-k+1)
                combine(result,primitive('f',(a,k),tail,Z),coefficient/m)
            combine(result,primitive('f',(b,1),tail,Z),(-1)**m/(a-b)**m/m)
    return result


def expression(s,Z):
    parts=partials(s);letters=soft_log_letters(s);out=defaultdict(mp.mpc)
    for a,ca in letters[0]:
        for kind,data,co in parts[0]:combine(out,primitive(kind,data,(a,),Z),ca*co)
        for j in (1,2):
            for b,cb in letters[j]:
                for kind,data,co in parts[j]:
                    combine(out,primitive(kind,data,(a,b),Z),ca*cb*co)
                    combine(out,primitive(kind,data,(b,a),Z),ca*cb*co)
    return dict(out)


def evaluate_gpl_ode(expr,Z):
    import numpy as np
    from scipy.integrate import solve_ivp
    words=set()
    for word in expr:
        while word:
            words.add(word);word=word[1:]
    words=sorted(words,key=lambda word:(len(word),str(word)));index={word:i for i,word in enumerate(words)}
    def rhs(z,values):
        result=np.empty(len(words),dtype=complex)
        for word,i in index.items():
            a=complex(word[0]);tail=word[1:]
            numerator=values[index[tail]] if tail else 1
            if z==0 and a==0:
                if len(tail)==1:result[i]=-1/complex(tail[0])
                elif len(tail)>1:result[i]=0
                else:raise ValueError('unexpected logarithm at origin')
            else:result[i]=numerator/(z-a)
        return result
    sol=solve_ivp(rhs,[0,float(Z)],np.zeros(len(words),dtype=complex),method='DOP853',rtol=2e-13,atol=2e-14)
    assert sol.success
    total=mp.mpc(0)
    for word,co in expr.items():total+=co*(complex(sol.y[index[word],-1]) if word else 1)
    return total,len(words),len(sol.t)


def evaluate_gpl_taylor(expr,Z):
    """Arbitrary-precision analytic continuation of convergent GPL series."""
    words=set()
    for word in expr:
        while word:
            words.add(word);word=word[1:]
    words=sorted(words,key=lambda word:(len(word),str(word)))
    values={word:mp.mpc(0) for word in words}
    letters=set(word[0] for word in words)
    order=int(mp.ceil((mp.mp.dps+12)/(-mp.log10(mp.mpf('.4')))))
    center=mp.mpf(0);steps=0
    while center<Z:
        distance=min(abs(center-a) for a in letters if center!=a)
        step=min(Z-center,mp.mpf('.4')*distance)
        series={():[mp.mpf(1)]+[mp.mpf(0)]*order}
        for word in words:
            a=word[0];tail=series[word[1:]]
            co=[values[word]]
            if center==a:
                assert a==0 and center==0 and tail[0]==0
                co.extend(tail[n]/n for n in range(1,order+1))
            else:
                den=center-a
                for n in range(order):co.append((tail[n]-n*co[n])/(den*(n+1)))
            series[word]=co
        for word in words:values[word]=mp.polyval(list(reversed(series[word])),step)
        center+=step;steps+=1
        if steps>1000:raise RuntimeError('series continuation failed')
    return sum(co*(values[word] if word else 1) for word,co in expr.items()),len(words),steps


def check(s,fraction,taylor=False):
    mp.mp.dps=60;s=mp.mpf(s);W=mp.sqrt(s)+mp.sqrt(s-1);Z=(W-1)*mp.mpf(fraction)
    p=partials(s)
    # Verify partial fractions independently from the explicit rational kernels.
    errors=[]
    for frac in [mp.mpf('.2'),mp.mpf('.7')]:
        z=Z*frac
        targets=soft_rational_kernels(s,z+1)
        for terms,target in zip(p,targets):
            value=mp.mpc(0)
            for kind,data,co in terms:
                if kind=='p':value+=co*z**data
                else:a,n=data;value+=co/(z-a)**n
            errors.append(abs(value-target))
    expr=expression(s,Z);value,nwords,nsteps=(evaluate_gpl_taylor(expr,Z) if taylor else evaluate_gpl_ode(expr,Z))
    letters=soft_log_letters(s)
    def original(z):
        ls=[sum(co*mp.log(1-z/a) for a,co in terms) for terms in letters]
        rs=soft_rational_kernels(s,z+1)
        return rs[0]*ls[0]+rs[1]*ls[0]*ls[1]+rs[2]*ls[0]*ls[2]
    direct=mp.quad(original,[0,Z/2,Z])
    print('s',s,'fraction',fraction,'partial-fraction error',mp.nstr(max(errors),7),flush=True)
    print('GPL expression terms',len(expr),'required words',nwords,'ODE steps',nsteps,flush=True)
    print('GPL',mp.nstr(value,23),'direct check',mp.nstr(direct,23),'difference',mp.nstr(value-direct,8),flush=True)
    x=s-1;beta=x/(s+1);w=Z+1;U=(w-1/w)/(2*mp.sqrt(x));delta=1-U
    b1=(x**3-2*x*x-16*x-8)/x**2
    h0=(x*x+2*x+4)/(2*s)-x**3/(4*s*s)+2*(x+2)/x
    tail=4*mp.atanh(beta)/beta*delta*(-2*mp.log(delta)+2+h0+2*mp.log(2)+b1*mp.log(s)/2)
    corrected=value+tail-mp.atanh(beta)/beta*born_angular(s)
    fullcheck=soft_angular_reduced(s)
    print('endpoint delta',mp.nstr(delta,8),'tail-corrected GPL',mp.nstr(corrected,23),'full check',mp.nstr(fullcheck,23),'difference',mp.nstr(corrected-fullcheck,8),flush=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--s',default='10');p.add_argument('--fraction',default='.9');p.add_argument('--taylor',action='store_true');args=p.parse_args()
    check(args.s,args.fraction,args.taylor)
