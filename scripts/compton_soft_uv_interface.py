"""Production analytic soft angular coefficient: GPL plus convergent endpoint series.

No numerical quadrature is called. The GPLs themselves are evaluated by
convergent Taylor series with analytic continuation along the physical path.
"""
import mpmath as mp
from compton_soft_analytic import born_angular
from compton_soft_analytic_gpl import expression,evaluate_gpl_taylor


def endpoint_series(s,degree):
    s=mp.mpf(s);x=s-1;beta=x/(s+1);E=(s+1)/(2*mp.sqrt(s))
    a1=-(x**3-6*x*x-16*x-8)/x**2
    b1=(x**3-2*x*x-16*x-8)/x**2
    def analytic_parts(d):
        U=1-d;k=mp.sqrt((1+x*U*U)/s);D=x*(1-U*U)/s
        a2=D*(x*x+2*x+4)/x;a3=-x*D*D;b2=4*D*(x+2)/x**2
        h0=U*((x*x+2*x+4)/(2*s)-x**3/(4*s*s)-3*x**3*(1-U*U)/(8*s*s)+2*(x+2)/(x*k*k))
        h1=a1+a2/2+3*a3/8;hk=b1/k+b2/(2*k**3)
        pref=4/beta*mp.atanh(beta*U)/(E*E*(1-beta*beta*U*U))
        A=-pref*(h1+hk)/2
        B=pref*(h0+(h1-hk)*mp.log(2-d)/2+hk*mp.log(k+U)+hk*mp.log(s)/2)
        return A,B
    aa=mp.taylor(lambda d:analytic_parts(d)[0],0,degree)
    bb=mp.taylor(lambda d:analytic_parts(d)[1],0,degree)
    return aa,bb


def soft_angular_closed(s,dps=40,diagnostics=False):
    """Return C_K(s)=integral dc B K via its analytic GPL endpoint expression.

The endpoint tail is integrated termwise analytically to degree eight. Its
radius-scaled cutoff includes fifteen guard digits; two additional analytic
series terms provide an independent truncation estimate. `diagnostics=True`
returns (value, metadata). The physical domain is dimensionless s>1.
"""
    if dps<15:raise ValueError('use at least 15 decimal digits')
    requested_s=mp.mpf(s)
    if requested_s<=1:raise ValueError('physical s must exceed 1')
    logscale=max(0,int(mp.ceil(mp.log10(requested_s))))
    with mp.workdps(dps+30+4*logscale):
        ss=mp.mpf(s);x=ss-1;beta=x/(ss+1)
        degree=8
        radius=1/(ss+1)
        delta=radius*mp.power(10,-mp.mpf(dps+15)/(degree+1))
        U=1-delta;Z=mp.sqrt(x)*U+mp.sqrt(1+x*U*U)-1
        expr=expression(ss,Z);partial,words,steps=evaluate_gpl_taylor(expr,Z)
        aa,bb=endpoint_series(ss,degree+2)
        def tail(nmax):
            return sum(delta**(n+1)/(n+1)*(aa[n]*(mp.log(delta)-1/mp.mpf(n+1))+bb[n]) for n in range(nmax+1))
        t=tail(degree);correction=tail(degree+2)-t
        raw_value=partial+t-mp.atanh(beta)/beta*born_angular(ss)
        imaginary_residual=abs(mp.im(raw_value))
        imaginary_tolerance=mp.power(10,-dps)*max(1,abs(raw_value))/100
        if imaginary_residual>imaginary_tolerance:
            raise ArithmeticError("soft angular GPL has an unexpected imaginary part")
        value=mp.re(raw_value)
        error=abs(correction)
        target=mp.power(10,-dps)*max(1,abs(value))
        if error>target/100:raise ArithmeticError('endpoint series guard digits insufficient')
        metadata={'imaginary_residual':+imaginary_residual,'imaginary_tolerance':+imaginary_tolerance,'series_tail_estimate':+error,'endpoint_delta':+delta,'gpl_words':words,'continuation_steps':steps,'working_dps':mp.mp.dps}
        result=+value
    return (result,metadata) if diagnostics else result

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--s',default='10');p.add_argument('--digits',type=int,default=35);a=p.parse_args()
    value,info=soft_angular_closed(a.s,a.digits,True)
    print(mp.nstr(value,a.digits))
    for k,v in info.items():print(k,mp.nstr(v,8) if isinstance(v,mp.mpf) else v)
