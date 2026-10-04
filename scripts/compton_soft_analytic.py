"""Independent dilogarithmic soft function; quadrature is verification only."""
import mpmath as mp

def li2real(x):
    return mp.re(mp.polylog(2,x))

def log_pair(Y,g,h):
    L=mp.log(abs(1-Y/g))
    if g==h:
        return L*L/2
    return mp.log(abs(1-g/h))*L-li2real((Y-g)/(h-g))+li2real(-g/(h-g))

def soft_closed(s,t):
    s,t=mp.mpf(s),mp.mpf(t)
    E=(s+1)/(2*mp.sqrt(s)); q=(s-1)/(2*mp.sqrt(s)); beta=q/E
    a=mp.sqrt(-t)/2
    if a==0:
        return mp.mpf(0)
    if t == -4*q*q:
        w=mp.atanh(beta)
        J=(li2real(beta)-li2real(-beta))/2+w*mp.log(2*mp.cosh(w))-w*w/2-li2real(-mp.exp(-2*w))/2-mp.pi**2/24
        return (1+2*q*q)*J/(beta*E*E)-w/beta
    k=mp.sqrt(q*q-a*a)/E
    d=mp.sqrt(1-k*k)
    roots=[(-1-d)/k, -k/(1+d), k/(1+d), (1+d)/k]
    signs=[1,1,-1,-1]
    Y=(q+a)/mp.sqrt(q*q-a*a)-1
    L0=mp.log((1+k)/(1-k))
    result=mp.mpf(0)
    for z in roots:
        g=z-1
        residue=1/(4-2*k*k*(z*z+1))
        result+=residue*(L0*mp.log(abs(1-Y/g))+sum(sign*log_pair(Y,g,r-1) for r,sign in zip(roots,signs)))
    return (1-t/2)*result/(a*E)-mp.atanh(beta)/beta

def soft_quad(s,t):
    E=(s+1)/(2*mp.sqrt(s)); beta=(s-1)/(s+1)
    def integrand(x):
        M=1-t*x*(1-x); v=mp.sqrt(1-M/E**2)
        return (mp.atanh(v)/v if v else 1)/M
    return (1-t/2)*mp.quad(integrand,[0,mp.mpf('.5'),1])-mp.atanh(beta)/beta


def born(s,t):
    x=s-1; y=x+t
    return 2*(x/y+y/x-4*(1/y-1/x)+4*(1/y-1/x)**2)

def ir_angular_closed(s):
    s=mp.mpf(s); x=s-1; A=x-4-8/x
    B0=2*(s+5+4/x)/x
    roots=[(s+1+mp.sqrt((s-1)*(s+3)))/2,2/(s+1+mp.sqrt((s-1)*(s+3)))]
    def power(r,n,log):
        if n==-1:
            return mp.log(r)**2/2 if log else mp.log(r)
        return r**(n+1)*(mp.log(r)/(n+1)-1/mp.mpf(n+1)**2) if log else r**(n+1)/(n+1)
    def primitive(r,sign,log):
        coeff={1:-2/x,0:B0,-3:-sign*2/x,-2:sign*B0,-1: -2/x+sign*(-2/x-2*A)}
        total=sum(co*power(r,n,log) for n,co in coeff.items())
        for z in roots:
            yp=-1+1/z**2; ypp=-2/z**3
            f=1+sign/z**2; fp=-2*sign/z**3
            d=8*f/yp**2
            c=f*(2*A/yp-8*ypp/yp**3)+8*fp/yp**2
            if log:
                p1=mp.log(r)*mp.log(abs(1-r/z))+li2real(r/z)
                p2=-mp.log(r)/(r-z)+(mp.log(abs(r-z))-mp.log(r))/z
            else:
                p1=mp.log(abs(r-z)); p2=-1/(r-z)
            total+=c*p1+d*p2
        return total
    return 2*s/x**2*((primitive(s,1,True)-primitive(1,1,True))-(primitive(s,-1,False)-primitive(1,-1,False)))

def ir_angular_quad(s):
    s=mp.mpf(s);q=(s-1)/(2*mp.sqrt(s))
    def f(c):
        t=-2*q*q*(1-c); z=1-t/2
        b=mp.acosh(z)*z/mp.sqrt(z*z-1) if z!=1 else 1
        return born(s,t)*(b-1)
    return mp.quad(f,[-1,0,1])

if __name__=='__main__':
    mp.mp.dps=70
    for s in map(mp.mpf,['1.6','10','100','19570.511835591835']):
        q=(s-1)/(2*mp.sqrt(s))
        for fraction in map(mp.mpf,['0.0001','.2','.8','.9999']):
            t=-4*q*q*fraction
            answer=soft_closed(s,t); check=soft_quad(s,t)
            print('s=',s,'fraction=',fraction,'K=',mp.nstr(answer,25),'difference=',mp.nstr(answer-check,5),flush=True)

    for s in map(mp.mpf,['1.6','10','100','19570.511835591835']):
        value=ir_angular_closed(s); check=ir_angular_quad(s)
        print('IR angular s=',s,'value=',mp.nstr(value,25),'difference=',mp.nstr(value-check,5),flush=True)
        q=(s-1)/(2*mp.sqrt(s)); t=-4*q*q
        print('backward s=',s,'difference=',mp.nstr(soft_closed(s,t)-soft_quad(s,t),5),flush=True)

def soft_angular_reduced(s):
    """One-dimensional check of analytically exchanged angular integral."""
    s=mp.mpf(s);x=s-1;E=(s+1)/(2*mp.sqrt(s));beta=x/(s+1)
    def integrand(U):
        k=mp.sqrt((1+x*U*U)/s);D=1-k*k
        a1=-(x**3-6*x*x-16*x-8)/x**2
        a2=D*(x*x+2*x+4)/x;a3=-x*D*D
        b1=(x**3-2*x*x-16*x-8)/x**2;b2=4*D*(x+2)/x**2
        f11=mp.atanh(U);f21=U/(2*(1-U*U))+f11/2
        f31=U/(4*(1-U*U)**2)+3*f21/4
        f1k=mp.atanh(U/k)/k;f2k=U/(2*k*k*(k*k-U*U))+f1k/(2*k*k)
        H=a1*f11+a2*f21+a3*f31+b1*f1k+b2*f2k
        M=E*E*(1-beta*beta*U*U)
        return mp.atanh(beta*U)*H/M
    return 4/beta*mp.quad(integrand,[0,.5,1])-mp.atanh(beta)/beta*born_angular(s)

def born_angular(s):
    s=mp.mpf(s);x=s-1;A=x-4-8/x
    def p(y):return y*y/(2*x)+(4/x+4/x**2)*y+A*mp.log(y)-4/y
    return 4*s/x**2*(p(x)-p(x/s))

def soft_angular_check(s):
    s=mp.mpf(s);q=(s-1)/(2*mp.sqrt(s))
    def f(c):
        t=-2*q*q*(1-c)
        return born(s,t)*soft_closed(s,t)
    return mp.quad(f,[-1,0,1])

def soft_rational_kernels(s,w):
    """Explicit rational kernels in the final weight-three endpoint formula."""
    s=mp.mpf(s);x=s-1;beta=x/(s+1);E=(s+1)/(2*mp.sqrt(s))
    U=(w-1/w)/(2*mp.sqrt(x));k=(w+1/w)/(2*mp.sqrt(s))
    D=1-k*k
    a1=-(x**3-6*x*x-16*x-8)/x**2
    a2=D*(x*x+2*x+4)/x;a3=-x*D*D
    b1=(x**3-2*x*x-16*x-8)/x**2;b2=4*D*(x+2)/x**2
    h0=a2*U/(2*(1-U*U))+a3*(U/(4*(1-U*U)**2)+3*U/(8*(1-U*U)))+b2*U/(2*k*k*(k*k-U*U))
    h1=a1+a2/2+3*a3/8;hk=b1/k+b2/(2*k**3)
    pref=4/beta*(1+1/w**2)/(2*mp.sqrt(x))/(E**2*(1-beta**2*U*U))
    return tuple(pref*h for h in [h0,h1,hk])

def soft_log_letters(s):
    """GPL letters and weights, with origin shifted from w=1 to z=0."""
    s=mp.mpf(s);x=s-1;beta=x/(s+1);W=mp.sqrt(s)+mp.sqrt(x)
    A=mp.sqrt(x)/beta;r=1/(mp.sqrt(1+A*A)+A)
    def pairs(roots,weights):return [(root-1,mp.mpf(weight)/2) for root,weight in zip(roots,weights)]
    return (pairs([r,-1/r,1/r,-r],[1,1,-1,-1]),
            pairs([1/W,-W,W,-1/W],[1,1,-1,-1]),
            pairs([1/W,-1/W,W,-W],[1,1,-1,-1]))

def check_soft_letters(s):
    """Pure evaluations of the factored logarithms, without integration."""
    s=mp.mpf(s);x=s-1;beta=x/(s+1);W=mp.sqrt(s)+mp.sqrt(x)
    letters=soft_log_letters(s)
    errors=[]
    for fraction in map(mp.mpf,['.01','.4','.9','.999']):
        w=1+(W-1)*fraction;z=w-1
        U=(w-1/w)/(2*mp.sqrt(x));k=(w+1/w)/(2*mp.sqrt(s))
        target=[mp.atanh(beta*U),mp.atanh(U),mp.atanh(U/k)]
        for terms,value in zip(letters,target):
            errors.append(abs(sum(weight*mp.log(1-z/a) for a,weight in terms)-value))
    return max(errors)

def soft_rest_closed(t):
    """Soft finite function with the final electron at rest, m=1."""
    t=mp.mpf(t);ze=1-t/2
    if t==0:return mp.mpf(0)
    eta=mp.acosh(ze);v=mp.tanh(eta)
    J=(li2real(v)-li2real(-v))/2+eta*mp.log(2*mp.cosh(eta))-eta*eta/2-li2real(-mp.exp(-2*eta))/2-mp.pi**2/24
    b=eta/v
    return J/v-(b+1)/2

def soft_frame_change(s,t):
    """(soft_rest(sqrt(s)*epsilon)-soft_CM(epsilon))/alpha."""
    s=mp.mpf(s);t=mp.mpf(t);ze=1-t/2
    b=mp.acosh(ze)*ze/mp.sqrt(ze*ze-1) if t else mp.mpf(1)
    return 2/mp.pi*((b-1)*mp.log(s)/2+soft_closed(s,t)-soft_rest_closed(t))

def soft_frame_angular_check(s,c,order=100):
    """Independent angular check of the Lorentz-frame cutoff conversion."""
    import numpy as np
    s=float(s);c=float(c);E=(s+1)/(2*np.sqrt(s));q=(s-1)/(2*np.sqrt(s))
    ze=E*E-q*q*c
    nodes,weights=np.polynomial.legendre.leggauss(order)
    phi=np.pi*(nodes+1);phiw=np.pi*weights
    total=0.
    for mu,weight in zip(nodes,weights):
        den=E-q*mu
        denp=E-q*(c*mu+np.sqrt(1-c*c)*np.sqrt(1-mu*mu)*np.cos(phi))
        minusj=2*ze/(den*denp)-1/den**2-1/denp**2
        total+=weight*np.sum(phiw*minusj*np.log(np.sqrt(s)/denp))
    return total/(4*np.pi**2)
