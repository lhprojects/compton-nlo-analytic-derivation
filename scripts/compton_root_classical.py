"""Classical primitives derived from the rationalized square-root kernel.
No reference cross section is used. Quadrature occurs only in the check.
"""
from compton_artifacts import write_artifact
import mpmath as mp
import json
from pathlib import Path

def J(z,rho):
    if z == 0: return mp.mpf(0)
    if z == rho == 1: return -2*mp.zeta(3)
    l=mp.log(z)
    if rho == 1 and mp.im(z) == 0 and mp.re(z)>1:
        return l*l*mp.log(z-1)+2*l*mp.re(mp.polylog(2,z))-2*mp.re(mp.polylog(3,z))
    return l*l*mp.log(1-z/rho)+2*l*mp.polylog(2,z/rho)-2*mp.polylog(3,z/rho)

def primitive(q):
    a=mp.exp(mp.j*mp.pi/3); b=a.conjugate()
    A=1-a*q; B=1-b*q; R=A/B; l=mp.log(q)
    return (2*J(A,1)+2*J(B,1)-J(R,1)+J(R,a/b)
            +2*l*(mp.polylog(2,a*q)+mp.polylog(2,b*q))
            -2*(mp.polylog(3,a*q)+mp.polylog(3,b*q))+l**3/3)

def root_m1_m1(X):
    X=mp.mpf(X)
    q=2/(X+2+mp.sqrt(X*(X+4)))
    return -(primitive(q)-primitive(mp.mpf(1)))/2

def mixed_pole(q,a):
    """Primitive of log(q) log(1-a*q)/(q-1), including its finite endpoint."""
    if q == 1:
        b=1-a
        return (-2*mp.zeta(3)+mp.log(b)**2*mp.log(-a)
                +4*mp.zeta(2)*mp.log(b)+J(1/b,-1/a))/2
    B=1-a*q; R=q/B
    return (J(q,1)+J(B,1-a)-J(R,1/(1-a))+J(R,-1/a))/2

def mixed_zero(q,a):
    """Primitive of log(1-q) log(1-a*q)/q."""
    B=1-a*q; R=(1-q)/B
    return (J(1-q,1)+J(B,1)-J(R,1)+J(R,1/a))/2

def primitive_r0m1(q):
    a=mp.exp(mp.j*mp.pi/3); b=a.conjugate(); l=mp.log(q)
    X=(1-q)**2/q
    N=sum(-l*mp.polylog(2,c*q)+mp.polylog(3,c*q) for c in [a,b])
    return l*mp.polylog(2,-X)+2*sum(mixed_pole(q,c) for c in [a,b])-N-2*J(q,1)+l**3/3

def primitive_shuffle(q):
    a=mp.exp(mp.j*mp.pi/3); b=a.conjugate(); l=mp.log(q)
    N=sum(-l*mp.polylog(2,c*q)+mp.polylog(3,c*q) for c in [a,b])
    return (-2*sum(mixed_zero(q,c) for c in [a,b])+N
            +2*(-l*mp.polylog(2,q)+mp.polylog(3,q))-l**3/3)

def root_0_m1(X):
    X=mp.mpf(X); q=2/(X+2+mp.sqrt(X*(X+4)))
    return primitive_r0m1(q)-primitive_r0m1(mp.mpf(1))

def root_m1_0(X):
    X=mp.mpf(X); q=2/(X+2+mp.sqrt(X*(X+4)))
    return primitive_shuffle(q)-primitive_shuffle(mp.mpf(1))-root_0_m1(X)

def mixed_general(q,a,c):
    """Primitive of log(q) log(1-a*q)/(q-c), for c=+i or -i."""
    B=1-a*q; R=q/B
    return (J(q,c)+J(B,1-a*c)-J(R,c/(1-a*c))+J(R,-1/a))/2

def mixed_plus(q,a):
    """Primitive of log(1+q) log(1-a*q)/q."""
    B=1-a*q; R=(1+q)/B
    return (J(1+q,1)+J(B,1)-J(R,1)+J(R,-1/a))/2

def groot(X):
    q=2/(X+2+mp.sqrt(X*(X+4)))
    return mp.log(q)**2/2+mp.polylog(2,-q**3)/3-mp.polylog(2,-q)-mp.pi**2/18

def zero_root_m1(X):
    X=mp.mpf(X); q=2/(X+2+mp.sqrt(X*(X+4)))
    return mp.log(X)*groot(X)-primitive_shuffle(q)+primitive_shuffle(mp.mpf(1))

def primitive_outer4(q):
    a=mp.exp(mp.j*mp.pi/3); b=a.conjugate(); l=mp.log(q)
    N=sum(-l*mp.polylog(2,c*q)+mp.polylog(3,c*q) for c in [a,b])
    return (2*sum(mixed_plus(q,c) for c in [a,b])-N
            -2*(-l*mp.polylog(2,-q)+mp.polylog(3,-q))+l**3/3
            +mp.log(4)*sum(mp.polylog(2,c*q) for c in [a,b])+mp.log(4)*l*l/2)

def m4_root_m1(X):
    X=mp.mpf(X); q=2/(X+2+mp.sqrt(X*(X+4)))
    return mp.log((X+4)/4)*groot(X)+primitive_outer4(q)-primitive_outer4(mp.mpf(1))

def primitive_rm2m1(q):
    a=mp.exp(mp.j*mp.pi/3); b=a.conjugate(); l=mp.log(q); X=(1-q)**2/q
    D=mp.polylog(2,-1-X)+mp.log(1+X)*mp.log(2+X)+mp.pi**2/12
    N=sum(-l*mp.polylog(2,c*q)+mp.polylog(3,c*q) for c in [a,b])
    return (-l*D+sum(sum(mixed_general(q,d,c) for d in [a,b])-J(q,c)
                    for c in [mp.j,-mp.j])-N+l**3/3)

def root_m2_m1(X):
    X=mp.mpf(X); q=2/(X+2+mp.sqrt(X*(X+4)))
    return primitive_rm2m1(q)-primitive_rm2m1(mp.mpf(1))

def all_roots(X):
    """Evaluate the eight root words with shared endpoint differences.

    Use the caller's precision and the original complex primitives, including
    their original branch choices and endpoint values. No constant is dropped.
    The shuffle identity supplies G(0,r,-1) from the same differences used for
    G(r,0,-1) and G(r,-1,0).
    """
    X=mp.mpf(X)
    if X<=0:
        raise ValueError('X must be positive')
    q=2/(X+2+mp.sqrt(X*(X+4))); l=mp.log(q); one=mp.mpf(1)
    g=groot(X)
    d0=primitive_r0m1(q)-primitive_r0m1(one)
    ds=primitive_shuffle(q)-primitive_shuffle(one)
    return {
        ('r','-1'):g,
        ('0','r','-1'):mp.log(X)*g-ds,
        ('-4','r','-1'):mp.log((X+4)/4)*g+primitive_outer4(q)-primitive_outer4(one),
        ('r','-1','-1'):-(primitive(q)-primitive(one))/2,
        ('r','-2','-1'):primitive_rm2m1(q)-primitive_rm2m1(one),
        ('r','0','-1'):d0,
        ('r','-1','0'):ds-d0,
        ('r','r','-1'):-l**3/6-mp.polylog(3,-q**3)/9+mp.polylog(3,-q)
            +mp.pi**2*l/18+2*mp.zeta(3)/3,
    }

if __name__=='__main__':
    mp.mp.dps=70
    rows=[]
    for xs in ['0.000001','0.01','0.6','1','10','10000']:
        x=mp.mpf(xs); q=2/(x+2+mp.sqrt(x*(x+4)))
        value=root_m1_m1(x)
        # t=4 sinh(y/2)^2 removes the original root kernel exactly.
        top=2*mp.asinh(mp.sqrt(x)/2)
        check=mp.quad(lambda y: mp.log(1+4*mp.sinh(y/2)**2)**2/2,[0,top])
        deriv=mp.diff(primitive,q)-mp.log((1-q+q*q)/q)**2/q
        rows.append({'word':'r,-1,-1','X':xs,'classical':mp.nstr(mp.re(value),60),
                     'quadrature_difference':mp.nstr(abs(value-check),8),
                     'derivative_residual':mp.nstr(abs(deriv),8),
                     'imaginary_residual':mp.nstr(abs(mp.im(value)),8)})
        assert abs(value-check)<mp.mpf('1e-55')*max(1,abs(check)), rows[-1]
        assert abs(deriv)<mp.mpf('1e-55')*max(1,abs(mp.log((1-q+q*q)/q)**2/q)),rows[-1]
    for xs in ['0.000001','0.01','0.6','1','10','10000']:
        x=mp.mpf(xs); top=2*mp.asinh(mp.sqrt(x)/2)
        for name,closed,inner in [
            ('r,0,-1',root_0_m1,lambda t: -mp.polylog(2,-t)),
            ('r,-1,0',root_m1_0,lambda t: mp.log(1+t)*mp.log(t)+mp.polylog(2,-t))]:
            value=closed(x)
            check=mp.quad(lambda y: inner(4*mp.sinh(y/2)**2),[0,top])
            error=abs(value-check)
            row={'word':name,'X':xs,'classical':mp.nstr(mp.re(value),60),
                 'quadrature_difference':mp.nstr(error,8),'imaginary_residual':mp.nstr(abs(mp.im(value)),8)}
            rows.append(row)
            assert error<mp.mpf('1e-55')*max(1,abs(check)),row
    for xs in ['0.000001','0.01','0.6','1','10','10000']:
        x=mp.mpf(xs); top=2*mp.asinh(mp.sqrt(x)/2)
        tests=[('0,r,-1',zero_root_m1,lambda y: groot(4*mp.sinh(y/2)**2)/mp.tanh(y/2)),
               ('-4,r,-1',m4_root_m1,lambda y: groot(4*mp.sinh(y/2)**2)*mp.tanh(y/2)),
               ('r,-2,-1',root_m2_m1,lambda y: mp.polylog(2,-1-4*mp.sinh(y/2)**2)
                +mp.log(1+4*mp.sinh(y/2)**2)*mp.log(2+4*mp.sinh(y/2)**2)+mp.pi**2/12)]
        for name,closed,kernel in tests:
            value=closed(x); check=mp.quad(kernel,[0,top]); error=abs(value-check)
            row={'word':name,'X':xs,'classical':mp.nstr(mp.re(value),60),
                 'quadrature_difference':mp.nstr(error,8),'imaginary_residual':mp.nstr(abs(mp.im(value)),8)}
            rows.append(row)
            assert error<mp.mpf('1e-55')*max(1,abs(check)),row
    derivatives=[]
    for xs in ['0.000001','0.6','1','10','10000']:
        x=mp.mpf(xs); root=mp.sqrt(x*(x+4))
        tests=[('r,-1,-1',root_m1_m1,mp.log(1+x)**2/(2*root)),
               ('r,0,-1',root_0_m1,-mp.polylog(2,-x)/root),
               ('r,-1,0',root_m1_0,(mp.log(1+x)*mp.log(x)+mp.polylog(2,-x))/root),
               ('0,r,-1',zero_root_m1,groot(x)/x),
               ('-4,r,-1',m4_root_m1,groot(x)/(x+4)),
               ('r,-2,-1',root_m2_m1,(mp.polylog(2,-1-x)+mp.log(1+x)*mp.log(2+x)+mp.pi**2/12)/root)]
        for name,closed,expected in tests:
            residual=abs(mp.diff(closed,x)-expected)
            assert residual<mp.mpf('1e-55')*max(1,abs(expected)),(name,x,residual)
            derivatives.append({'word':name,'X':xs,'residual':mp.nstr(residual,8)})
    path=Path(__file__).resolve().parents[1]/'root_classical_check.json'
    write_artifact(path , json.dumps({'dps':70,'method':'classical endpoints versus defining one-dimensional integrals; quadrature is only a check','rows':rows,'derivative_checks':derivatives},indent=2)+'\n')
    print(json.dumps(rows,indent=2))
