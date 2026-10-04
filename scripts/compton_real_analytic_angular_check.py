"""Numerical checks ONLY of the analytically completed direction integral."""
from compton_artifacts import read_artifact
from pathlib import Path
import json
import numpy as np
import mpmath as mp
import sympy as S
from numpy.polynomial.legendre import leggauss
H,x,y,U,V=S.symbols('H x y U V')
records=json.loads(read_artifact(Path(__file__).with_name('compton_real_analytic_angles.json')))
coeff=[(tuple(r['powers']),S.lambdify((H,x,y),S.sympify(r['coefficient']),'numpy')) for r in records]


def average_monomial_double(e,a):
    """Finite-precision diagnostic, unreliable for small Gram determinants."""
    neg=[i for i,n in enumerate(e) if n<0];pos=[i for i,n in enumerate(e) if n>0]
    if not neg:
        if not pos:return 1.
        assert sum(e)==1
        return a[pos[0]][0]
    i=neg[0];av=a[i];A=av[0]**2-av[1:]@av[1:];am=np.linalg.norm(av[1:])
    L=np.log((av[0]+am)/(av[0]-am))/(2*am)
    if len(neg)==1:
        if not pos:return L if e[i]==-1 else 1/A
        j=pos[0];bv=a[j];B=bv[0]**2-bv[1:]@bv[1:]
        d=av[0]*bv[0]-av[1:]@bv[1:]
        h=(av[0]*bv[0]-d)/(am*am);g=bv[0]-h*av[0]
        if e[j]==1:return g*L+h if e[i]==-1 else g/A+h*L
        assert e[j]==2 and e[i]==-1
        perp=bv[0]**2-B-(av[0]*bv[0]-d)**2/(am*am)
        c0=g*g-perp*A/(2*am*am);c1=2*g*h+perp*av[0]/(am*am);c2=h*h-perp/(2*am*am)
        return c0*L+c1+c2*av[0]
    j=neg[1];bv=a[j];B=bv[0]**2-bv[1:]@bv[1:];d=av[0]*bv[0]-av[1:]@bv[1:]
    D=d*d-A*B;root=np.sqrt(D);J=np.log((d+root)/(d-root))/(2*root);at,bt=av[0],bv[0]
    if e[i]==-1 and e[j]==-1:return J
    if e[i]==-2 and e[j]==-1:return ((bt*d-at*B)*J+at*d/A-bt)/D
    if e[i]==-1 and e[j]==-2:return ((at*d-bt*A)*J+bt*d/B-at)/D
    assert e[i]==e[j]==-2
    return ((2*at*bt*(A*B+2*d*d)-d*(D+3*A*bt*bt+3*B*at*at))*J+D+(2*A*B+d*d)*(at*at/A+bt*bt/B)-6*at*bt*d)/(D*D)


def vectors(h,xx,yy):
    s=h+1;z=xx+yy-h;c=1-2*s*z/(xx*yy)
    a=np.array([xx*(s+1)/(2*s),0,0,-xx*h/(2*s)])
    b=np.array([yy*(s+1)/(2*s),-yy*h*np.sqrt(1-c*c)/(2*s),0,-yy*h*c/(2*s)])
    w=a+b-np.array([z,0,0,0]);return [a,b,w]


def closed_double(h,xx,yy):
    """Double-precision diagnostic only, not a reference near threshold."""
    a=vectors(h,xx,yy)
    return sum(f(h,xx,yy)*average_monomial_double(e,a) for e,f in coeff)


def quadrature_double(h,xx,yy,n=100):
    """Double-precision diagnostic only; coefficients can cancel severely."""
    a=vectors(h,xx,yy);c,w=leggauss(n);phi=(np.arange(2*n)+.5)*np.pi/n
    dirs=np.stack(np.broadcast_arrays(np.sqrt(1-c[:,None]**2)*np.cos(phi),np.sqrt(1-c[:,None]**2)*np.sin(phi),c[:,None]),axis=-1)
    denominators=[v[0]-dirs@v[1:] for v in a]
    value=sum(f(h,xx,yy)*Sprod([d**p for d,p in zip(denominators,e)]) for e,f in coeff)
    return np.sum(w[:,None]*value)/(4*n)

def Sprod(vals):
    a=1
    for v in vals:a=a*v
    return a


coeff_mp=[(tuple(r['powers']),S.lambdify((H,x,y),S.sympify(r['coefficient']),'mpmath')) for r in records]

def _spatial(a,b):return sum(a[j]*b[j] for j in range(1,4))
def average_monomial(e,a):
 neg=[i for i,n in enumerate(e) if n<0];pos=[i for i,n in enumerate(e) if n>0]
 if not neg:return 1 if not pos else a[pos[0]][0]
 i=neg[0];av=a[i];A=av[0]**2-_spatial(av,av);am=mp.sqrt(_spatial(av,av));L=mp.log((av[0]+am)/(av[0]-am))/(2*am)
 if len(neg)==1:
  if not pos:return L if e[i]==-1 else 1/A
  j=pos[0];bv=a[j];B=bv[0]**2-_spatial(bv,bv);d=av[0]*bv[0]-_spatial(av,bv);h=(av[0]*bv[0]-d)/(am*am);g=bv[0]-h*av[0]
  if e[j]==1:return g*L+h if e[i]==-1 else g/A+h*L
  perp=bv[0]**2-B-(av[0]*bv[0]-d)**2/(am*am);c0=g*g-perp*A/(2*am*am);c1=2*g*h+perp*av[0]/(am*am);c2=h*h-perp/(2*am*am)
  return c0*L+c1+c2*av[0]
 j=neg[1];bv=a[j];B=bv[0]**2-_spatial(bv,bv);d=av[0]*bv[0]-_spatial(av,bv);D=d*d-A*B;root=mp.sqrt(D);J=mp.log((d+root)/(d-root))/(2*root);at,bt=av[0],bv[0]
 if e[i]==e[j]==-1:return J
 if e[i]==-2 and e[j]==-1:return ((bt*d-at*B)*J+at*d/A-bt)/D
 if e[i]==-1 and e[j]==-2:return ((at*d-bt*A)*J+bt*d/B-at)/D
 return ((2*at*bt*(A*B+2*d*d)-d*(D+3*A*bt*bt+3*B*at*at))*J+D+(2*A*B+d*d)*(at*at/A+bt*bt/B)-6*at*bt*d)/(D*D)

def vectors_mp(h,xx,yy):
    s=h+1;z=xx+yy-h;c=1-2*s*z/(xx*yy)
    a=[xx*(s+1)/(2*s),0,0,-xx*h/(2*s)]
    b=[yy*(s+1)/(2*s),-yy*h*mp.sqrt(1-c*c)/(2*s),0,-yy*h*c/(2*s)]
    return [a,b,[a[j]+b[j]-(z if j==0 else 0) for j in range(4)]]


def _closed_mp(h,xx,yy):
    a=vectors_mp(h,xx,yy)
    return mp.fsum(f(h,xx,yy)*average_monomial(e,a) for e,f in coeff_mp)


def closed(h,xx,yy,dps=50):
    """Controlled-precision master sum, checked at two working precisions.

    Near threshold, small Gram determinants and cancellation in the full
    sum make the double diagnostic unreliable. Repetition verifies the
    numerical evaluation, without choosing any coefficient or branch.
    """
    if not isinstance(dps,int) or dps<1:
        raise ValueError('dps must be a positive integer')
    with mp.workdps(dps+30):
        energy=mp.mpf(h)
        if not mp.isfinite(energy) or energy<=0:
            raise ValueError('H must be finite and positive')
        guards=30+int(mp.ceil(8*abs(mp.log10(energy))))
    previous=None
    for extra in (0,30,60):
        with mp.workdps(dps+guards+extra):
            hh,xxx,yyy=map(mp.mpf,(h,xx,yy))
            s=hh+1
            if not (0<xxx<hh and hh-xxx<yyy<s*(hh-xxx)/(s-xxx)):
                raise ValueError('x,y must lie inside the physical Dalitz region')
            value=_closed_mp(hh,xxx,yyy)
            if previous is not None and mp.isfinite(value) and mp.isfinite(previous):
                tolerance=mp.mpf(10)**(-dps)*max(abs(value),abs(previous))
                if abs(value-previous)<=tolerance:
                    return +value
            previous=value
    raise ArithmeticError('direction master sum did not converge to the requested precision')


def quadrature(h,xx,yy,n=20,dps=60):
    """Independent sphere quadrature with arbitrary-precision nodes and sums.

    dps controls arithmetic, and n controls angular integration refinement.
    Increasing dps alone does not bound the quadrature truncation error.
    """
    with mp.workdps(dps):
        hh,xxx,yyy=map(mp.mpf,(h,xx,yy));a=vectors_mp(hh,xxx,yyy)
        cv=[(e,f(hh,xxx,yyy)) for e,f in coeff_mp]
        nodes,weights=mp.gauss_quadrature(n,'legendre');rows=[]
        for j in range(n):
            ct=nodes[j];st=mp.sqrt(1-ct*ct);samples=[]
            for k in range(2*n):
                phi=(k+mp.mpf('.5'))*mp.pi/n
                direction=[st*mp.cos(phi),st*mp.sin(phi),ct]
                denominators=[v[0]-sum(v[l+1]*direction[l] for l in range(3)) for v in a]
                samples.append(mp.fsum(co*mp.fprod(v**power for v,power in zip(denominators,e)) for e,co in cv))
            rows.append(weights[j]*mp.fsum(samples)/(4*n))
        return +mp.fsum(rows)

if __name__=='__main__':
    for h,xx in [(.6,.2),(3.,1.3),(10.,3.),(40.,10.)]:
        lo=h-xx;hi=(h+1)*(h-xx)/(h+1-xx);yy=.4*lo+.6*hi
        v=closed(h,xx,yy);q=quadrature_double(h,xx,yy,160)
        print(h,xx,yy,v,q,v/q-1)
        assert abs(v/q-1)<2e-8

    with mp.workdps(80):
        h=mp.mpf('.02');xx=mp.mpf('.006')
        lo=h-xx;hi=(h+1)*(h-xx)/(h+1-xx);yy=mp.mpf('.4')*lo+mp.mpf('.6')*hi
        reference=closed(h,xx,yy,dps=50)
        coarse=quadrature(h,xx,yy,n=12,dps=80)
        fine=quadrature(h,xx,yy,n=20,dps=80)
        for value in (coarse,fine):
            assert abs(value-reference)<mp.mpf('1e-30')*abs(reference)
        print('H=.02 high-precision direction sum / independent quadrature',
              mp.nstr(reference,40),'relative difference',mp.nstr(abs(fine-reference)/abs(reference),10),flush=True)
