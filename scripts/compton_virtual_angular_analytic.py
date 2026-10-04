"""Analytic elastic-angle integration of vertex and self-energy interference.

The integration engine reduces rational coefficients times GPLs by partial
fractions and integration by parts. It never calls numerical quadrature.
"""
from compton_artifacts import read_artifact, write_artifact
from functools import lru_cache
import json
from pathlib import Path
import sympy as S
import mpmath as mp
from compton_vertex_analytic import vertex_moments_mp

x,y=S.symbols('x y')
cz=S.Symbol('cz')
G=S.Function('G')

def g(word):return G(S.Tuple(*word),y) if word else S.S.One

@lru_cache(None)
def primitive(rational,word):
    rational=S.cancel(rational)
    if rational==0:return S.S.Zero
    apart=S.apart(rational,y)
    out=0
    for term in S.Add.make_args(apart):
        numerator,denominator=S.fraction(term)
        poly=S.Poly(denominator,y)
        if poly.degree()==0:
            for (power,),coefficient in S.Poly(term,y).terms():
                p=coefficient*y**(power+1)/S.Rational(power+1)
                out+=p*g(word)
                if word:out-=primitive(S.cancel(p/(y-word[0])),word[1:])
        else:
            roots=S.roots(poly.as_expr(),y)
            if len(roots)!=1:raise ValueError((term,roots))
            root,power=next(iter(roots.items()))
            coefficient=S.cancel(term*(y-root)**power)
            if power==1:
                out+=coefficient*g((root,)+word)
            else:
                p=-coefficient/(power-1)/(y-root)**(power-1)
                out+=p*g(word)
                if word:out-=primitive(S.cancel(p/(y-word[0])),word[1:])
    return S.expand(out)


def crossed_integrand():
    ell,d2=S.symbols('ell d2');A=1-y;D=-y
    P=[d2,d2+(1-1/A)*ell,2*A-2*d2+(-A+1/A)*ell-1,
       d2+(S.Rational(3,2)-2/A+1/(2*A*A))*ell+1/(2*A),
       S.Rational(5)*A/2-3*d2+(-A-S.Rational(3,2)+3/A-1/(2*A*A))*ell-1-1/(2*A)]
    bounds=[S.pi**2/6,S.pi**2/6,1-S.pi**2/3,S.Rational(1,2)+S.pi**2/6,1-S.pi**2/2]
    J=[(p-b)/D**power for p,b,power in zip(P,bounds,[1,1,2,1,2])]
    F0=-ell/A;F1=(1-y*F0)/A;F2=(S.Rational(1,2)-y*F1)/A
    params=dict(zip(['uJ00','uJ10','uJ11','uJ20','uJ21'],J))
    params.update(ua=-S.Rational(1,2)-2*A*F1+A*F2,ub=2+4*A*F1,cz=cz,constant=1)
    data=json.loads(read_artifact(Path(__file__).with_name('compton_virtual_trace_analytic.json')))
    expr=S.expand(sum(S.sympify(v)*params[k] for key,rec in data.items() if key.startswith('u_') for k,v in rec.items()))
    rD=S.factor(expr.coeff(d2));rL=S.factor(expr.coeff(ell))
    r0=S.factor(expr.subs({d2:0,ell:0})+S.pi**2*rD/6)
    return r0,rL,-rD


def derivative(expr):
    # Differentiate independent GPL atoms first, then apply their recursion.
    atoms=list(expr.atoms(G));dummy=S.symbols('h:'+str(len(atoms)))
    f=expr.xreplace(dict(zip(atoms,dummy)))
    result=S.diff(f,y)
    for atom,h in zip(atoms,dummy):
        word=tuple(atom.args[0]);result+=S.diff(f,h)*g(word[1:])/(y-word[0])
    return S.expand(result.xreplace(dict(zip(dummy,atoms))))


@lru_cache(None)
def direct_coefficients():
    data=json.loads(read_artifact(Path(__file__).with_name('compton_virtual_trace_analytic.json')))
    totals={}
    for name,record in data.items():
        if name.startswith('s_'):
            for key,value in record.items():
                totals[key]=totals.get(key,0)+S.sympify(value)
    integrated={}
    for key,value in totals.items():
        result=0
        for term in S.Add.make_args(S.expand(value)):
            power=term.as_powers_dict().get(y,0)
            coefficient=term/y**power
            result+=coefficient*(S.log(1+x) if power==-1 else
                (x**(power+1)-(x/(1+x))**(power+1))/(power+1))
        integrated[key]=S.lambdify(x,S.factor(result),'mpmath')
    return integrated


def crossed_primitive_mp(x,y,cz):
    ell=mp.log(y)
    # For y>1, imaginary parts on the upper lip cancel exactly in H010.
    H10=mp.pi**2/6-mp.polylog(2,1-y)
    H010=mp.re(2*mp.polylog(3,y)-ell*mp.polylog(2,y))
    Pregular=-4*(x+8)/y+(8*x+12)*y/x+(2*cz*(x*x-3*x-4)+2*x*x+10*x)/x
    Plog=Pregular*ell-2*(ell/(y-1) if y!=1 else 1)
    Q=-4*x/y**2+16/y+(2*x*x-18*x+20)/x
    R=4*(-2*cz+x-8)/y+(2*cz-12*x-16)*y/x
    return -2*(x*x+8)/x*ell**2+8*(x+1)/x*H010+Plog+Q*H10+R


def coefficient_closed(s,lam):
    """Coefficient of alpha^3/m^2 for the vertex+self-energy subset."""
    s=mp.mpf(s);lam=mp.mpf(lam);x=s-1;cz=-4-4*mp.log(lam)
    js=vertex_moments_mp(s)
    logB=mp.log(x)-mp.j*mp.pi
    F0=-logB/s;F1=(1+x*F0)/s;F2=(mp.mpf(1)/2+x*F1)/s
    values=dict(zip(['sJ00','sJ10','sJ11','sJ20','sJ21','sJ22'],js))
    values.update(sa=-mp.mpf(1)/2-2*s*F1+s*F2,sb=2+4*s*F1,cz=cz,constant=1)
    Is=sum(fun(x)*values[key] for key,fun in direct_coefficients().items())
    Iu=crossed_primitive_mp(x,x,cz)-crossed_primitive_mp(x,x/(1+x),cz)
    return mp.re(Is+Iu)/(2*x*x)


def main():
    r0,rL,rG=crossed_integrand()
    F=primitive(r0,())+primitive(rL,(0,))+primitive(rG,(1,0))
    for atom in list(F.atoms(G)):
        F=S.collect(F,atom)
    assert S.cancel(derivative(F)-(r0+rL*g((0,))+rG*g((1,0))))==0
    coefficients={str(atom):str(S.factor(F.coeff(atom))) for atom in sorted(F.atoms(G),key=str)}
    coefficients['constant']=str(S.factor(F.subs({atom:0 for atom in F.atoms(G)})))
    write_artifact(Path(__file__).with_suffix('.json') , json.dumps(coefficients,indent=2)+'\n')
    print(json.dumps(coefficients,indent=2))

if __name__=='__main__':main()
