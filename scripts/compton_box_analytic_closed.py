"""Closed analytic Compton spanning contribution integrated over the angle.

The public function performs finite Laurent algebra and GPL endpoint
 evaluation. No angular or loop-parameter quadrature is called. Positive
crossed-channel GPL poles are passed above, as fixed by the Feynman +i0.
"""
from compton_artifacts import read_artifact
import json,math
from functools import lru_cache
from pathlib import Path
from collections import defaultdict
import sympy as S
import mpmath as mp
from compton_box_analytic_angle import x,v,ll,lx,ds,du
from compton_box_analytic_continuation import GPLContinuation

LOCAL={'x':x,'v':v,'ll':ll,'lx':lx,'ds':ds,'du':du}

@lru_cache(None)
def _specifications(channel):
    data=json.loads(read_artifact(Path(__file__).with_name('compton_box_analytic_angle.json')))[channel]['integrand']
    specifications=[]
    for row in data:
        word=tuple(S.sympify(a,locals=LOCAL) for a in row['word'])
        expr=S.sympify(row['coefficient'],locals=LOCAL)
        p,q=S.fraction(S.cancel(expr));P,Q=S.Poly(p,v),S.Poly(q,v)
        numerator=S.lambdify((x,ll,lx,ds,du),list(reversed(P.all_coeffs())),'mpmath')
        denominator=S.lambdify((x,ll,lx,ds,du),list(reversed(Q.all_coeffs())),'mpmath')
        poles=[]
        for f,m in S.factor_list(q,v)[1]:
            poly=S.Poly(f,v)
            if poly.degree()==1:
                a,b=poly.all_coeffs();roots=[-b/a]
            elif poly.degree()==2:
                a,b,c=poly.all_coeffs();rad=S.sqrt(b*b-4*a*c)
                roots=[(-b+rad)/(2*a),(-b-rad)/(2*a)]
            else:raise ArithmeticError('Unexpected GPL pole polynomial')
            poles.extend((S.simplify(a),int(m)) for a in roots)
        assert sum(m for _,m in poles)==Q.degree()
        specifications.append((word,numerator,denominator,poles))
    return specifications


def _taylor_polynomial(coefficients,a,n):
    """Coefficient of (v-a)^n in a polynomial, a finite exact formula."""
    return mp.fsum(coefficients[k]*math.comb(k,n)*a**(k-n) for k in range(n,len(coefficients)))


class _Primitive:
    def __init__(self,evaluator):
        self.evaluator=evaluator
        self.terms=defaultdict(lambda:mp.mpc(0))

    def add(self,coefficient,word,kind='poly',power=0,pole=None):
        if coefficient:
            self.terms[(word,kind,power,pole)]+=coefficient

    def integrate_poly(self,coefficient,n,word):
        coefficient/=n+1
        self.add(coefficient,word,'poly',n+1)
        if not word:return
        b=word[0];bb=self.evaluator.letter(b)
        # v^(n+1)/(v-b) is divided exactly into a polynomial and one pole.
        for k in range(n+1):
            self.integrate_poly(-coefficient*bb**k,n-k,word[1:])
        self.integrate_pole(-coefficient*bb**(n+1),b,1,word[1:])

    def integrate_pole(self,coefficient,a,j,word):
        if j==1:
            self.add(coefficient,(a,)+word)
            return
        primitive=-coefficient/(j-1)
        self.add(primitive,word,'pole',j-1,a)
        if not word:return
        b=word[0];aa=self.evaluator.letter(a);bb=self.evaluator.letter(b)
        n=j-1
        if aa==bb:
            self.integrate_pole(-primitive,a,n+1,word[1:])
        else:
            delta=aa-bb
            for k in range(1,n+1):
                factor=(-1)**(n-k)/delta**(n-k+1)
                self.integrate_pole(-primitive*factor,a,k,word[1:])
            self.integrate_pole(-primitive*(-1)**n/delta**n,b,1,word[1:])

    def integrate_rational(self,word,p,q,poles):
        # Polynomial division has a known finite number of steps.
        p=[mp.mpc(c) for c in p];q=[mp.mpc(c) for c in q]
        if len(p)>=len(q):
            for n in range(len(p)-len(q),-1,-1):
                c=p[n+len(q)-1]/q[-1]
                self.integrate_poly(c,n,word)
                for k,qk in enumerate(q):p[n+k]-=c*qk
            p=p[:len(q)-1]
        merged={}
        for a,m in poles:
            aa=self.evaluator.letter(a)
            if aa in merged:merged[aa][1]+=m
            else:merged[aa]=[S.S.Zero if aa==0 else a,m]
        for aa,(a,m) in merged.items():
            qj=[_taylor_polynomial(q,aa,m+j) for j in range(m)]
            pj=[_taylor_polynomial(p,aa,j) for j in range(m)]
            h=[]
            for j in range(m):
                h.append((pj[j]-mp.fsum(qj[k]*h[j-k] for k in range(1,j+1)))/qj[0])
                self.integrate_pole(h[-1],a,m-j,word)

    def endpoint(self,upper):
        words={key[0] for key in self.terms}
        values,steps=self.evaluator.values(list(words),upper)
        result=mp.mpc(0)
        lower=mp.mpc(0)
        for (word,kind,power,a),c in self.terms.items():
            if kind=='poly':
                rational=upper**power
                low=1 if power==0 and not word else 0
            else:
                aa=self.evaluator.letter(a)
                rational=(upper-aa)**(-power)
                if aa!=0:low=(-aa)**(-power) if not word else 0
                elif not word or word==(S.S.Zero,):low=0
                else:low=self.evaluator.coefficients(word)[power]
            result+=c*rational*values[word]
            lower+=c*low
        return result-lower,steps


def box_angular_components(s,lam,dps=50):
    """Return the two unnormalized analytic dy integrals, including branches."""
    with mp.workdps(dps+30):
        ss=mp.mpf(s);xx=ss-1;lambda_=mp.mpf(lam)
        if xx<=0 or lambda_<=0:raise ValueError('Require s>1 and lambda>0 in m=1 units')
        scale=max(mp.mpf(1),xx,1/xx,1/abs(xx-1) if xx!=1 else 1)
        guard=35+8*int(mp.ceil(mp.log10(scale)))
        mp.mp.dps=dps+guard
        ss=mp.mpf(s);xx=ss-1;lambda_=mp.mpf(lam)
        values=(xx,mp.log(lambda_),mp.log(xx),
                -mp.re(mp.polylog(2,1+xx)-mp.pi**2/6)/xx,
                mp.polylog(2,1-xx)-mp.pi**2/6)
        order=max(140,int((mp.mp.dps+15)/.49)+1)
        result=[]
        for channel in ['s','u']:
            evaluator=GPLContinuation(xx,order)
            primitive=_Primitive(evaluator)
            for word,numerator,denominator,poles in _specifications(channel):
                word=tuple(S.S.Zero if evaluator.letter(a)==0 else a for a in word)
                primitive.integrate_rational(word,numerator(*values),denominator(*values),poles)
            value,_=primitive.endpoint(xx/(xx+2))
            tolerance=mp.power(10,-dps)*max(1,abs(value))
            if abs(mp.im(value))>tolerance:
                raise ArithmeticError('Causal GPL imaginary parts did not cancel: '+mp.nstr(value,20))
            result.append(+value)
        return tuple(result)


def box_angular_closed(s,lam,dps=50):
    """Return C_box=(J_s+J_u)/(2*(s-1)^2), by closed GPL evaluation only."""
    with mp.workdps(dps+10):
        pieces=box_angular_components(s,lam,dps)
        return +mp.re(sum(pieces)/(2*(mp.mpf(s)-1)**2))

if __name__=='__main__':
    for point in ['1.6','10']:
        print(point,mp.nstr(box_angular_closed(point,'1e-5',50),52),flush=True)
