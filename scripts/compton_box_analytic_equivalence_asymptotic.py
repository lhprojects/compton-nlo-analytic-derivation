"""Exact low-energy series of our independent virtual coefficient.

The GPL recursion is integrated as rational Puiseux-log series at X=0.
No numerical fitting and no published asymptotic coefficients are used.
"""
from compton_artifacts import read_artifact, write_artifact
from functools import lru_cache
from collections import defaultdict
from pathlib import Path
import json
import sympy as S
X,Y,ll,L,t=S.symbols('X Y ll L t')
ROOT=Path(__file__).resolve().parents[1]
MAX=16


@lru_cache(None)
def gseries(word):
    if not word:return {(0,0):S.S.One}
    if all(a==0 for a in word):return {(0,len(word)):1/S.factorial(len(word))}
    a=word[0];tail=gseries(word[1:]);out=defaultdict(lambda:S.S.Zero)
    if a==0:kernel={-2:S.S.One}
    elif a=='r':kernel={2*j-1:S.binomial(-S.Rational(1,2),j)/(2*4**j) for j in range(MAX//2+2)}
    else:kernel={2*j:-S.Rational(1,int(a)**(j+1)) for j in range(MAX//2+2)}
    for (n,k),c in tail.items():
        for p,d in kernel.items():
            power=n+p+2
            if power>MAX:continue
            if power==0:out[0,k+1]+=c*d/(k+1)
            else:
                exponent=S.Rational(power,2)
                for j in range(k+1):out[power,k-j]+=c*d*(-1)**j*S.factorial(k)/S.factorial(k-j)/exponent**(j+1)
    return {key:S.factor(value) for key,value in out.items() if value!=0}


def high_partial(rows,return_values=False):
    p=S.pi;Z=S.zeta(3);a,b,c,d,e,f=S.symbols('K_a K_b K_c K_d K_e K_f')
    values={():1,(-1,):L,(0,):L,('r',-1):L**2/2-p*p/18,
      (-1,-1):L**2/2,(-1,0):L**2/2-p*p/6,(-2,-1):L**2/2-p*p/12,
      (0,-1):L**2/2+p*p/6,(1,0):L**2/2+p*p/3,
      ('r','r',-1):L**3/6-p*p*L/18+2*Z/3,
      ('r',-1,-1):L**3/6+a,('r',-1,0):L**3/6-p*p*L/6+b,
      ('r',0,-1):L**3/6+p*p*L/6+c,('r',-2,-1):L**3/6-p*p*L/12+d,
      (-4,'r',-1):L**3/6-p*p*L/18+e,(0,'r',-1):L**3/6-p*p*L/18+f,
      (-1,-1,-1):L**3/6,(-1,-1,0):L**3/6-p*p*L/6+Z,
      (-1,-2,-1):L**3/6-p*p*L/12+3*Z/2,(-1,0,-1):L**3/6+p*p*L/6-2*Z,
      (0,-1,-1):L**3/6+Z,(0,-1,0):L**3/6-p*p*L/6,(0,1,0):L**3/6+p*p*L/3}
    if return_values:return values
    answer=0
    for row in rows:
        word=tuple(a if a=='r' else int(a) for a in row['word'])
        coefficient=S.limit(X*S.sympify(row['coefficient']).subs(Y,1),X,S.oo)
        answer+=coefficient*values[word]
    return S.collect(S.expand(answer),L)/X


def build():
    folder=ROOT
    rows=json.loads(read_artifact(folder/'compton_box_analytic_equivalence_virtual_x.json'))['coefficients']
    total=defaultdict(lambda:S.S.Zero)
    for row in rows:
        word=tuple(a if a=='r' else int(a) for a in row['word'])
        coefficient=S.sympify(row['coefficient']).subs({X:t*t,Y:t/S.sqrt(t*t+4)})
        pref=S.series(coefficient,t,0,9).removeO().expand()
        terms=[]
        for term in S.Add.make_args(pref):
            n=int(term.as_powers_dict().get(t,0));terms.append((n,term/t**n))
        for n,c in terms:
            for (p,k),d in gseries(word).items():
                if n+p<=6:total[n+p,k]+=c*d
    total={key:S.factor(c) for key,c in total.items() if S.factor(c)!=0}
    assert all(n%2==0 for n,k in total)
    expression=S.expand(sum(c*X**(S.Rational(n,2))*L**k for (n,k),c in total.items()))
    print('Independent CV low X:',S.collect(expression,X),flush=True)
    # The Born total follows directly from the Klein--Nishina formula.
    r=X/2
    born=2*((1+r)/r**2*(2*(1+r)/(1+2*r)-S.log(1+2*r)/r)+S.log(1+2*r)/(2*r)-(1+3*r)/(1+2*r)**2)
    born_low=S.series(born,X,0,4)
    inv=S.symbols('inv',positive=True)
    born_high=S.series(born.subs(X,1/inv),inv,0,3).removeO().subs(inv,1/X).expand()
    print('Born low:',born_low,flush=True)
    print('Born high:',born_high,flush=True)
    out=dict(source='Independent virtual X-word coefficients, exact iterated-kernel series',normalization='sigma_virtual=alpha_T^3/m^2*CV, X=s/m^2-1',virtual_low=str(S.collect(expression,X)),log_symbols={'L':'log(X)','ll':'log(lambda/m)'},virtual_remainder='O(X^4*log(X)) at fixed lambda/m',born_normalization='sigma_Born=pi*alpha_T^2/m^2*B(X)',born_low=str(born_low),born_high=str(born_high),born_high_remainder='O(log(X)/X^3)')
    before=high_partial(rows)
    constants=json.loads(read_artifact(folder/'total_energy_high.json'))['constants']
    evaluated=S.expand(before.subs({S.Symbol(k):S.sympify(v) for k,v in constants.items()}))
    assert not evaluated.has(S.Symbol('C'))
    out['virtual_high_before_constants']=str(before)
    out['virtual_high']=str(S.collect(evaluated*X,L)/X)
    out['high_status']='Complete, using the independently evaluated endpoints in total_energy_high.json.'
    write_artifact(folder/'compton_box_analytic_equivalence_asymptotic.json' , json.dumps(out,indent=2)+'\n')


if __name__=='__main__':build()
