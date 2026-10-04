"""Exact soft residues and z-cut matching for the remaining real sectors.

Every residue is extracted algebraically from the independent angular trace.
No fit, published total, or numerical integration is used.
"""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
import json
import sympy as S
H,x,z,r,e=S.symbols('H x z r e')


def leading(expr,second=False):
    num,den=S.fraction(expr)
    change={x:H+e*(r-1) if second else e,z:e*r}
    def first(poly):
        P=S.Poly(S.expand(poly.xreplace(change)),e)
        for (power,),c in reversed(P.terms()):
            c=S.factor(c)
            if c!=0:return power,c
        return 1000,S.S.Zero
    pn,cn=first(num);pd,cd=first(den)
    return pn-pd,S.factor(cn/cd)


def main():
    data=json.loads(read_artifact(Path(__file__).with_name('compton_real_analytic_regular.json')))
    s=H+1;y=H-x+z
    sectors={'rational':S.sympify(data['rational_regular']),
             'L0':S.factor(S.sympify(data['L0'])*s/(H*x)),
             'L1':S.factor(S.sympify(data['L1'])*s/(H*y)),
             'L2':S.sympify(data['L2'])}
    output={}
    for name,expr in sectors.items():
        values=[]
        for second in [False,True]:
            power,value=leading(expr,second)
            print(name,'Y' if second else 'X',power,str(value),flush=True)
            values.append(dict(power=power,coefficient=str(value)))
        output[name]=values
        write_artifact(Path(__file__).with_suffix('.json') , json.dumps(output,indent=2)+'\n')




def integrate_residue(expr):
    """Return exact integral f(r) and f(r)*log(r), 0<r<H/(H+1).

    Only poles r=1 and r=-4/H occur; neither lies in this interval.
    Derivatives with respect to the pole give all repeated-pole moments.
    """
    R=H/(H+1)
    apart=S.apart(expr,r)
    A=M=S.S.Zero
    a=S.Symbol('a')
    for term in S.Add.make_args(apart):
        num,den=S.fraction(term)
        if not den.has(r):
            for (n,),c in S.Poly(term,r).terms():
                A+=c*R**(n+1)/(n+1)
                M+=c*R**(n+1)*(S.log(R)/(n+1)-S.Rational(1,(n+1)**2))
            continue
        factors=S.factor_list(den,r)[1]
        poles=[(q,n) for q,n in factors if q.has(r)]
        assert len(poles)==1 and S.degree(poles[0][0],r)==1
        q,n=poles[0];root=S.solve(q,r)[0]
        c=S.cancel(term*(r-root)**n)
        assert not c.has(r)
        A1=S.log(1-R/a)
        M1=S.log(R)*A1+S.polylog(2,R/a)
        A+=c*S.diff(A1,a,n-1).subs(a,root)/S.factorial(n-1)
        M+=c*S.diff(M1,a,n-1).subs(a,root)/S.factorial(n-1)
    return S.expand_func(A),S.expand_func(M)


def integrated_matching():
    data=json.loads(read_artifact(Path(__file__).with_suffix('.json')))
    output={}
    for name,keys in [('rational',['rational']),('logs',['L0','L1'])]:
        f=S.factor(sum(S.sympify(v['coefficient']) for k in keys for v in data[k]))
        assert all(v['power']==-2 for k in keys for v in data[k])
        A,M=integrate_residue(f)
        factor=S.log(H+1) if name=='logs' else 1
        output[name]={'residue_sum':str(f),'pole':str(factor*A),
                      'log_moment':str(factor*M),
                      'finite_shift':str(factor*(A*S.log(H)-M))}
        print('integrated',name,flush=True)
    assert all(v['power']==-1 for v in data['L2'])
    output['LW']={'pole':'0','log_moment':'0','finite_shift':'0'}
    write_artifact(Path(__file__).with_name('compton_real_soft_matching_integrated.json') , json.dumps(output,indent=2)+'\n')


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--integrate',action='store_true',help='Integrate saved exact soft residues')
    args=parser.parse_args()
    if args.integrate: integrated_matching()
    else: main()
