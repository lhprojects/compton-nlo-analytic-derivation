"""Independent exact traces for integrated Compton vertex/self-energy terms.

Uses the Clifford trace recursion from the independent real-amplitude trace.
All scalar form factors remain symbols; coefficients are exact rational
functions of x=s-1 and y=1-u, in electron-mass units.
"""
from compton_artifacts import write_artifact
from collections import defaultdict
from pathlib import Path
import json
import sympy as S
from compton_real_analytic_trace import trace,canon


def add(*terms):
    out=defaultdict(lambda:S.S.Zero)
    for a in terms:
        for w,c in a.items(): out[w]+=c
    return {w:c for w,c in out.items() if c!=0}


def scale(a,f):return {w:c*f for w,c in a.items()}

def mul(*terms):
    out={():S.S.One}
    for a in terms:
        nxt=defaultdict(lambda:S.S.Zero)
        for w,c in out.items():
            for v,d in a.items():nxt[w+v]+=c*d
        out=dict(nxt)
    return out


def slash(i):return {(i,):S.S.One}
I={():S.S.One}
def num(i):return add(slash(i),I)
def contract(a):return {(-3,)+w+(-3,):c for w,c in a.items()}


def vertex(left,right,eps,on,off,J,cz):
    j00,j10,j11,j20,j21,j22=J
    P1=num(right);P2=num(left);E=slash(eps)
    H=slash(on);Z=add(slash(off),scale(H,-1))
    N0=contract(mul(P2,E,P1))
    U0=scale(contract(add(mul(H,E,P1),mul(P2,E,H))),-1)
    U1=scale(contract(add(mul(Z,E,P1),mul(P2,E,Z))),-1)
    V0=contract(mul(H,E,H))
    V1=contract(add(mul(H,E,Z),mul(Z,E,H)))
    V2=contract(mul(Z,E,Z))
    D=x if off==2 else -y
    return add(scale(N0,-j00),scale(U0,-j10),scale(U1,-j11),
               scale(V0,-j20),scale(V1,-j21),scale(V2,-j22),
               scale(E,-S.Rational(3,2)+j20+D*j21+cz))


x,y=S.symbols('x y',positive=True)
cz=S.Symbol('cz',real=True)
gram=S.Matrix([[1,x/2,y/2],[x/2,0,(x-y)/2],[y/2,(x-y)/2,0]])
vec=[S.Matrix(v) for v in [(1,0,0),(1,1,-1),(1,1,0),(1,0,-1)]]
dots={(i,j):S.expand((vec[i].T*gram*vec[j])[0]) for i in range(4) for j in range(i,4)}


def traced(a):
    out=0
    for w,c in a.items():
        factor=sum(n*S.prod(dots[p] for p in mon) for mon,n in trace(canon(w)))
        out+=c*factor/4
    return S.expand(out)


def main():
    born=add(scale(mul(slash(-2),num(2),slash(-1)),1/x),
             scale(mul(slash(-1),num(3),slash(-2)),-1/y))
    conjugate={tuple(reversed(w)):c for w,c in born.items()}
    results={}
    for channel,Q,A,B,D in [('s',2,-1,-2,x),('u',3,-2,-1,-y)]:
        J=S.symbols(' '.join(channel+'J'+k for k in ['00','10','11','20','21','22']))
        a,b=S.symbols(channel+'a '+channel+'b')
        prop=scale(num(Q),1/D)
        vleft=vertex(1,Q,B,1,Q,J,cz)
        vright=vertex(Q,0,A,0,Q,J,cz)
        selfenergy=add(scale(slash(Q),a-cz),scale(I,b-4+cz))
        kernels={'vertex_left':mul(vleft,prop,slash(A)),
                 'vertex_right':mul(slash(B),prop,vright),
                 'self_energy':mul(slash(B),prop,selfenergy,prop,slash(A))}
        for name,kernel in kernels.items():
            result=traced(mul(num(1),kernel,num(0),conjugate))
            variables=list(J)+[a,b,cz]
            coeffs={str(v):str(S.factor(result.coeff(v))) for v in variables if result.coeff(v)!=0}
            constant=result.subs({v:0 for v in variables})
            if constant!=0:coeffs['constant']=str(S.factor(constant))
            results[channel+'_'+name]=coeffs
            print(channel,name,len(str(coeffs)),flush=True)
    write_artifact(Path(__file__).with_suffix('.json') , json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
