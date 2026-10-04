"""Collect all independently integrated sectors in one energy GPL basis."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
from collections import defaultdict
from functools import lru_cache
import sympy as S

folder=Path(__file__).resolve().parents[1]
X,Y,s,ll=S.symbols('X Y s ll');z2=S.pi**2/6
@lru_cache(None)
def shuffle(a,b):
    if not a:return (b,)
    if not b:return (a,)
    return tuple((a[0],)+w for w in shuffle(a[1:],b))+tuple((b[0],)+w for w in shuffle(a,b[1:]))
def add(out,data,scale=1):
    for w,c in data.items():out[w]+=scale*c
def multiply(a,b):
    out=defaultdict(lambda:S.S.Zero)
    for w,c in a.items():
        for v,d in b.items():
            for word in shuffle(w,v):out[word]+=c*d
    return dict(out)
def power(a,n):
    out={():S.S.One}
    for _ in range(n):out=multiply(out,a)
    return out
def read(name,key):
    return {tuple(r['word']):S.sympify(r['coefficient']).subs({S.Symbol('x'):X,ll:0})
            for r in json.loads(read_artifact(folder/name))[key]}

if __name__=='__main__':
    parts={}
    parts['virtual_uv']=read('total_simplification_virtual_uv.json','coefficients')
    parts['uw']={w:c/(4*X*(1+X)) for w,c in read('uw_inverse_basis.json','coefficients_X').items()}
    parts['soft']={w:-c/(1+X) for w,c in read('ck_energy_collected.json','terms').items()}
    L={('-1',):1};M={('0',):1};J={('-2',):1}
    D={():z2,('-1','0'):1,('-1','-1'):-1}
    V={():z2,('-2','0'):1,('-2','-2'):-1}
    T={('-1',):z2,('-1','-1','0'):1,('-1','-1','-1'):-1}
    c={k:S.sympify(v).subs(s,1+X) for k,v in json.loads(read_artifact(folder/'real_four_matched_coefficients.json')).items()}
    real=defaultdict(lambda:S.S.Zero)
    for i in range(4):add(real,power(L,i),c[f'A_{i}'])
    add(real,{():1},6*z2*c['B_0']);add(real,L,6*z2*c['B_1'])
    add(real,D,c['C_0']);add(real,multiply(L,D),c['C_1'])
    add(real,T,c['K']);add(real,V,c['N']);add(real,power(J,2),c['E']);add(real,J,c['F'])
    for i in range(3):add(real,multiply(M,power(L,i)),c[f'H_{i}'])
    add(real,multiply(M,J),c['H_3'])
    parts['real_four']={w:a/(8*X*(1+X)) for w,a in real.items()}
    P=X**5-22*X**4-184*X**3-444*X**2-424*X-144
    Q=X**4+2*X**3-24*X**2-80*X-48
    R=5*X**3+60*X**2+100*X+48
    ir={('-1','-1'):-4*(X**2-4*X-6)/X**3,
        ('-1',):-2*P/(X**3*(1+X)**2*(X+4)),
        ('r','-1'):12*Y*Q/(X**4*(X+4)),():-R/(X**2*(1+X)**2)}
    parts['regulator_conversion']=multiply(ir,{('0',):-1,('-1',):-S.Rational(1,2)})
    out=defaultdict(lambda:S.S.Zero)
    for data in parts.values():add(out,data)
    total={w:S.factor(S.cancel(a)) for w,a in out.items() if S.cancel(a)!=0}
    vanished=[list(w) for w,a in out.items() if S.cancel(a)==0]
    report={'definition':'complete C1(X) multiplying alpha_T^3/m^2; X>0, Y=sqrt(X/(X+4))',
            'source':'own virtual+UV, UW, four matched real sectors, finite soft energy equations, and infrared conversion',
            'parts':{n:[{'word':list(w),'coefficient':str(S.factor(a))} for w,a in data.items()] for n,data in parts.items()},
            'terms':[{'word':list(w),'coefficient':str(a)} for w,a in sorted(total.items(),key=lambda p:(len(p[0]),p[0]))],
            'vanishing_words':vanished,'word_count':len(total),
            'numerical_check':'pending'}
    write_artifact(folder/'total_energy_collected.json' , json.dumps(report,indent=2)+'\n')
    print('Complete total:',len(total),'words,',len(vanished),'exact cancellations',flush=True)
    for w,a in total.items():print(w,a,flush=True)
