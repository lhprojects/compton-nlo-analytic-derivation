"""Collect products and remove redundant GPL words using exact shuffles."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
from collections import defaultdict
import sympy as S
from compton_total_energy_collect import shuffle

folder=Path(__file__).resolve().parents[1]
X,Y,L,M=S.symbols('X Y L M')
rows=json.loads(read_artifact(folder/'total_energy_collected.json'))['terms']
c={tuple(r['word']):S.sympify(r['coefficient']) for r in rows}
out=defaultdict(lambda:S.S.Zero)
for w,a in c.items():
    if w==('-1',):out[()]+=a*L
    elif w==('0',):out[()]+=a*M
    elif w==('-1','-1'):out[()]+=a*L*L/2
    elif w==('-1','-1','-1'):out[()]+=a*L**3/6
    elif w==('-1','0'):
        out[()]+=a*L*M;out[('0','-1')]-=a
    elif w==('-1','-1','0'):
        out[()]+=a*M*L*L/2
        out[('0','-1')]-=a*L
        out[('0','-1','-1')]+=a
    elif w==('-1','0','-1'):
        out[('0','-1')]+=a*L
        out[('0','-1','-1')]-=2*a
    elif w==('0','r','-1'):
        out[('r','-1')]+=a*M
        out[('r','0','-1')]-=a
        out[('r','-1','0')]-=a
    else:out[w]+=a
out={w:S.factor(a) for w,a in out.items() if S.cancel(a)!=0}
# Re-expand every log product into words to certify the whole replacement.
back=defaultdict(lambda:S.S.Zero)
for w,a in out.items():
    for (l,m),v in S.Poly(S.expand(a),L,M).terms():
        words={w:S.S.One}
        for letter,n in [('-1',l),('0',m)]:
            for _ in range(n):
                result=defaultdict(lambda:S.S.Zero)
                for word,co in words.items():
                    for q in shuffle((letter,),word):result[q]+=co
                words=result
        for word,co in words.items():back[word]+=co*v
assert all(S.cancel(back[w]-c.get(w,0))==0 for w in set(back)|set(c))
report={'definition':'C1=sum coefficient(X,Y,L,M)*G(word;X), L=ln(1+X), M=ln(X)',
        'terms':[{'word':list(w),'coefficient':str(v)} for w,v in out.items()],
        'function_count':sum(bool(w) for w in out),
        'certificate':'complete shuffle reconstruction equals all original 23 energy-word coefficients exactly'}
write_artifact(folder/'total_energy_compact.json' , json.dumps(report,indent=2)+'\n')
print('nontrivial GPL functions:',report['function_count'])
for w,a in out.items():print(w,a)
