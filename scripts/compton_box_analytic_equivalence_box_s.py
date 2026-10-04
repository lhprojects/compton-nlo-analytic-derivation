"""Exact constant-letter reduction of the direct spanning-box channel.

All coefficients come from the independent box trace.  Moving endpoint
GPLs are converted recursively by their total differential and their
q=1 boundary, not by comparing values with a published cross section.
"""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
from collections import defaultdict
from functools import lru_cache
import json,itertools,time
import sympy as S
from compton_virtual_angular_analytic import primitive,g,G,y
from compton_gpl_algebra import multiply,kernel_map

ROOT=Path(__file__).resolve().parents[1]
b,q,x,v=S.symbols('b q x v');ll,lx,ds=S.symbols('ll lx ds')
B=(1-q)/(1+q);V=(1-q)**2/(1+q*q)


def clean(d):return {w:S.factor(c) for w,c in d.items() if S.cancel(c)!=0}


def add_to(out,data,scale=1):
    for w,c in data.items():out[w]+=scale*c


@lru_cache(None)
def dlog_letters(f):
    f=S.factor(f)
    if f==0:return {}
    n,d=S.fraction(f);out=defaultdict(lambda:S.S.Zero)
    for poly,sign in [(n,1),(d,-1)]:
        for factor,m in S.factor_list(poly,q)[1]:
            roots=S.roots(factor,q)
            assert sum(roots.values())==S.degree(factor,q)
            for root,mult in roots.items():out[S.simplify(root)]+=sign*m*mult
    result=clean(out)
    assert S.cancel(S.diff(f,q)/f-sum(c/(q-a) for a,c in result.items()))==0
    return result


@lru_cache(None)
def endpoint_word_general(word,endpoint):
    """GPL(word;endpoint) pulled to a q=1 tangential basepoint.

    Letters and endpoint are rational functions of q.  The nonzero-letter
    endpoint/a limit must vanish, which is checked rather than assumed.
    """
    if not word:return {():S.S.One}
    if all(a==0 for a in word):
        residues=dlog_letters(endpoint);order=residues.get(S.Integer(1),0)
        leading=S.limit(endpoint/(1-q)**order,q,1)
        assert leading!=0 and leading not in (S.oo,-S.oo,S.zoo)
        logE={(a,):c for a,c in residues.items()}
        if leading!=1:logE[()]=S.log(leading)
        out={():S.S.One}
        for _ in word:out=multiply(out,logE)
        return clean({w:c/S.factorial(len(word)) for w,c in out.items()})
    for a in word:
        if a!=0:assert S.limit(endpoint/a,q,1)==0,(word,endpoint,a)
    pieces=[(word[1:],endpoint-word[0],1)]
    if word[-1]!=0:pieces.append((word[:-1],-word[-1],-1))
    for j in range(len(word)-1):
        if word[j]==word[j+1]:continue
        difference=word[j]-word[j+1]
        pieces.append((word[:j+1]+word[j+2:],difference,1))
        pieces.append((word[:j]+word[j+1:],difference,-1))
    out=defaultdict(lambda:S.S.Zero)
    for lower,f,sign in pieces:
        for a,c in dlog_letters(f).items():
            for w,d in endpoint_word_general(lower,endpoint).items():out[(a,)+w]+=sign*c*d
    out=clean(out)
    assert not any(w and all(a==1 for a in w) for w in out),(word,out)
    return out


@lru_cache(None)
def endpoint_word(word):
    letters=tuple(S.factor(S.sympify(a).subs(b,B)) for a in word)
    return endpoint_word_general(letters,V)


def angular_primitive(rows=None,label='s'):
    source=ROOT/'scripts/compton_box_analytic_angle.json'
    if rows is None:rows=json.loads(read_artifact(source))['s']['integrand']
    total=0
    for j,row in enumerate(rows):
        c=S.sympify(row['coefficient']).subs(x,4*b*b/(1-b*b)).subs(v,y)
        c=S.factor(c)
        word=[]
        for value in row['word']:
            a=S.sympify(value)
            if a.has(x):
                assert S.cancel(a*a-x/(x+4))==0,a
                a=-b if a.could_extract_minus_sign() else b
            word.append(a)
        word=tuple(word)
        t=time.time();total+=primitive(c,word)
        print('integrated',label,'box word',word,'seconds',time.time()-t,flush=True)
    atoms=total.atoms(G);out=defaultdict(lambda:S.S.Zero)
    for atom in atoms:
        word=tuple(atom.args[0]);out[word]+=S.factor(total.coeff(atom))
    out[()]+=S.factor(total.subs({a:0 for a in atoms}))
    return clean(out)


def build():
    primitive_data=angular_primitive()
    path=ROOT/'compton_box_analytic_equivalence_box_s_primitive.json'
    write_artifact(path , json.dumps([dict(word=list(map(str,w)),coefficient=str(c)) for w,c in primitive_data.items()],indent=2)+'\n')
    out=defaultdict(lambda:S.S.Zero)
    upper=2*b*b/(1+b*b)
    for w,c in primitive_data.items():
        value=S.factor(c.subs(y,upper).subs(b,B))
        add_to(out,endpoint_word(w),value)
        if not w:out[()]-=S.factor(c.subs(y,0).subs(b,B))
    out=clean(out)
    angular_out=out
    X=(1-q)**2/q;maps=kernel_map();ds_words=defaultdict(lambda:S.S.Zero)
    for a,c in maps[-1].items():
        for d,e in maps[0].items():ds_words[(S.simplify(a),S.simplify(d))]+=c*e/X
    scalar_words={ds:dict(ds_words),lx:{(S.Integer(1),):2,(S.Integer(0),):-1}}
    full=defaultdict(lambda:S.S.Zero)
    for word,c in angular_out.items():
        for powers,coef in S.Poly(S.expand(c),ds,lx).terms():
            term={word:coef}
            for symbol,n in zip((ds,lx),powers):
                for _ in range(n):term=multiply(term,scalar_words[symbol])
            add_to(full,term)
    out=clean(full)
    assert all(not c.has(ds,lx,b,x) for c in out.values())
    assert max(map(len,out))<=3
    report=dict(source='Independent direct spanning box',normalization='J_s, unnormalized angular integral',q='(sqrt(x+4)-sqrt(x))/(sqrt(x+4)+sqrt(x))',coefficients=[dict(word=list(map(str,w)),coefficient=str(c)) for w,c in out.items()],C_s_coefficients=[dict(word=list(map(str,w)),coefficient=str(S.factor(c/(2*X*X)))) for w,c in out.items()],status='Complete direct-box reduction, including scalar factors. Only log(lambda/m) remains as ll.')
    write_artifact(ROOT/'compton_box_analytic_equivalence_box_s.json' , json.dumps(report,indent=2)+'\n')
    print('s box constant-letter words:',len(out),flush=True)


if __name__=='__main__':build()
