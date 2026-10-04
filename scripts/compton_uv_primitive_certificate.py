"""Verify each UV angular primitive before summing large rational expressions."""
from compton_artifacts import read_artifact, write_artifact
import json,time
from pathlib import Path
from collections import defaultdict
import sympy as S
from compton_soft_uv_constant_basis import q,H,v,LH,X,multiply,primitive,G,y,ROOT
folder=ROOT

def check(component):
    rows=json.loads(read_artifact(ROOT/'scripts/compton_real_uv_analytic.json'))[component]
    logv={(-S.S.One,):S.S.One}
    logk={(q-1,):S.S.One,(1/q-1,):S.S.One,(-S.S.One,):-S.S.One}
    kernels=defaultdict(lambda:S.S.Zero)
    for row in rows:
        co=S.sympify(row['coefficient']).subs(S.log(H),LH).subs(H,X).subs(v,1+y)
        words={():S.S.One}
        for power,base in zip(row['powers'],[logv,logk]):
            for _ in range(power):words=multiply(words,base)
        for word,c in words.items():kernels[word]+=co*c
    completed=[]
    for word,co in kernels.items():
        start=time.time()
        co=S.factor(co)
        f=primitive(co,word)
        atoms=f.atoms(G)
        terms={tuple(a.args[0]):S.factor(f.coeff(a)) for a in atoms}
        terms[()]=S.factor(f.xreplace({a:0 for a in atoms}))
        residual=defaultdict(lambda:S.S.Zero)
        for w,c in terms.items():
            residual[w]=S.cancel(residual[w]+S.diff(c,y))
            if w:residual[w[1:]]=S.cancel(residual[w[1:]]+c/(y-w[0]))
        residual[word]-=co
        for w,c in residual.items():
            assert S.cancel(S.together(c))==0,(component,word,w)
        completed.append({'word':list(map(str,word)),'residual_coefficients':len(residual),'all_zero':True,'seconds':time.time()-start})
        print(component,word,'exact zero',round(time.time()-start,2),flush=True)
        write_artifact(folder/('uv_primitive_certificate_'+component+'.json') , json.dumps({'component':component,'complete':len(completed)==len(kernels),'method':'differentiate each independent primitive and cancel every rational GPL coefficient separately','checks':completed},indent=2)+'\n')

if __name__=='__main__':
    check('pole')
    check('finite')
