"""Emit the exact angular coefficient list into the owned supplement source."""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
import json
import sympy as S
H,x,y,z=S.symbols('H x y z')
root=Path(__file__).resolve().parents[1]
p=root/'tex/real_analytic.tex'
marker='% BEGIN EXACT REAL ANGULAR COEFFICIENTS'
full=read_artifact(p)
s=full.split(marker)[0]
end_marker='% END EXACT REAL ANGULAR COEFFICIENTS'
suffix=full.split(end_marker,1)[1] if end_marker in full else ''
lines=[marker,r'\subsubsection{Coefficients of the integrated direction moments}',
       'The following coefficients use $y=H-x+z$ to shorten the polynomials.',
       'Every exponent triple not displayed has zero coefficient.']
for row in json.loads(read_artifact(Path(__file__).with_name('compton_real_analytic_angles.json'))):
    e=row['powers'];name='C_{'+','.join(str(a) for a in e)+'}'
    c=S.factor(S.sympify(row['coefficient']).subs(y,H-x+z));n,d=S.fraction(c)
    terms=S.Add.make_args(S.expand(n))
    chunks=[];chunk=''
    for term in terms:
        tx=S.latex(term);tx=(' '+tx if tx.startswith('-') else ' + '+tx)
        if len(chunk)+len(tx)>105 and chunk:chunks.append(chunk);chunk=tx
        else:chunk+=tx
    if chunk:chunks.append(chunk)
    if len(chunks)==1:
        lines.extend([r'\begin{equation}',name+'='+S.latex(c)+'.',r'\end{equation}'])
    else:
        pname='P_{'+','.join(str(a) for a in e)+'}'
        lines.extend([r'\begin{equation}',name+'=\\frac{'+pname+'}{'+S.latex(d)+'}.',r'\end{equation}',r'\begin{equation}',r'\begin{aligned}'])
        for i,ch in enumerate(chunks):
            if i==0:ch=ch.removeprefix(' + ')
            lines.append((pname+'={}&' if i==0 else '&')+ch+(r'\\' if i<len(chunks)-1 else '.'))
        lines.extend([r'\end{aligned}',r'\end{equation}'])
write_artifact(p , s+'\n'.join(lines)+'\n'+end_marker+suffix)
