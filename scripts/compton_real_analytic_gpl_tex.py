"""Write explicit rational coefficients of the L_W GPL primitive."""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
import sympy as S,json
root=Path(__file__).resolve().parents[1];p=root/'tex/real_analytic.tex'
mark='% BEGIN REAL W GPL COEFFICIENTS';end='% END REAL W GPL COEFFICIENTS';whole=read_artifact(p);start=whole.split(mark)[0];suffix=whole.split(end,1)[1] if end in whole else ''
d=json.loads(read_artifact(Path(__file__).with_name('compton_real_analytic_single_gpl.json')))['coefficients'];r,y,u=S.symbols('r y u')
items=[('00','G((0, 0), y)'),('01','G((0, 1), y)'),('10','G((1, 0), y)'),('s0','G((r**2, 0), y)'),('0','G((0,), y)'),('1','G((1,), y)'),('\\mathrm c','constant')]
lines=[mark]
for label,key in items:
 c=S.factor(S.sympify(d[key]).subs(y,u));n,den=S.fraction(c);name='A_{'+label+'}'
 chunks=[];chunk=''
 for t in S.Add.make_args(S.expand(n)):
  tx=S.latex(t);tx=(' '+tx if tx.startswith('-') else ' + '+tx)
  if chunk and len(chunk)+len(tx)>105:chunks.append(chunk);chunk=tx
  else:chunk+=tx
 if chunk:chunks.append(chunk)
 if len(chunks)==1:lines.extend([r'\begin{equation}',name+'='+S.latex(c)+'.',r'\end{equation}'])
 else:
  pname='N_{'+label+'}'
  lines.extend([r'\begin{equation}',name+'=\\frac{'+pname+'}{'+S.latex(den)+'}.',r'\end{equation}',r'\begin{equation}',r'\begin{aligned}'])
  for i,ch in enumerate(chunks):
   if i==0:ch=ch.removeprefix(' + ')
   lines.append((pname+'={}&' if i==0 else '&')+ch+(r'\\' if i<len(chunks)-1 else '.'))
  lines.extend([r'\end{aligned}',r'\end{equation}'])
write_artifact(p , start+'\n'.join(lines)+'\n'+end+suffix)
