"""Print the exact independently merged coefficients in reviewable TeX."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1];F=ROOT
X,Y,ll=S.symbols('X Y ll')
def tex(e):return S.latex(e).replace('ll',r'\ell_\lambda')
rows=json.loads(read_artifact(F/'total_simplification_virtual_uv.json'))['coefficients']
out=[r'\paragraph{Explicit partial sum.}',r'Write $G_{a_1,\ldots,a_n}=G(a_1,\ldots,a_n;X)$ and $G_{\varnothing}=1$,',r'with the kernels above, $\omega_{-2}=\dd X/(X+2)$ and $\omega_1=\dd X/(X-1)$.',r'These iterated integrals have base point zero:',r'\begin{equation}',r' G_{a,w}(X)=\int_0^X\omega_a(t)G_w(t),\qquad G_{\underbrace{\scriptstyle 0,\ldots,0}_{n}}(X)=\frac{\ln^nX}{n!}.',r'\end{equation}',r'The words ending in zero that contain the letter $1$ are understood by their continuous real extension at $X=1$.',r'Let $\ell_\lambda=\ln(\lambda_{\mathrm{phys}}/m)$.  The following expression includes the complete virtual contribution and the finite UV real contribution:',r'\begin{equation}',r' C_{V+\mathrm{UV}}(X;\ell_\lambda)=\sum_{w\in\mathcal W}c_w(X;\ell_\lambda)G_w(X).',r'\end{equation}',r'The set $\mathcal W$ and every coefficient are specified by the following equations.']
for row in rows:
 label=','.join(row['word']) or r'\varnothing';co=S.sympify(row['coefficient']);num,den=S.fraction(co)
 terms=S.Add.make_args(S.expand(num));terms=sorted(terms,key=lambda t:(-t.as_powers_dict().get(X,0),str(t)))
 chunks=[S.Add(*terms[i:i+3]) for i in range(0,len(terms),3)]
 lhs=tex(den)+r'\,c_{'+label+'}' if den!=1 else r'c_{'+label+'}'
 out.extend([r'\begin{equation}',r'\begin{aligned}'])
 for i,c in enumerate(chunks):
  lead=lhs+'&={}' if i==0 else '&\quad{}+'
  out.append(lead+'\\bigl('+tex(c)+r'\bigr)'+(r'\\' if i<len(chunks)-1 else '.'))
 out.extend([r'\end{aligned}',r'\end{equation}'])
write_artifact(F/'tex/total_simplification_coefficients.tex' , '\n'.join(out)+'\n')
