"""Generate complete readable coefficient definitions of our combined primitive."""
from compton_artifacts import read_artifact, write_artifact
import json,sympy as S
from pathlib import Path
B=Path(__file__).resolve().parent;D=B.parent
r,y,u=S.symbols('r y u');data=json.loads(read_artifact(D/'real_equivalence_combined.json'))['original'];co=data['after_weight_three_shuffle']['coefficients'];products=data['log_product_coefficients']
lines=[r'\subsection{Explicit combined rational and \texorpdfstring{$W$}{W} primitive}',r'Put $r=\sqrt{s}>1$. All GPLs below have basepoint $r$ and endpoint $u$. The following finite sum specifies the complete combined primitive, including its rational and elementary terms.',r'\begin{equation}',r' \mathcal F_{R_0+W}(u)=\sum_{j=0}^{37}c_j(u)\,\mathcal G_j(u)+\sum_{j=1}^3 d_j(u)\mathcal L_j(u)',r' +K\left[G(0)\{G(0,s)-G(0,1)\}-G(0)^3/3\right].',r'\end{equation}',r'Here $K=32s(2s+1)/(s-1)^2$. The two weight-three GPL terms already occur in the coefficient table, so the last term contains only their shuffle products. Write $q_\pm=(r^2-3\pm\sqrt{r^4-10r^2+9})/2$ and $a_\pm=1\pm\sqrt{1-r^2}$, $b_\pm=-1\pm\sqrt{1-r^2}$.']
coefficient_start=len(lines)
repl={ (r*r-3-S.sqrt(r**4-10*r*r+9))/2:S.Symbol('q_-'),(r*r-3+S.sqrt(r**4-10*r*r+9))/2:S.Symbol('q_+'),1-S.sqrt(1-r*r):S.Symbol('a_-'),1+S.sqrt(1-r*r):S.Symbol('a_+'),-1-S.sqrt(1-r*r):S.Symbol('b_-'),-1+S.sqrt(1-r*r):S.Symbol('b_+')}
def latex(e):return S.latex(e.subs(y,u))
def word(w):
 if not w:return '1'
 return 'G('+','.join(S.latex(a.xreplace(repl)) for a in w)+';u)'
polydefs={}
def compact_term(term):
 for factor in list(term.atoms(S.Add)):
  if len(latex(factor))>65:
   key=S.Symbol("P_{"+str(len(polydefs)+1)+"}")
   existing=next((k for k,v in polydefs.items() if v==factor),None)
   if existing is not None:key=existing
   else:polydefs[key]=factor
   term=term.xreplace({factor:key})
 return term
def define(label,expr):
 expr=S.factor(S.cancel(S.sympify(expr)).subs(r,S.sqrt(S.Symbol('s'))))
 if len(latex(expr))<180:
  lines.extend([r'\begin{equation}',label+'='+latex(expr)+'.',r'\end{equation}']);return
 # Keep proper fractions instead of expanding a large common numerator.
 parts=S.Add.make_args(S.apart(expr,y))
 lines.extend([r'\begin{align}'])
 for j,term in enumerate(parts):
  st=latex(compact_term(S.factor(term)))
  if j and not st.startswith('-'):st='+'+st
  lines.append((label+'={}&' if j==0 else '&')+st+(r'\\' if j<len(parts)-1 else '.'))
 lines.append(r'\end{align}')
# Small groups keep each display to ordinary page width.
items=list(co.items());assert len(items)==38
for i in range(0,len(items),4):
 lines.extend([r'\begin{align}'])
 for j,(key,val) in enumerate(items[i:i+4],i):
  w=tuple(S.sympify(key));lines.append('\\mathcal G_{'+str(j)+'}&='+word(w)+(r',\\' if j<min(i+4,len(items))-1 else '.'))
 lines.append(r'\end{align}')
 for j,(key,val) in enumerate(items[i:i+4],i):define('c_{'+str(j)+'}',val)
for j,(key,val) in enumerate(products.items(),1):
 a,b=tuple(S.sympify(key));lines.extend([r'\begin{equation}',r'\mathcal L_{'+str(j)+'}='+word((a,))+word((b,))+'.',r'\end{equation}']);define('d_{'+str(j)+'}',val)
polynomial_start=len(lines)
for key,expr in list(polydefs.items()):
 terms=S.Add.make_args(S.expand(expr));chunks=[];part=''
 for term in terms:
  t=latex(term)
  if part and not t.startswith('-'):t='+'+t
  if len(part)+len(t)>60 and part:chunks.append(part);part=t
  else:part+=t
 if part:chunks.append(part)
 lines.extend([r'\begin{align}'])
 for j,c in enumerate(chunks):lines.append((latex(key)+'={}&' if j==0 else '&')+c+(r'\\' if j<len(chunks)-1 else '.'))
 lines.append(r'\end{align}')
lines=lines[:coefficient_start]+['Define the polynomial abbreviations before listing the coefficients.']+lines[polynomial_start:]+lines[coefficient_start:polynomial_start]
lines.append('All coefficients follow by collecting the original two primitives and applying the displayed shuffle identities. The generating script checks their exact coefficient reconstruction before this table is used. Endpoint constants and the common regulator conversion are applied to the combined primitive as before.')
lines += ['', r'\input{real_grouped_endpoint}']
write_artifact(D/'tex/real_combined_coefficients.tex' , '\n'.join(lines)+'\n')
print('wrote complete table',len(lines),'lines')
