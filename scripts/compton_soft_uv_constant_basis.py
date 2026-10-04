"""Independent UV-sector reduction to constant-letter GPLs based at q=1.
Only the independently integrated UV Laurent coefficients enter.  No
reference total or difference against a reference fixes any coefficient.
"""
from compton_artifacts import read_artifact, write_artifact
import json,time
from pathlib import Path
from functools import lru_cache
from collections import defaultdict
import sympy as S
from compton_virtual_angular_analytic import primitive,G,y
from compton_box_analytic_equivalence_box_s import dlog_letters,clean,add_to,endpoint_word_general
from compton_gpl_algebra import multiply
ROOT=Path(__file__).resolve().parents[1]
q,H,v,LH=S.symbols('q H v LH');X=(1-q)**2/q

@lru_cache(None)
def pullback(word):
 if not word:return {():S.S.One}
 if all(a==0 for a in word):
  logX={(S.Integer(1),):2,(S.Integer(0),):-1};out={():S.S.One}
  for _ in word:out=multiply(out,logX)
  return {w:c/S.factorial(len(word)) for w,c in out.items()}
 pieces=[(word[1:],X-word[0],1)]
 if word[-1]!=0:pieces.append((word[:-1],-word[-1],-1))
 for j in range(len(word)-1):
  if word[j]==word[j+1]:continue
  f=word[j]-word[j+1]
  pieces += [(word[:j+1]+word[j+2:],f,1),(word[:j]+word[j+1:],f,-1)]
 out=defaultdict(lambda:S.S.Zero)
 for lower,f,sign in pieces:
  for a,c in dlog_letters(f).items():
   for w,d in pullback(lower).items():out[(a,)+w]+=sign*c*d
 out=clean(out)
 assert not any(w and all(a==1 for a in w) for w in out),(word,out)
 return out

@lru_cache(None)
def origin_series(word,N):
 if not word:return [S.S.One]+[S.S.Zero]*N
 assert word[-1]!=0,word
 tail=origin_series(word[1:],N);a=word[0];out=[S.S.Zero]*(N+1)
 for n in range(1,N+1):
  out[n]=tail[n]/n if a==0 else -sum(tail[j]/a**(n-j) for j in range(n))/n
 return out

def build(component):
 rows=json.loads(read_artifact(ROOT/'scripts/compton_real_uv_analytic.json'))[component]
 logv={(-S.S.One,):S.S.One}
 logk={(q-1,):S.S.One,(1/q-1,):S.S.One,(-S.S.One,):-S.S.One}
 kernels=defaultdict(lambda:S.S.Zero)
 for row in rows:
  co=S.sympify(row['coefficient']).subs(S.log(H),LH).subs(H,X).subs(v,1+y)
  wds={():S.S.One}
  for power,base in zip(row['powers'],[logv,logk]):
   for _ in range(power):wds=multiply(wds,base)
  for word,c in wds.items():kernels[word]+=co*c
 total=S.S.Zero
 for word,co in kernels.items():
  start=time.time();total+=primitive(S.factor(co),word)
  print('UV',component,'primitive',word,'seconds',time.time()-start,flush=True)
 atoms=total.atoms(G);terms={tuple(a.args[0]):S.factor(total.coeff(a)) for a in atoms}
 terms[()]=S.factor(total.xreplace({a:0 for a in atoms}))
 terms=clean(terms)
 derivative=defaultdict(lambda:S.S.Zero)
 for word,co in terms.items():
  derivative[word]+=S.diff(co,y)
  if word:derivative[word[1:]]+=co/(y-word[0])
 for word,co in kernels.items():derivative[word]-=co
 for word,co in derivative.items():assert S.cancel(co)==0,('primitive derivative',word)
 print('all primitive derivative coefficients zero',flush=True)
 for word in terms:
  a=pullback(word);b=endpoint_word_general(word,X)
  assert all(S.cancel(a.get(w,0)-b.get(w,0))==0 for w in set(a)|set(b))
 print('independent general pullback agrees',flush=True)
 moving=dict(component=component,kernels=[dict(word=list(map(str,w)),coefficient=str(c)) for w,c in kernels.items()],primitive=[dict(word=list(map(str,w)),coefficient=str(c)) for w,c in terms.items()])
 write_artifact(ROOT/('soft_uv_moving_basis_'+component+'.json') , json.dumps(moving,indent=2)+'\n')
 # Normalize the combined lower endpoint by its exact finite Taylor coefficient.
 lower=0
 for word,co in terms.items():
  if word and all(a==0 for a in word):
   assert S.cancel(co)==0,('uncancelled logarithm',word,co)
  if not word:
   lower+=sum(S.residue(term/y,y,0) for term in S.Add.make_args(S.apart(co,y)))
  else:
   num,den=S.fraction(co);dp=S.Poly(den,y)
   m=min(k[0] for k,c in dp.terms())
   series=origin_series(word,m)
   for n,c in enumerate(series):
    if c:lower+=c*sum(S.residue(term*y**(n-1),y,0) for term in S.Add.make_args(S.apart(co,y)))
 lower=S.factor(lower)
 print('lower finite constant',lower,flush=True)
 out=defaultdict(lambda:S.S.Zero)
 for word,co in terms.items():add_to(out,pullback(word),S.factor(co.subs(y,X)))
 out[()]-=lower
 # The coefficient log(H) uses the same left-tangential log as pure-zero words.
 result=defaultdict(lambda:S.S.Zero)
 for word,co in out.items():
  for (power,),c in S.Poly(S.expand(co),LH).terms():
   piece={word:c}
   for _ in range(power):piece=multiply(piece,{(S.Integer(1),):2,(S.Integer(0),):-1})
   add_to(result,piece)
 result=clean(result)
 report=dict(component=component,source='Independent UV inner Laurent coefficients',basepoint='q=1 left tangential',H=str(X),normalization='before 1/(8*s*H)',lower_constant=str(lower),coefficients=[dict(word=list(map(str,w)),coefficient=str(c)) for w,c in result.items()])
 write_artifact(ROOT/('soft_uv_constant_basis_'+component+'.json') , json.dumps(report,indent=2)+'\n')
 print('constant-letter terms',len(result),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--component',default='pole',choices=['pole','finite']);a=p.parse_args();build(a.component)
