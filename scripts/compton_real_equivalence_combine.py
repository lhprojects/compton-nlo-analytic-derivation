"""Combine our R0+W primitives and apply exact weight-two shuffle reduction."""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
from collections import defaultdict
import sympy as S,json
B=Path(__file__).resolve().parent;r,y=S.symbols('r y');G=S.Function('G')
def read(path):
 d=json.loads(read_artifact(path))['coefficients']
 return {(() if k=='constant' else tuple(S.sympify(k,locals={'G':G}).args[0])):S.sympify(c) for k,c in d.items()}
W=read(B/'compton_real_analytic_single_gpl.json')
def simplify(R):
 co=defaultdict(lambda:S.S.Zero)
 for inp in (R,W):
  for w,c in inp.items():co[w]+=c
 co={w:S.factor(c) for w,c in co.items() if S.factor(c)!=0};before=len([w for w in co if w]);products={}
 # G(a,b)+G(b,a)=G(a)G(b); diagonal G(a,a)=G(a)^2/2.
 for w in list(co):
  if w not in co or len(w)!=2:continue
  a,b=w
  if a==b:
   products[w]=co.pop(w)/2
  elif str(a)>str(b) and (b,a) in co:
   rev=(b,a);c=co.pop(w);co[rev]=S.factor(co.get(rev,0)-c);products[rev]=c
 co={w:c for w,c in co.items() if c!=0}
 return co,products,before
out={}
for name,path in [('original',B/'compton_real_analytic_rational_gpl.json'),('qpair',B.parent/'real_equivalence_qpair.json')]:
 co,p,n=simplify(read(path))
 out[name]={'combined_words_before_shuffle':n,'remaining_GPL_words':len([w for w in co if w]),'log_products':len(p),'coefficients':{str(w):str(c) for w,c in co.items()},'log_product_coefficients':{str(w):str(c) for w,c in p.items()}}
 print(name,n,'->',out[name]['remaining_GPL_words'],'GPL +',len(p),'log products',flush=True)
# Verify reconstruction of every coefficient after the shuffle substitution.
for name,path in [('original',B/'compton_real_analytic_rational_gpl.json'),('qpair',B.parent/'real_equivalence_qpair.json')]:
 R=read(path);co,p,n=simplify(R);res=defaultdict(lambda:S.S.Zero,co)
 for (a,b),c in p.items():res[(a,b)]+=c;res[(b,a)]+=c
 for inp in (R,W):
  for w,c in inp.items():res[w]-=c
 assert all(S.cancel(c)==0 for c in res.values())
 out[name]['reconstruction']='every GPL coefficient exactly zero residual'
# The only weight-three terms originate in W. Reduce the equal-coefficient pairs.
orig=out['original']; co={S.sympify(w):S.sympify(c) for w,c in orig['coefficients'].items()}
weight3_products={}
for b in (S.Integer(1),r*r):
 c=co[(0,0,b)];assert S.cancel(c-co[(0,b,0)])==0
 co.pop((0,b,0));co[(0,0,b)]=-c
 weight3_products[str(((0,),(0,b)))]=str(c)
cube=co.pop((0,0,0))/6
orig['after_weight_three_shuffle']={'GPL_words':len([w for w in co if w]),'coefficients':{str(w):str(c) for w,c in co.items()},'log_times_dilog_coefficients':weight3_products,'log_cube_coefficient':str(cube),'certificate':'G(0)G(0,b)=2G(0,0,b)+G(0,b,0), G(0,0,0)=G(0)^3/6; equal coefficients verified exactly for b=1,s.'}
out['decision']='Retain original combined expression: qpair replacement is a branch-structure identity, not a net total-function reduction. Log products are elementary; remaining GPL count is not a proof of functional independence.'
write_artifact(B.parent/'real_equivalence_combined.json' , json.dumps(out,indent=2)+'\n')
