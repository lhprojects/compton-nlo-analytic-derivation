"""Add independently derived finite matching shifts to the shortened real sum."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S
from compton_dilog_identities import certify, ordinary as certified_ordinary, quadratic as certified_quadratic
folder=Path(__file__).resolve().parents[1]
s,x,Ls,Lh,Lp,c=S.symbols('s x Ls Lh Lp c');Li=S.Function('L2');Li3=S.Function('Li3')
d=json.loads(read_artifact(folder/'real_four_finite_rootfree.json'))
atoms=[S.sympify(a) for a in d['atoms']]
expr=sum(S.sympify(t['coefficient'])*S.prod(a**p for a,p in zip(atoms,t['powers'])) for t in d['terms'])
Dinv=Li(1/s);Dplus=Li(1/(s+1));Dhalf=Li(2/(s+1))
ordinary=Dinv-Ls*Lh+Ls*Ls/2-S.pi**2/6
E=S.expand(2*(Lh-Lp)*(c-Lp)-2*(Lh-c-Ls)*(Lp-c-Ls)+(c+Ls-Lp)**2-(2*Lp-2*c-Ls)**2/2-2*c*(Lp-Ls)+c*c-S.pi**2/6)
quadratic=4*Dhalf-2*Dinv-2*Dplus+E.subs(c,S.log(2))
# Equality on the connected real domain follows from the exact derivative
# identity and the common analytic boundary value, not sampled agreement.
checks=certify()
assert S.expand(ordinary-certified_ordinary)==0
assert S.expand(quadratic-certified_quadratic.subs(c,S.log(2)))==0
m=json.loads(read_artifact(folder/'matching_simplified.json'))
basis=[Ls,Lh,Lp-S.log(2),ordinary,quadratic]
for t in m['terms']:
    expr+=S.sympify(t['coefficient']).subs(x,s-1)*S.prod(b**p for b,p in zip(basis,t['powers']))
expr=S.expand(expr)
ats=sorted(expr.atoms(Li,Li3)|{Ls,Lh,Lp,S.log(2),S.pi,S.zeta(3)},key=str)
vs=S.symbols('v:'+str(len(ats)));poly=S.Poly(expr.xreplace(dict(zip(ats,vs))),*vs)
terms=[]
for powers,co in poly.terms():
    co=S.factor(S.cancel(co))
    if co!=0:terms.append({'powers':powers,'coefficient':str(co)})
used=sorted({i for t in terms for i,p in enumerate(t['powers']) if p})
ats=[ats[i] for i in used]
for t in terms:t['powers']=[t['powers'][i] for i in used]
out={**d,'normalization':'F_R0^z+F_W^z+F_log^z+F_D^z+Delta_R0+Delta_log+Delta_D',
     'scope':'matched four-sector real finite part; UV, 2 UW and soft K are not included',
     'atoms':list(map(str,ats)),'terms':terms,'number_of_products':len(terms),
     'dilogarithm_count':sum(a.func==Li for a in ats),'quadratic_dilog_identity':str(quadratic),
     'dilogarithm_identity_certificates':checks}
write_artifact(folder/'real_four_matched_short.json' , json.dumps(out,indent=2)+'\n')
print('Matched four-sector sum:',len(terms),'products,',out['dilogarithm_count'],'dilogarithms')
for t in terms:print('*'.join(str(a)+'^'+str(p) for a,p in zip(ats,t['powers']) if p) or '1',t['coefficient'])
