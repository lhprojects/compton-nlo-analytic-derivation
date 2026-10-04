"""Collect all four real finite endpoints before adding UW and soft K.
No reference result is read. The expression follows the original endpoint
coefficients and the two root-pair identities derived in this supplement.
"""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S
from compton_real_qpair_endpoint import r,Li,Log,loga,l2,g,canonical,converted
folder=Path(__file__).resolve().parents[1]
d=json.loads(read_artifact(folder/'real_allfour_endpoint.json'))
s,ell=S.symbols('s ell');l=Log(r)
def cv(e):return S.sympify(e).subs({s:r*r,ell:2*l})
def const(row):return cv(row['series_through_constant'].get('0','0'))
# Real endpoint logs, including the constant term of G(s;u).
logs={S.S.Zero:l,S.S.One:loga(r+1),r*r:-loga(r*r-r),-r*r:loga(2*r/(r+1)),-S.S.One:loga((r*r+1)/(r+1))}
expr=const(d['coefficients']['()'])+cv(d['finite_upper_Taylor_correction'])
qplus=(s-3+S.sqrt(s*s-10*s+9))/2;qminus=(s-3-S.sqrt(s*s-10*s+9))/2
asym=[1+S.sqrt(1-s),1-S.sqrt(1-s)]
C=128*r*r*(r*r+5)/(3*(r*r-1));ca=None
weights={S.S.Zero:-1,S.S.One:1,s:-1,-s:2,-S.S.One:-2}
removed_q=removed_a=removed_three=0
K=32*r*r*(2*r*r+1)/(r*r-1)**2
three={(0,0,0):-2,(0,0,1):-1,(0,0,s):1,(0,1,0):-1,(0,s,0):1}
for key,row in d['coefficients'].items():
    w=S.sympify(key)
    if not w:continue
    c=const(row)
    if len(w)==3:
        assert S.cancel(c-K*three[w])==0;removed_three+=1;continue
    if len(w)==2 and any(S.simplify(w[0]-v)==0 for v in [qplus,qminus]):
        assert S.cancel(c-C*weights[w[1]])==0;removed_q+=1;continue
    if len(w)==2 and any(S.simplify(w[0]-v)==0 for v in asym):
        cc=c/({0:1,1:1,s:-1}[w[1]])
        if ca is None:ca=cc
        assert S.cancel(ca-cc)==0;removed_a+=1;continue
    w=tuple(a.subs(s,r*r) for a in w)
    expr+=c*(logs[w[0]] if len(w)==1 else g(*w))
for key,row in d['log_products'].items():
    a,b=S.sympify(key);expr+=const(row)*logs[a.subs(s,r*r)]*logs[b.subs(s,r*r)]
assert (removed_q,removed_a,removed_three)==(10,6,5)
Li3=S.Function('Li3')
expr+=C*converted+ca*(S.pi**2/4+l2(-1/r)-l2(1/r)+l*loga(r*r-1))
expr+=K*(S.zeta(3)-Li3(1/(r*r))-l*l2(1/(r*r))-S.pi**2*l/6-2*l**3/3)
expr-=cv(d['lower_endpoint'])
expr+=cv(d['simple_pole_integrand_coefficient']).subs(S.Symbol('log_2sqrt_s_over_splus1'),loga(2*r/(r*r+1)))*loga((r*r-1)/(r*r))
expr-=cv(d['double_pole_integrand_coefficient'])/(r*r*(r*r-1))
# Transform all remaining real dilogarithms to the same arguments as the q pair.
expr=S.expand(expr.xreplace({v:canonical(v.args[0]) for v in expr.atoms(Li)}))
# Treat logarithms and special functions as independent atoms for coefficient collection.
atoms=sorted(expr.atoms(Li,Li3,Log)|{S.log(2),S.zeta(3),S.pi},key=str)
syms=S.symbols('v:'+str(len(atoms)));forward=dict(zip(atoms,syms));back=dict(zip(syms,atoms))
poly=S.Poly(expr.xreplace(forward),*syms);terms=[]
for powers,c in poly.terms():
    c=S.factor(S.cancel(c))
    if c!=0:terms.append({'powers':powers,'coefficient':str(c)})
out={'source':'combined four real endpoint, lower boundary and z-cut conversion',
     'variable':'r=sqrt(s)>1','normalization':'unnormalized F_R0^z+F_W^z+F_log^z+F_D^z',
     'atoms':list(map(str,atoms)),'terms':terms,'number_of_products':len(terms),
     'scope':'z-cut finite part; matching shifts, UV, UW and soft K not included',
     'removed_GPL_groups':{'q_pair':removed_q,'a_pair':removed_a,'weight_three':removed_three}}
write_artifact(folder/'real_four_finite_collected.json' , json.dumps(out,indent=2)+'\n')
print('Finite collection',len(terms),'products',len(atoms),'atoms',flush=True)
