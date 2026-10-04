"""Shorten the complete four-sector finite part using real five-term identities."""
from compton_artifacts import read_artifact, write_artifact
from pathlib import Path
import json
import sympy as S
from compton_real_qpair_endpoint import r,Li,Log,loga,l2
folder=Path(__file__).resolve().parents[1]
d=json.loads(read_artifact(folder/'real_four_finite_collected.json'))
atoms=[S.sympify(a,locals={'r':r,'L2':Li,'logabs':Log}) for a in d['atoms']]
expr=sum(S.sympify(t['coefficient'],locals={'r':r})*S.prod(a**p for a,p in zip(atoms,t['powers'])) for t in d['terms'])
# Use canonical rational arguments without expanding inside special functions.
expr=expr.xreplace({v:Li(S.factor(v.args[0])) for v in expr.atoms(Li)})
def half_identity(z):
    return Li(1/z)+Li(1/(z+1))+S.pi**2/12-S.log(2)**2/2+S.log(2)*loga((z+1)/z)-Li(2/(z+1))
a=(r+1)/(r*r+1);b=(r-1)/(r*(r+1));u=1/r;v=1/(r*r+1);w=1/(r+1)
f=lambda z:loga(z)*loga(1-z)
E=S.expand((f(u)+f(v)-f(w)-f(a)+f(b))/2)
replacements={Li(S.factor((r+1)/(2*r))):half_identity(r),
              Li(S.factor((r*r+1)/(2*r*r))):half_identity(r*r),
              Li(S.factor(a)):Li(S.factor(b))+Li(u)+Li(v)-Li(w)+E}
short=S.expand(expr.xreplace(replacements))
# Differentiate identities from d Li2(z)=-log(1-z) dz/z.
def differential(e):
    atoms=sorted(e.atoms(Li,Log),key=str);vs=S.symbols('q:'+str(len(atoms)))
    p=e.xreplace(dict(zip(atoms,vs)));out=S.diff(p,r)
    for atom,var in zip(atoms,vs):
        z=atom.args[0]
        deriv=-loga(1-z)*S.diff(z,r)/z if atom.func==Li else S.diff(z,r)/z
        out+=S.diff(p,var)*deriv
    return S.factor(S.expand(out.xreplace(dict(zip(vs,atoms)))))
identities=[Li((r+1)/(2*r))-half_identity(r),Li(a)-Li(b)-Li(u)-Li(v)+Li(w)-E]
for identity in identities:assert differential(identity)==0
newatoms=sorted(short.atoms(Li,S.Function('Li3'),Log)|{S.log(2),S.zeta(3),S.pi},key=str)
vs=S.symbols('a:'+str(len(newatoms)));poly=S.Poly(short.xreplace(dict(zip(newatoms,vs))),*vs)
terms=[]
for powers,c in poly.terms():
    c=S.factor(S.cancel(c))
    if c!=0:terms.append({'powers':powers,'coefficient':str(c)})
used={i for t in terms for i,p in enumerate(t['powers']) if p}
indices=sorted(used)
newatoms=[newatoms[i] for i in indices]
for t in terms:t['powers']=[t['powers'][i] for i in indices]
result={**d,'atoms':list(map(str,newatoms)),'terms':terms,'number_of_products':len(terms),
        'reduction':'two real five-term identities, with exact zero derivative certificates and zero constants at r->infinity',
        'dilogarithm_count':sum(a.func==Li for a in newatoms),
        'trilogarithm_count':sum(a.func==S.Function('Li3') for a in newatoms)}
write_artifact(folder/'real_four_finite_short.json' , json.dumps(result,indent=2)+'\n')
print(result['dilogarithm_count'],'dilogarithms',result['trilogarithm_count'],'trilogarithms',len(terms),'products')
# Duplication plus the negative-argument Landen identity removes both
# square-root dilogarithm arguments in their common-coefficient difference.
root_identity=Li(1/r)-Li(1/(r*r))/2-(Log(r+1)-Log(r))**2/2
root_free=S.expand(short.xreplace({Li(1/(r+1)):root_identity}))
# All occurrences of Li2(1/r) must cancel, rather than being relabelled.
assert S.cancel(root_free.coeff(Li(1/r)))==0
s,ls,lh,lp=S.symbols('s Ls Lh Lp')
root_free=S.expand(root_free.subs({Log(r-1):lh-Log(r+1),Log(r):ls/2,Log(r*r+1):lp}))
assert S.cancel(root_free.coeff(Log(r+1),1))==0
assert S.cancel(root_free.coeff(Log(r+1),2))==0
root_free=S.expand(root_free)
ats=sorted(root_free.atoms(Li,S.Function('Li3'))|{ls,lh,lp,S.log(2),S.zeta(3),S.pi},key=str)
vs=S.symbols('c:'+str(len(ats)));p=S.Poly(root_free.xreplace(dict(zip(ats,vs))),*vs)
terms2=[]
for powers,c in p.terms():
    c=S.cancel(c);n,den=S.fraction(c)
    def even(poly):
        terms=S.Poly(poly,r).terms();assert all(k[0]%2==0 for k,v in terms)
        return sum(v*s**(k[0]//2) for k,v in terms)
    c=S.factor(S.cancel(even(n)/even(den)))
    if c!=0:terms2.append({'powers':powers,'coefficient':str(c)})
indices=sorted({i for t in terms2 for i,power in enumerate(t["powers"]) if power})
ats=[ats[i].subs(r*r,s) for i in indices]
for t in terms2:t["powers"]=[t["powers"][i] for i in indices]
result2={**result,'variable':'s>1, Ls=log(s), Lh=log(s-1), Lp=log(s+1)',
         'atoms':list(map(str,ats)),'terms':terms2,'number_of_products':len(terms2),
         'dilogarithm_count':sum(a.func==Li for a in ats),'square_root_cancellation':'Exact zero coefficient of Li2(1/r) and all log(r+1) powers after Landen and duplication',
         'scope':'complete four-sector z-cut finite part; matching shifts, UV, UW and soft K not included'}
write_artifact(folder/'real_four_finite_rootfree.json' , json.dumps(result2,indent=2)+'\n')
print('Root-free:',result2['dilogarithm_count'],'dilogarithms',len(terms2),'products')
