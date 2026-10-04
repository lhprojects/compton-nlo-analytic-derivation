"""Exact extra audit of the sector pole sum, no numerical quadrature."""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S
H,v,t=S.symbols('H v t');root=Path(__file__).parent
u=json.loads(read_artifact(root/'compton_real_uv_analytic.json'))['pole']
w=json.loads(read_artifact(root/'compton_real_analytic_UW_endpoints.json'))
R={tuple(a['powers']):S.sympify(a['coefficient']) for a in u}
uwL=S.sympify(w['LH1']).subs(t,(v-1)/H)/H
uwR=(-S.sympify(w['LH2'])*t/(1+H*t)-S.sympify(w['TH1'])).subs(t,(v-1)/H)/H
T=(v-1)**2/v;y=H-T
B=2*(H/y+y/H-4*(1/y-1/H)+4*(1/y-1/H)**2)
# physical target 8 H integral dc B (b-1), oriented v=1..s.
targetL=16*(H+1)/H*B*(v*v+1)/v**2
targetR=-16*(H+1)/H*B*(1-1/v**2)
L=S.factor(R[(1,0)]+2*uwL-targetL)
Q=S.factor(R[(0,0)]+2*uwR-targetR)
report=root.parent/'soft_pole_certificate.json'
# These rational candidates were obtained by Hermite reduction.  The proof
# checks their derivatives against the raw sector inputs, not their provenance.
d={k:S.sympify(val) for k,val in json.loads(read_artifact(report))['formulae'].items()}
Ar,Jr,C=d['Ar'],d['Jr'],d['C']
quad=v*v-(H+2)*v+1
cq=-2*(H+1)*(H**3-8*H**2-8*H+48)/H
c1=d['c1'];cp=d['cp']
assert S.cancel(L-S.diff(Ar,v)-C/v)==0
assert S.cancel(Q-Ar/v-S.diff(Jr,v)-c1/v-cp/(v+1)-cq*S.diff(quad,v)/quad)==0
print('Both rational derivative identities: exact zero',flush=True)
# Endpoint constants: all apparent v=1 singularities cancel in the combined
# primitive.  Ar has order two, so the Taylor series of log(v) fixes its
# finite part from precisely these two Laurent coefficients.
ar2=S.cancel((v-1)**2*Ar)
a2=S.factor(ar2.subs(v,1));a1=S.factor(S.diff(ar2,v).subs(v,1))
Jconst=sum(S.residue(term/(v-1),v,1) for term in S.Add.make_args(S.apart(Jr,v)))
low=S.factor(a1-a2/2+Jconst)
assert S.cancel(low-d['lower_rational'])==0
assert S.cancel(Jr.subs(v,H+1)-d['upper_rational'])==0
assert S.cancel(Ar.subs(v,H+1)-d['upper_log_coeff'])==0
# Persist the formerly one-off origin cancellation checks for UV and UW.
z=S.symbols('z')
lnseries=sum((-1)**(n+1)*z**n/S.Integer(n) for n in range(1,8))
for name,kernel in [('UV',R[(1,0)].subs(v,z+1)*lnseries+R[(0,0)].subs(v,z+1)),
                    ('UW',uwL.subs(v,z+1)*lnseries+uwR.subs(v,z+1))]:
 num,den=S.fraction(S.cancel(kernel));pn=S.Poly(num,z);pd=S.Poly(den,z)
 n0=min(k[0] for k,c in pn.terms() if c!=0)
 d0=min(k[0] for k,c in pd.terms() if c!=0)
 assert n0-d0>=1,(name,n0,d0)
 print(name,'combined endpoint starts at power',n0-d0,flush=True)
ls,l2=S.symbols('Ls L2')
E=d['upper_rational']-low+(d['upper_log_coeff']+c1)*ls+C*ls**2/2+cp*l2
raw=json.loads(read_artifact(root/'compton_real_soft_matching_integrated.json'))
bose=json.loads(read_artifact(root/'compton_real_bose_subtraction.json'))
other=sum(S.sympify(raw[k]['pole']) for k in ['rational','logs'])+S.sympify(bose['pole'])
changes={}
for log in other.atoms(S.log):
 arg=S.factor(log.args[0])
 if S.cancel(arg-(H+1))==0:changes[log]=ls
 elif S.cancel(arg-1/(H+1))==0:changes[log]=-ls
 elif S.cancel(arg-(H+2)**2/(4*(H+1)))==0:changes[log]=2*l2-ls
 else:raise ValueError('unexpected logarithm '+str(arg))
poly=S.Poly(S.expand(E+other.xreplace(changes)),ls,l2)
certificate={str(powers):str(S.factor(co)) for powers,co in poly.terms()}
assert all(S.cancel(co)==0 for _,co in poly.terms())
print('All endpoint logarithm polynomial coefficients:',certificate,flush=True)
output=json.loads(read_artifact(report))
output['certificate']=certificate
output['checks']=['both rational derivative identities zero','finite basepoint constants checked',
 'UV and UW pole densities vanish at lower endpoint','all final log polynomial coefficients zero']
write_artifact(report , json.dumps(output,indent=2)+'\n')
