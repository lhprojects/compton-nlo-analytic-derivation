"""Exact rational endpoint audit of the four-sector common primitive.
The algebraic coefficient of each power of b=s-u is extracted without
approximating any GPL. This identifies the required endpoint corrections.
"""
from compton_artifacts import read_artifact, write_artifact
import json
from pathlib import Path
import sympy as S
folder=Path(__file__).resolve().parents[1]
s,y,ell=S.symbols('s y ell')
data=json.loads(read_artifact(folder/'real_equivalence_allfour.json'))
result={}
for group in ['coefficients','log_products']:
    result[group]={}
    for word,c in data[group].items():
        c=S.cancel(S.sympify(c));n,d=S.fraction(c)
        order=0; dd=d
        while S.expand(dd.subs(y,s))==0:
            dd=S.diff(dd,y);order+=1
        # Removing the exact factor avoids a general-purpose series expansion.
        regular=S.cancel(c*(s-y)**order)
        terms={}
        for k in range(order+1):
            value=S.factor((-1)**k*S.diff(regular,y,k).subs(y,s)/S.factorial(k))
            if value!=0:terms[str(k-order)]=str(value)
        result[group][word]={'pole_order':order,'series_through_constant':terms}
    print(group,'maximum pole',max(v['pole_order'] for v in result[group].values()),flush=True)
result['variable']='b=s-u > 0'
result['scope']='Laurent coefficients only; singular coefficients must multiply the corresponding GPL Taylor terms before taking the finite part'
write_artifact(folder/'real_allfour_endpoint.json' , json.dumps(result,indent=2)+'\n')
# The only singular function coefficient multiplies G(0;u).
assert [w for w,v in result['coefficients'].items() if v['pole_order'] and w!='()']==['(0,)']
assert all(v['pole_order']==0 for v in result['log_products'].values())
pole0=S.sympify(result['coefficients']['()']['series_through_constant']['-1'])
poleg=S.sympify(result['coefficients']['(0,)']['series_through_constant']['-1'])
finite_taylor=S.factor(-poleg/s)
D=S.factor(pole0+poleg*ell/2)
def endpoint(group,word):
    return S.sympify(result[group][word]['series_through_constant'].get('0','0'))
k=S.Symbol('log_2sqrt_s_over_splus1')
assert S.cancel(endpoint('coefficients','(s, -s)')+endpoint('coefficients','(s, -1)'))==0
A=S.collect(S.factor(-endpoint('coefficients','(s,)')-endpoint('log_products','(0, s)')*ell/2-endpoint('coefficients','(s, -s)')*k),[ell,k])
r=S.Symbol('r',positive=True)
base=S.factor(S.sympify(data['coefficients']['()']).subs(s,r*r).subs(y,r))
base_s=-2*(s+1)*(s**4-3*s**3-41*s**2+27*s-8)/(3*(s-1))
assert S.cancel(base-base_s.subs(s,r*r))==0
# All GPL terms vanish at the lower endpoint. Products of two weight-one
# functions vanish quadratically, while their coefficients have at most a
# simple pole there.
for group in ['coefficients','log_products']:
    for word,c in data[group].items():
        if not S.sympify(word):continue
        e=S.cancel(S.sympify(c).subs(s,r*r));dd=S.denom(e);order=0
        while S.expand(dd.subs(y,r))==0:
            dd=S.diff(dd,y);order+=1
        assert order <= (1 if group=='log_products' else 0),(word,order)
result['finite_upper_Taylor_correction']=str(finite_taylor)
result['lower_endpoint']=str(base_s)
result['double_pole_integrand_coefficient']=str(D)
result['simple_pole_integrand_coefficient']=str(A)
result['log_2sqrt_s_over_splus1']='log(2*sqrt(s)/(s+1))'
result['zcut_conversion']='FP_b F - F(sqrt(s)) + A log((s-1)/s) - D/(s*(s-1))'
result['lower_endpoint_certificate']='all nonconstant primitive terms vanish at u=sqrt(s), including products with simple rational poles'
write_artifact(folder/'real_allfour_endpoint.json' , json.dumps(result,indent=2)+'\n')
print('finite Taylor correction',finite_taylor,flush=True)
print('lower endpoint',base_s,flush=True)
print('D',D,flush=True)
