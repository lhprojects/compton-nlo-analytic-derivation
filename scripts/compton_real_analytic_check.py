"""Cross-check exact invariant trace against transverse Dirac matrices."""
from compton_artifacts import read_artifact
import json,itertools
from pathlib import Path
import numpy as np
import sympy as sp
from check_compton_amplitudes import final,real,pol,spin,dot
D=json.loads(read_artifact(Path(__file__).with_name('compton_real_analytic_trace.json')))
H,x,y,U,V=sp.symbols('H x y U V')
expr=sum(t['multiplicity']*sp.sympify(t['numerator'])/sp.sympify(t['denominator']) for t in D['terms'])
f=sp.lambdify((H,x,y,U,V),expr,'numpy')
for w in (.15,.8,5,40):
 p,k,k1,k2,pp=final(w,.04*w);P=p+k
 analytic=f(2*dot(p,k),2*dot(P,k1),2*dot(P,k2),2*dot(p,k1),2*dot(p,k2))
 direct=sum(spin(real(p,k,k1,k2,*es),p,pp) for es in itertools.product(pol(k),pol(k1),pol(k2)))
 print(w,analytic,direct,analytic/direct-1)
 assert abs(analytic/direct-1)<1e-10
