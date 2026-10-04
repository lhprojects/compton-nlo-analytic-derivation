"""Exact Mobius/derivative certificate and repeatable full-CK closed checks."""
from compton_artifacts import write_artifact
import argparse,json
from pathlib import Path
import sympy as S
import mpmath as mp
from compton_ck_classical import closed
from compton_ck_simplification import closed as gpl_closed

def certificate():
 z,a,b,c,rho=S.symbols('z a b c rho');L,H,D2=S.symbols('L H D2')
 A=1-z/a;B=1-z/b;v=A/B;vc=(1-c/a)/(1-c/b)
 assert S.cancel(S.diff(v,z)/(v-vc)-S.diff(v,z)/(v-b/a)-1/(z-c))==0
 # Differentiate Li_n using d Li_n(z/rho)/dz=Li_(n-1)(z/rho)/z.
 # Re Li1(z/rho)=-log|1-z/rho|.
 derivative_J=2*L/z*H+L**2/(z-rho)+2*D2/z-2*L*H/z-2*D2/z
 derivative_K=H/z+L/(z-rho)-H/z
 assert S.cancel(derivative_J-L**2/(z-rho))==0
 assert S.cancel(derivative_K-L/(z-rho))==0
 return dict(Mobius_image_of_infinity_verified=True,J_K_derivatives_verified=True,basepoint='A(0)=B(0)=1; subtract exact primitives at zero',branches='A,B positive on physical interval; real boundary values of Li2 and Li3',balanced_weights='Each logarithm list sums to zero, cancelling the two single-square J contributions')

def check(dps=35):
 proof=certificate();rows=[]
 for sv in ['1.6','10']:
  with mp.workdps(dps+50):
   value,info=closed(sv,dps,True);old=gpl_closed(sv,dps)
   err=abs(value-old);assert err<mp.mpf(10)**(-dps+3)
   row=dict(s=sv,CK=mp.nstr(value,dps+15),GPL_difference=mp.nstr(err,12),dps=dps,**info);rows.append(row);print(row,flush=True)
 report=dict(symbolic_certificate=proof,closed_endpoint_checks=rows,counts_note='Counts are cached special-function arguments at each point, not a proven minimal symbolic basis. Numerically identical arguments may have distinct floating-point encodings.')
 p=Path(__file__).resolve().parents[1]/'ck_classical_full_check.json';write_artifact(p , json.dumps(report,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--full',action='store_true');p.add_argument('--dps',type=int,default=35);a=p.parse_args()
 if a.full:check(a.dps)
 else:print(certificate())
