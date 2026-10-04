"""Verify the new CK GPL primitive against the old closed primitive."""
from compton_artifacts import write_artifact
from pathlib import Path
import json
import mpmath as mp
from compton_ck_simplification import expression,certificate
from compton_soft_analytic_gpl import expression as old_expression,evaluate_gpl_taylor

def boundary(s,Z):
 x=s-1;b=x/(x+2);E2=(x+2)**2/(4*s);w=1+Z
 U=(w-1/w)/(2*mp.sqrt(x));k=mp.sqrt((1+x*U*U)/s)
 p=2*(x+2)*U*(1-U*U)/x
 return 4/(b*E2)*mp.atanh(b*U)/(1-b*b*U*U)*p*mp.atanh(U/k)/k

def check():
 certificate();rows=[]
 with mp.workdps(65):
  for s in map(mp.mpf,['1.6','10']):
   Z=(mp.sqrt(s)+mp.sqrt(s-1)-1)*mp.mpf('.7')
   new=expression(s,Z);old=old_expression(s,Z)
   nv,nwords,_=evaluate_gpl_taylor(new,Z);ov,owords,_=evaluate_gpl_taylor(old,Z)
   error=abs(nv+boundary(s,Z)-ov)
   assert error<mp.mpf('1e-50')
   # No imaginary GPL letter survives, including those from the IBP recurrence.
   letters={a for w in new for a in w};assert all(mp.im(a)==0 for a in letters)
   row=dict(s=str(s),fraction='0.7',new_terms=len(new),old_terms=len(old),new_family=nwords,old_family=owords,absolute_difference=mp.nstr(error,12),all_new_letters_real=True)
   rows.append(row);print(row,flush=True)
 folder=Path(__file__).resolve().parents[1]
 write_artifact(folder/'ck_simplification_check.json' , json.dumps(rows,indent=2)+'\n')
if __name__=='__main__':check()
