"""The independently derived CK endpoint using only real logs, Li2 and Li3.

Products of two real logarithms are reduced by the square-of-a-ratio identity.
The Mobius image of infinity is retained explicitly. No quadrature or GPL
 evaluator is called by this module.
"""
from functools import lru_cache
from collections import Counter
import mpmath as mp
from compton_ck_simplification import partials
from compton_soft_analytic import soft_log_letters,born_angular

class ClassicalPrimitive:
 def __init__(self,Z,balanced_products=False):self.Z=Z;self.balanced_products=balanced_products;self.calls=Counter();self.arguments={2:set(),3:set()}
 @lru_cache(None)
 def li(self,n,z):
  self.arguments[n].add(z);return mp.re(mp.polylog(n,z))
 @lru_cache(None)
 def J(self,z,rho):
  # rho=infinity is the constant primitive 0, not an omitted Mobius pole.
  if rho==mp.inf:return mp.mpf(0)
  l=mp.log(z)
  if rho==0:return l**3/3
  return l*l*mp.log(abs(1-z/rho))+2*l*self.li(2,z/rho)-2*self.li(3,z/rho)
 @lru_cache(None)
 def K(self,z,rho):
  l=mp.log(z)
  if rho==0:return l*l/2
  return l*mp.log(abs(1-z/rho))+self.li(2,z/rho)
 @lru_cache(None)
 def simple(self,c,logs):
  if len(logs)==1:
   a=logs[0];return self.K(1-self.Z/a,1-c/a)-self.K(mp.mpf(1),1-c/a)
  a,b=logs
  if a==b:
   return mp.mpf(0) if self.balanced_products else self.J(1-self.Z/a,1-c/a)-self.J(mp.mpf(1),1-c/a)
  ac=1-c/a;bc=1-c/b;zc=ac/bc if bc!=0 else mp.inf
  def endpoint(z):
   A=1-z/a;B=1-z/b
   assert A>0 and B>0
   return ((0 if self.balanced_products else self.J(A,ac)+self.J(B,bc))-self.J(A/B,zc)+self.J(A/B,b/a))/2
  return endpoint(self.Z)-endpoint(mp.mpf(0))
 @lru_cache(None)
 def integral(self,kind,data,logs):
  Z=self.Z;self.calls[(kind,len(logs))]+=1
  if not logs:
   if kind=='p':return Z**(data+1)/(data+1)
   c,n=data
   if n==1:return mp.log(abs((Z-c)/(-c)))
   return ((Z-c)**(1-n)-(-c)**(1-n))/(1-n)
  if kind=='f' and data[1]==1:return self.simple(data[0],logs)
  product=mp.fprod(mp.log(1-Z/a) for a in logs)
  if kind=='p':
   n=data;result=Z**(n+1)/(n+1)*product
   for i,a in enumerate(logs):
    rest=logs[:i]+logs[i+1:]
    result-=sum(a**j*self.integral('p',n-j,rest) for j in range(n+1))/(n+1)
    result-=a**(n+1)*self.integral('f',(a,1),rest)/(n+1)
   return result
  c,n=data;m=n-1;result=-product/(m*(Z-c)**m)
  for i,a in enumerate(logs):
   rest=logs[:i]+logs[i+1:]
   if a==c:result+=self.integral('f',(c,m+1),rest)/m
   else:
    diff=c-a
    result+=sum((-1)**(m-j)/diff**(m-j+1)*self.integral('f',(c,j),rest) for j in range(1,m+1))/m
    result+=(-1)**m/diff**m*self.integral('f',(a,1),rest)/m
  return result

def interior(s,Z,diagnostics=False):
 parts=partials(s);B,L,K=soft_log_letters(mp.mpf(s));engine=ClassicalPrimitive(Z,balanced_products=True);terms=[]
 assert all(sum(c for a,c in letters)==0 for letters in [B,L,K])
 for a,ca in B:
  for kind,data,c in parts[0]:terms.append(ca*c*engine.integral(kind,data,(a,)))
  for index,letters in [(1,L),(2,K)]:
   for b,cb in letters:
    for kind,data,c in parts[index]:terms.append(ca*cb*c*engine.integral(kind,data,tuple(sorted((a,b)))))
 for a,ca in K:
  for kind,data,c in parts[3]:terms.append(ca*c*engine.integral(kind,data,(a,)))
 value=mp.fsum(terms)
 info={'Li2_arguments':len(engine.arguments[2]),'Li3_arguments':len(engine.arguments[3]),'finite_rational_log_terms':len(terms),'primitive_cache':sum(engine.calls.values())}
 return (value,info) if diagnostics else value

def closed(s,dps=40,diagnostics=False):
 # The Taylor coefficients are derived from the original rational/log kernel.
 # This reuses no GPL value, only the existing analytic endpoint expansion.
 from compton_soft_uv_interface import endpoint_series
 with mp.workdps(dps+50):
  ss=mp.mpf(s);x=ss-1;b=x/(x+2);E2=(x+2)**2/(4*ss);degree=8
  delta=mp.power(10,-mp.mpf(dps+20)/(degree+1))/(ss+1)
  U=1-delta;k=mp.sqrt((1+x*U*U)/ss);Z=mp.sqrt(x)*U+mp.sqrt(1+x*U*U)-1
  value,info=interior(ss,Z,True)
  p=2*(x+2)*U*(1-U*U)/x
  boundary=4/(b*E2)*mp.atanh(b*U)/(1-b*b*U*U)*p*mp.atanh(U/k)/k
  aa,cc=endpoint_series(ss,degree+2)
  def tail(n):return sum(delta**(j+1)/(j+1)*(aa[j]*(mp.log(delta)-1/mp.mpf(j+1))+cc[j]) for j in range(n+1))
  value+=boundary+tail(degree)-mp.atanh(b)/b*born_angular(ss)
  change=abs(tail(degree+2)-tail(degree));assert change<mp.mpf(10)**(-dps-5)*max(1,abs(value))
  info['tail_change_estimate']=mp.nstr(change,10)
  return (+value,info) if diagnostics else +value
