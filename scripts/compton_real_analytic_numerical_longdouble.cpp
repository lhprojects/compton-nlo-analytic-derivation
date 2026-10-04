// Extended-precision audit of the independent six-diagram numerical integral.
// Independent six-diagram long double-Compton integral with local soft subtraction.
// Units m=1. Returns the coefficient of alpha^3, never a published total.
#include <complex>
#include <cmath>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "cuba.h"
using C=std::complex<long double>;
constexpr long double pi=3.1415926535897932384626433832795L;
struct V {long double a[4];long double &operator[](int i){return a[i];}long double operator[](int i)const{return a[i];}};
V operator+(V a,V b){for(int i=0;i<4;i++)a[i]+=b[i];return a;}
V operator-(V a,V b){for(int i=0;i<4;i++)a[i]-=b[i];return a;}
V operator*(long double s,V a){for(int i=0;i<4;i++)a[i]*=s;return a;}
long double dot(V a,V b){return a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];}
struct M {C a[16]{};C &operator()(int i,int j){return a[4*i+j];}C operator()(int i,int j)const{return a[4*i+j];}};
M operator+(M a,const M&b){for(int i=0;i<16;i++)a.a[i]+=b.a[i];return a;}
M operator*(const M&a,const M&b){M c;for(int i=0;i<4;i++)for(int k=0;k<4;k++)for(int j=0;j<4;j++)c(i,j)+=a(i,k)*b(k,j);return c;}
M operator*(long double x,M a){for(auto &z:a.a)z*=x;return a;}
M slash(V p){M a; a(0,0)=a(1,1)=p[0];a(2,2)=a(3,3)=-p[0];
 a(0,2)=-p[3];a(0,3)=C(-p[1],p[2]);a(1,2)=C(-p[1],-p[2]);a(1,3)=p[3];
 a(2,0)=p[3];a(2,1)=C(p[1],-p[2]);a(3,0)=C(p[1],p[2]);a(3,1)=-p[3];return a;}
M completeness(V p){M a=slash(p);for(int i=0;i<4;i++)a(i,i)+=1;return a;}
M adj(const M&a){M b;for(int i=0;i<4;i++)for(int j=0;j<4;j++)b(i,j)=(i<2?1:-1)*(j<2?1.L:-1.L)*std::conj(a(j,i));return b;}
long double squared(V p,V pp,V k,V k1,V k2,V e1a,V e1b,V e2a,V e2b){
 V ks[3]={k,-1*k1,-1*k2};M Si[3],Sf[3];
 for(int j=0;j<3;j++){Si[j]=(1/(2*dot(p,ks[j])))*completeness(p+ks[j]);Sf[j]=(-1/(2*dot(pp,ks[j])))*completeness(pp-ks[j]);}
 M Es[3][2]={{slash({{0,1,0,0}}),slash({{0,0,1,0}})},{slash(e1a),slash(e1b)},{slash(e2a),slash(e2b)}};
 long double val=0;M Pin=completeness(p),Pout=completeness(pp);
 for(int h0=0;h0<2;h0++)for(int h1=0;h1<2;h1++)for(int h2=0;h2<2;h2++){
  M es[3]={Es[0][h0],Es[1][h1],Es[2][h2]},L[3],R[3],K;
  for(int j=0;j<3;j++){L[j]=es[j]*Sf[j];R[j]=Si[j]*es[j];}
  for(int a=0;a<3;a++)for(int b=0;b<3;b++)if(a!=b){int c=3-a-b;K=K+L[c]*es[b]*R[a];}
  M T=Pout*K*Pin*adj(K);for(int j=0;j<4;j++)val+=T(j,j).real()/4;
 }return val;
}
struct Setup {long double s,E,q,beta,logR;};
void polar(V n,V &a,V &b){
 long double st=std::hypot(n[1],n[2]);long double cp=st>1e-15?n[1]/st:1,sp=st>1e-15?n[2]/st:0;
 a={{0,n[3]*cp,n[3]*sp,-st}};b={{0,-sp,cp,0}};
}
int integrand(const int*,const double x[],const int*,double f[],void *data){
 auto d=*(Setup*)data;long double s=d.s,E=d.E,q=d.q;
 long double ell=std::exp(d.logR*x[0])/(E+q); // E-q=1/(E+q)
 long double c=(ell-E)/q;c=std::clamp(c,-1.L,1.L);long double st=std::sqrt(std::max(0.L,1-c*c));long double jc=ell*d.logR/q;
 V p={{E,0,0,-q}},k={{q,0,0,q}},nb={{1,st,0,c}},pb={{E,-q*st,0,-q*c}};
 V e2a={{0,c,0,-st}},e2b={{0,0,1,0}};
 long double a=(s-1)/2,b=q*ell,dd=1/a-1/b,Born=2*(a/b+b/a+2*dd+dd*dd);
 long double zc=(1-d.beta)*std::exp(d.logR*x[1]);long double cg=(1-zc)/d.beta;
 cg=std::clamp(cg,-1.L,1.L);long double sg=std::sqrt(std::max(0.L,1-cg*cg)),phi=2*pi*x[2];
 long double sum=0,logs=0;
 for(int ch=0;ch<2;ch++){
  V n;
  if(ch==0)n={{1,sg*std::cos(phi),sg*std::sin(phi),-cg}};
  else n={{1,-cg*st+sg*std::cos(phi)*c,sg*std::sin(phi),-cg*c-sg*std::cos(phi)*st}};
  long double ci=-n[3],co=-st*n[1]-c*n[3];
  long double density=.5*d.beta/(2*pi*d.logR)*(1/(1-d.beta*ci)+1/(1-d.beta*co));
  long double c12=n[1]*st+n[3]*c;
  long double ommax=(s-1)/(2*std::sqrt(s)+std::sqrt(4*s-2*(1-c12)*(s-1)));
  long double om=ommax*x[3],D=std::sqrt(s)-om*(1-c12);
  long double w2=(s-1-2*std::sqrt(s)*om)/(2*D);
  V k1=om*n,k2=w2*nb,pp=p+k-k1-k2,e1a,e1b;polar(n,e1a,e1b);
  long double sq=squared(p,pp,k,k1,k2,e1a,e1b,e2a,e2b);
  V j=(1/dot(pb,n))*pb-(1/dot(p,n))*p;
  long double A=-dot(j,j)*Born*q/std::sqrt(s);
  sum+=(ommax*om*w2/D*sq-A/x[3])/density/2;
  logs+=A*std::log(ommax)/density/2;
 }
 f[0]=jc*(sum+logs)/(4*pi*(s-1));
 return 0;
}
int main(int argc,char**argv){
 long double s=argc>1?std::strtold(argv[1],nullptr):1.6,eps=argc>2?std::strtold(argv[2],nullptr):1e-5;int maxeval=argc>3?std::atoi(argv[3]):10000000;
 Setup d{s,(s+1)/(2*std::sqrt(s)),(s-1)/(2*std::sqrt(s)),(s-1)/(s+1),std::log(s)};
 int regions,neval,fail;double result[1],error[1],prob[1];
 Cuhre(4,1,integrand,&d,1,eps,eps*1e-2,0,0,maxeval,0,nullptr,nullptr,&regions,&neval,&fail,result,error,prob);
 std::printf("{\"s_over_m2\":%.17g,\"real_subtracted\":%.17g,\"error\":%.8g,\"neval\":%d,\"fail\":%d}\n",double(s),result[0],error[0],neval,fail);
}
