// Scalar integrals only. No Compton amplitude or cross-section formula.
#include <qcdloop/qcdloop.h>
extern "C" void compton_box(double S,double t,double lambda,double *out) {
 static ql::QCDLoop<ql::qcomplex,ql::qdouble,ql::qdouble> integral;
 std::vector<ql::qcomplex> r(3);
 std::vector<ql::qdouble> m={lambda < 0 ? 0 : (__float128)lambda*lambda,1,1,1};
 std::vector<ql::qdouble> p={1,0,0,1,S,t};
 integral.integral(r,1,m,p);
 for(int i=0;i<3;++i){out[2*i]=(double)crealq(r[i]);out[2*i+1]=(double)cimagq(r[i]);}
}
