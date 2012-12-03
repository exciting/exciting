#include "array.h"
#include "plugin.h"

class MorseC3ForceConstant: public ArrayFunctionArray<Real> {
  Array<Real> param;
public:
  MorseC3ForceConstant(void): ArrayFunctionArray<Real>(), param() {}
  void eval(Array<Real> *py, const Array<Real> &x);
  int read(istream &file);
};

int MorseC3ForceConstant::read(istream &file) {
  file >> param;
  if (param.get_size()!=10 && param.get_size()!=20) {
    ERRORQUIT("Number of parameters in morsec3 input file must be 10 or 20");
  }
  return 1;
}

void MorseC3ForceConstant::eval(Array<Real> *py, const Array<Real> &x) {
  Real *a=param.get_buf()-1;
  Real r=x(0);
  Real r2=r*r;
  Real r3=r2*r;
  Real c=x(2);
  Real c2=c*c;
  Real c3=c*c*c;
  Real s=(a[1]+a[2]*c+a[3]*c2+a[4]*c3)*(2.0+2.0*(a[5]+a[6]*c+a[7]*c2+a[8]*c3)*a[9]*r+pow(a[5]+a[6]*c+a[7]*c2+a[8]*c3,2.0)*a[9]*a[9]*r2)/(a[5]+a[6]*c+a[7]*c2+a[8]*c3-1.0)/(r3)*exp(-(a[5]+a[6]*c+a[7]*c2+a[8]*c3)*a[9]*(r-a[10]))/2.0+(a[1]+a[2]*c+a[3]*c2+a[4]*c3)*(-2.0*a[5]-2.0*a[6]*c-2.0*a[7]*c2-2.0*a[8]*c3-2.0*(a[5]+a[6]*c+a[7]*c2+a[8]*c3)*a[9]*r-(a[5]+a[6]*c+a[7]*c2+a[8]*c3)*a[9]*a[9]*r2)/(a[5]+a[6]*c+a[7]*c2+a[8]*c3-1.0)/(r3)*exp(-a[9]*(r-a[10]))/2.0;
  if (x.get_size()==20) {
    a=param.get_buf()+10-1;
  }
  Real b=-(a[1]+a[2]*c+a[3]*c*c+a[4]*c*c*c)*(1.0+(a[5]+a[6]*c+a[7]*c*c+a[8]*c*c*c)*a[9]*r)/(a[5]+a[6]*c+a[7]*c*c+a[8]*c*c*c-1.0)/(r*r*r)*exp(-(a[5]+a[6]*c+a[7]*c*c+a[8]*c*c*c)*a[9]*(r-a[10]))/2.0-(a[1]+a[2]*c+a[3]*c*c+a[4]*c*c*c)*(-a[5]-a[6]*c-a[7]*c*c-a[8]*c*c*c-(a[5]+a[6]*c+a[7]*c*c+a[8]*c*c*c)*a[9]*r)/(a[5]+a[6]*c+a[7]*c*c+a[8]*c*c*c-1.0)/(r*r*r)*exp(-a[9]*(r-a[10]))/2.0;
  py->resize(2);
  (*py)(0)=s;
  (*py)(1)=b;
}

SpecificPlugIn<ArrayFunctionArray<Real>,MorseC3ForceConstant> MorseC3ForceConstantRegister("morsec3"); // register plug-in under name "morsec3";

