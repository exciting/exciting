#include "array.h"
#include "plugin.h"

class MorseForceConstant: public ArrayFunctionArray<Real> {
  Array<Real> param;
public:
  MorseForceConstant(void): ArrayFunctionArray<Real>(), param() {}
  void eval(Array<Real> *py, const Array<Real> &x);
  int read(istream &file);
};

int MorseForceConstant::read(istream &file) {
  file >> param;
  if (param.get_size()!=8 && param.get_size()!=16) {
    ERRORQUIT("Number of parameters in morse input file must be 8 or 16");
  }
  return 1;
}

void MorseForceConstant::eval(Array<Real> *py, const Array<Real> &x) {
  Real *a=param.get_buf()-1;
  Real r=x(0);
  Real r2=r*r;
  Real r3=r2*r;
  Real c=x(2);
  Real c2=c*c;
  Real s=0.5*(a[1]+a[2]*c+a[3]*c2)*(2.+2.*(a[4]+a[5]*c+a[6]*c2)*a[7]*r+ipow(a[4]+a[5]*c+a[6]*c2,2)*ipow(a[7],2)*r2)/(a[4]+a[5]*c+a[6]*c2-1)/r3*exp(-(a[4]+a[5]*c+a[6]*c2)*a[7]*(r-a[8]))+0.5*(a[1]+a[2]*c+a[3]*c2)*(-2*a[4]-2*a[5]*c-2*a[6]*c2-2*(a[4]+a[5]*c+a[6]*c2)*a[7]*r-(a[4]+a[5]*c+a[6]*c2)*ipow(a[7],2)*r2)/(a[4]+a[5]*c+a[6]*c2-1)/r3*exp(-a[7]*(r-a[8]));
  if (x.get_size()==16) {
    a=param.get_buf()+8-1;
  }
  Real b=-0.5*(a[1]+a[2]*c+a[3]*c2)*(1+(a[4]+a[5]*c+a[6]*c2)*a[7]*r)/(a[4]+a[5]*c+a[6]*c2-1)/r3*exp(-(a[4]+a[5]*c+a[6]*c2)*a[7]*(r-a[8]))-0.5*(a[1]+a[2]*c+a[3]*c2)*(-a[4]-a[5]*c-a[6]*c2-(a[4]+a[5]*c+a[6]*c2)*a[7]*r)/(a[4]+a[5]*c+a[6]*c2-1)/r3*exp(-a[7]*(r-a[8]));
  py->resize(2);
  (*py)(0)=s;
  (*py)(1)=b;
}

SpecificPlugIn<ArrayFunctionArray<Real>,MorseForceConstant> MorseForceConstantRegister("morse"); // register plug-in under name "morse";
