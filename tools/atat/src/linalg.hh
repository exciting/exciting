#ifndef _LINALG_H_
#define _LINALG_H_

#include <complex.h>
#include "linsolve.h"
#include "arraylist.h"

void transpose_matrix(Array2d<Real> *trans,const Array2d<Real> &a);
void diff(Array<Real> *result, const Array<Real> &a, const Array<Real> &b);
void diff(Array2d<Real> *result, const Array2d<Real> &a, const Array2d<Real> &b);
void sum(Array<Real> *result, const Array<Real> &a, const Array<Real> &b);
void sum(Array2d<Real> *result, const Array2d<Real> &a, const Array2d<Real> &b);
void sum(Array2d<Complex> *result, const Array2d<Complex> &a, const Array2d<Complex> &b);
Real mean(const Array<Real> &a);
void product(Array2d<Real> *result, const Array2d<Real> &a, const Array2d<Real> &b);
void product_diag(Array2d<Real> *result, const Array<Real> &diag, const Array2d<Real> &b);
void product_diag(Array2d<Real> *result, const Array2d<Real> &b, const Array<Real> &diag);
void product_diag(Array<Real> *result, const Array<Real> &diag, const Array<Real> &b);
void product(Array<Real> *result, const Array2d<Real> &a, const Array<Real> &b);
void product(Array<Real> *result, const Array<Real> &a, const Array2d<Real> &b);
void product(Array2d<Real> *result, const Array2d<Real> &a, Real b);
void product(Array<Real> *result, const Array<Real> &a, Real b);
Real inner_product(const Array<Real> &a, const Array<Real> &b);
void inner_product(Array2d<Real> *result, const Array2d<Real> &a, const Array2d<Real> &b);
void inner_product(Array<Real> *result, const Array2d<Real> &a, const Array<Real> &b);
Real normalize(Array<Real> *a);
void outer_product(Array2d<Real> *result, const Array2d<Real> &a, const Array2d<Real> &b);
void outer_product(Array2d<Real> *result, const Array<Real> &a, const Array<Real> &b);
void outer_product_diag(Array<Real> *result, const Array2d<Real> &a, const Array2d<Real> &b);
Real quadratic_form(const Array2d<Real> &a, const Array<Real> &x);
Real quadratic_form(const Array2d<Complex> &a, const Array<Complex> &x);
Real trace(const Array2d<Real> &a);
void set_identity(Array2d<Real> *result, int n);
void to_real(Array<Real> *p_r, const Array<int> &i);
Real max_norm(const Array<Real> &v);
void invert_matrix_tall(Array2d<Real> *inv, const Array2d<Real> &a);
void invert_matrix_wide(Array2d<Real> *inv, const Array2d<Real> &a);

template<class T>
class LinearInterpolator {
  Array<T> a;
  Real x0;
  Real x1;
 public:
  LinearInterpolator(void): a() {
    x0=0.; x1=0.;
  }
  LinearInterpolator(Real _x0, Real _x1, const Array<T> &_a) {
    init(_x0,_x1,_a);
  }
  LinearInterpolator(const LinearInterpolator<T> &li) {
    init(li.x0,li.x1,li.a);
  }

  void init(Real _x0, Real _x1, const Array<T> &_a) {
#ifdef DEBUG
    if (_a.get_size()<2) ERRORQUIT("array size<2 in LinearInterpolator");
#endif
    a=_a;
    x0=_x0;
    x1=_x1;
  }
  T operator ()(Real x) {
    int n=a.get_size();
    Real y=(Real)(n-1)*(x-x0)/(x1-x0);
    int iy=(int)floor(y);
    if (iy>n-2) iy=n-2;
    if (iy<0) iy=0;
    Real fy=y-(Real)iy;
    return (a(iy)*(1.-fy)+a(iy+1)*fy);
  }
};

template<class T>
class LinearInterpolatorBig {
  Array<T> a;
  Real x0;
  Real x1;
 public:
  LinearInterpolatorBig(void): a() {
    x0=0.; x1=0.;
  }
  LinearInterpolatorBig(Real _x0, Real _x1, const Array<T> &_a) {
    init(_x0,_x1,_a);
  }
  LinearInterpolatorBig(const LinearInterpolatorBig<T> &li) {
    init(li.x0,li.x1,li.a);
  }
  void init(Real _x0, Real _x1, const Array<T> &_a) {
#ifdef DEBUG
    if (_a.get_size()<2) ERRORQUIT("array size<2 in LinearInterpolator");
#endif
    a=_a;
    x0=_x0;
    x1=_x1;
  }
  void interpol(T *py, Real x) const {
    int n=a.get_size();
    Real y=(Real)(n-1)*(x-x0)/(x1-x0);
    int iy=(int)floor(y);
    if (iy>n-2) iy=n-2;
    if (iy<0) iy=0;
    Real fy=y-(Real)iy;
    T tmp1,tmp2;
    product(&tmp1,a(iy),1.-fy);
    product(&tmp2,a(iy+1),fy);
    sum(py,tmp1,tmp2);
  }
};

template<class T>
class PolyInterpolatorBig {
  Array<T> a;
  Real x0;
  Real x1;
 public:
  PolyInterpolatorBig(void): a() {
    x0=0.; x1=0.;
  }
  PolyInterpolatorBig(Real _x0, Real _x1, const Array<T> &_a): a() {
    init(_x0,_x1,_a);
  }
  PolyInterpolatorBig(const PolyInterpolatorBig<T> &li): a(li.a) {
    x0=li.x0; x1=li.x1;
  }
  void init(Real _x0, Real _x1, const Array<T> &_a) {
#ifdef DEBUG
    if (_a.get_size()<3 && _a.get_size()!=1) ERRORQUIT("array size too small in PolyInterpolator");
#endif
    x0=_x0;
    x1=_x1;
    if (_a.get_size()==1) {
      a=_a;
    }
    else {
      a.resize(_a.get_size()+2);
      for (int i=0; i<_a.get_size(); i++) {
	a(1+i)=_a(i);
      }
      diff(&a(0),a(1),a(2));
      sum(&a(0),a(0),a(1));
      int n=a.get_size();
      diff(&a(n-1),a(n-2),a(n-3));
      sum(&a(n-1),a(n-1),a(n-2));
    }
  }
  void interpol(T *py, Real x) const {
    int n=a.get_size();
    if (n==1) {
      *py=a(0);
    }
    else {
      Real y=1+(Real)(n-3)*(x-x0)/(x1-x0);
      int iy=(int)floor(y)-1;
      if (iy<0) {iy=0;}
      if (iy>n-4) {iy=n-4;}
      Real fy=y-(Real)iy;
      Real c[4];
      c[0]= 2.0+(-4.0+( 2.5-0.5*fy)*fy)*fy;
      c[1]=-3.0+( 9.5+(-7.0+1.5*fy)*fy)*fy;
      c[2]= 3.0+(-8.0+( 6.5-1.5*fy)*fy)*fy;
      c[3]=-1.0+(+2.5+(-2.0+0.5*fy)*fy)*fy;
      T tmp;
      product(py,a(iy),c[0]);
      for (int i=1; i<4; i++) {
	product(&tmp,a(iy+i),c[i]);
	sum(py,*py,tmp);
      }
    }
  }
  Real min(void) const {return x0;}
  Real max(void) const {return x1;}
};

#endif
