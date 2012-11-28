#ifndef __TENSOR_H__
#define __TENSOR_H__

#include <math.h>
#include "misc.h"
#include "linalg.h"

template<class T>
class Tensor {
protected:
  Array<int> dim;
  Array<T> x;
public:
  Tensor(void): dim(), x() {}
  Tensor(const Tensor<T> &a): dim(a.dim), x(a.x) {}
  void operator=(const Tensor<T> &a) {
    dim=a.dim;
    x=a.x;
  }
  Tensor(const Array<int> &_dim): dim(), x() {
    resize(_dim);
  }
  void resize(const Array<int> &_dim) {
    dim=_dim;
    int m=1;
    for (int i=0; i<dim.get_size(); i++) {m*=dim(i);}
    x.resize(m);
  }
  void zero(void) {
    zero_array(&x);
  }
  T& operator()(const Array<int> &k) {
    int j=0;
    for (int i=0; i<dim.get_size(); i++) {
      j=j*dim(i)+k(i);
    }
    return x(j);
  }
  const T& operator()(const Array<int> &k) const {
    int j=0;
    for (int i=0; i<dim.get_size(); i++) {
      j=j*dim(i)+k(i);
    }
    return x(j);
  }
  const Array<int> & get_size(void) const {return dim;}
  void operator+=(const Tensor<T> &a) {
    sum(&x,x,a.x);
  }
  void operator-=(const Tensor<T> &a) {
    diff(&x,x,a.x);
  }
  void operator*=(Real a) {
    for (int i=0; i<x.get_size(); i++) {
      x(i)*=a;
    }
  }
  void operator/=(Real a) {
    for (int i=0; i<x.get_size(); i++) {
      x(i)/=a;
    }
  }
  friend Real norm(const Tensor<T> &t) {
    Real sum=0.;
    for (int i=0; i<t.x.get_size(); i++) {
      sum+=sqr(t.x(i));
    }
    return sqrt(sum);
  }
  T * get_buf(void) {
    return x.get_buf();
  }
  Array<T> & vectorize(void) {
    return x;
  }
  Array<T> const & vectorize(void) const {
    return x;
  }
};

#define rTensor Tensor<Real>

#endif
