#include "chull.h"
#include "arraylist.h"
#include "linalg.h"

class MultiDimPoly {
  int dim;
  Array<Real> maxpow;
  LinkedList<LinearInequality> bnd;
public:
  MultiDimPoly(void): maxpow(), bnd() {dim=0;}
  MultiDimPoly(const Array<Real> &_maxpow, const LinkedList<LinearInequality> &_bnd);
  void init(const Array<Real> &_maxpow, const LinkedList<LinearInequality> &_bnd);
  Real eval(const Array<Real> &a, const Array<Real> &x) const;
  void eval_powers(Array<Real> *ppowers, const Array<Real> &x) const;
  int get_dim_param(void) const {
    return dim;
  }
  int get_dim_arg(void) const {
    return maxpow.get_size();
  }
};
