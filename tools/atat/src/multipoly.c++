#include "multipoly.h"

MultiDimPoly::MultiDimPoly(const Array<Real> &_maxpow, const LinkedList<LinearInequality> &_bnd): maxpow(), bnd() {
  init(_maxpow,_bnd);
}

void MultiDimPoly::init(const Array<Real> &_maxpow, const LinkedList<LinearInequality> &_bnd) {
  maxpow=_maxpow;
  add_copy(&bnd,_bnd);
  Array<Real> zeros(maxpow.get_size());
  zero_array(&zeros);
  MultiDimIterator<Array<Real> > p(zeros,maxpow);
  dim=0;
  for (; p; p++) {
    LinkedListIterator<LinearInequality> ibnd(bnd);
    int keep=1;
    for (; ibnd; ibnd++) {
      if (inner_product((Array<Real> &)p,ibnd->v) > ibnd->c) {
	keep=0;
	break;
      }
    }
    if (keep) {
      //cout << (Array<Real> &)p << endl;
      dim++;
    }
  }
}

Real MultiDimPoly::eval(const Array<Real> &a, const Array<Real> &x) const {
  Array<Real> zeros(maxpow.get_size());
  zero_array(&zeros);
  MultiDimIterator<Array<Real> > p(zeros,maxpow);
  Real acc=0.;
  int i=0;
  for (; p; p++) {
    Array<Real> &ap=(Array<Real> &)p;
    int keep=1;
    LinkedListIterator<LinearInequality> ibnd(bnd);
    for (; ibnd; ibnd++) {
      if ( inner_product(ap,ibnd->v) > ibnd->c) {
	keep=0;
	break;
      }
    }
    if (keep) {
      Real xp=1.;
      for (int j=0; j<ap.get_size(); j++) {
	if (ap(j)>0) {xp*=pow(x(j),ap(j));}
      }
      acc+=a(i)*xp;
      i++;
    }
  }
  return acc;
}

void MultiDimPoly::eval_powers(Array<Real> *ppowers, const Array<Real> &x) const {
  Array<Real> zeros(maxpow.get_size());
  zero_array(&zeros);
  MultiDimIterator<Array<Real> > p(zeros,maxpow);
  ppowers->resize(get_dim_param());
  int i=0;
  for (; p; p++) {
    Array<Real> &ap=(Array<Real> &)p;
    int keep=1;
    LinkedListIterator<LinearInequality> ibnd(bnd);
    for (; ibnd; ibnd++) {
      if ( inner_product(ap,ibnd->v) > ibnd->c) {
        keep=0;
        break;
      }
    }
    if (keep) {
      Real xp=1.;
      for (int j=0; j<ap.get_size(); j++) {
        if (ap(j)>0) {xp*=pow(x(j),ap(j));}
      }
      (*ppowers)(i)=xp;
      i++;
    }
  }
}

