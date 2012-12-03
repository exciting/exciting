#ifndef __CHULL_H__
#define __CHULL_H__

#include <iostream.h>
#include "arraylist.h"
#include "stringo.h"

class PolytopeFace {
public:
  Array<int> pts;
  Array<Real> up;
  Real c;
};

class LinearInequality {
  public:
   Array<Real> v;
   Real c;
  LinearInequality(int size=0): v(size) {c=0;}
  LinearInequality(const LinearInequality &a): v(a.v) {c=a.c;}
  LinearInequality(const Array<Real> &_v, Real _c): v(_v) {c=_c;}
};

void de_mean(Array<Array<Real> > *px, const Array<Array<Real> > &org_x);

void calc_convex_hull(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &org_x, const Array<Real> &ground);

void calc_convex_hull_degenerate(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &org_x, const Array<Real> &ground);

void update_normals(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &org_x);

int flag_outside_hull(Array<int> *poutside, LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, int flag_plane_too);

int is_point_in_hull(const Array<Real> &x, const LinkedList<LinearInequality> &ineq_list);

void clip_hull(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, const LinkedList<LinearInequality> &ineq);

void calc_formation(Array<Real> *pfe, Array<Real> *ppure, const Array<Array<Real> > &x, const Array<Real> &e);

void read_inequalities(LinkedList<LinearInequality> *ineq_list, const Array<AutoString> &label, istream &s);
void read_equalities(LinkedList<LinearInequality> *ineq_list, const Array<AutoString> &label, istream &s);
void ineq_list_to_Array(Array2d<Real> *a, Array<Real> *c, const LinkedList<LinearInequality> &l);

// dim n+1 versions;

void calc_convex_hull_p1(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, const Array<Real> &e);

void update_normals_p1(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, const Array<Real> &e);

int flag_outside_hull_p1(Array<int> *poutside, LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, const Array<Real> &e, int flag_plane_too);

template<class T>
void paste_row_vector(Array<Array<T> > *pout, const Array<Array<T> > &x, const Array<Array<T> > &y) {
  int total_dim=x(0).get_size()+y(0).get_size();
  pout->resize(x.get_size(),total_dim);
  for (int i=0; i<pout->getsize(); i++) {
    int j=0;
    for (int k=0; k<x(i).get_size(); k++) {
      (*pout)(i)(j)=x(i)(k);
      j++;
    }
    for (int k=0; k<y(i).get_size(); k++) {
      (*pout)(i)(j)=y(i)(k);
      j++;
    }
  }
}

template<class T>
void paste_col_vector(Array<Array<T> > *pout, const Array<Array<T> > &x, const Array<T> &y) {
  if (pout->get_size()>0) {
    if ((*pout)(0).get_size()==x(0).get_size()) {
      int j=x(0).get_size()-1;
      for (int i=0; i<pout->get_size(); i++) {
	(*pout)(i)(j)=y(i);
      }
      return;
    }
  }
  pout->resize(x.get_size());
  for (int i=0; i<pout->get_size(); i++) {
    (*pout)(i).resize(x(0).get_size()+1);
    int j=0;
    for (j=0; j<x(i).get_size(); j++) {
      (*pout)(i)(j)=x(i)(j);
    }
    (*pout)(i)(j)=y(i);
  }
}

#endif
