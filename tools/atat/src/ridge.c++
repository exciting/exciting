#include "ridge.h"
#include "arraylist.h"
#include "lstsqr.h"

int list_nonredundant_columns_regul(Array<int> *which_col, const Array2d<Real> &x, const Array<Real> &r) {
  int nb_in_basis=0;
  Array<ArrayReal> basis;
  LinkedList<int> cols;
  for (int i=0; i<x.get_size()(1); i++) {
    if (r(i)!=0) {
      cols << new int(i);
    }
    else {
      Array<Real> cur_col;
      extract_column(&cur_col,x,i);
      if (build_basis(&basis, &nb_in_basis, cur_col)) {
	cols << new int(i);
	/*
	int j;
	for (j=0; j<i; j++) {
	  int nb_diff=0;
	  for (int l=0; l<x.get_size()(0); l++) {
	    if (fabs(x(l,j)-x(l,i))>zero_tolerance) nb_diff++;
	  }
	  if (nb_diff<=1) break;
	}
	if (j==i) {
	  cols << new int(i);
	}
	*/
      }
    }
  }
  LinkedList_to_Array(which_col,cols);
  return (which_col->get_size()==x.get_size()(1));
}

void calc_ols_regul(Array<Real> *pb, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &r) {
  Array2d<Real> xx_i;
  inner_product(&xx_i, x,x);
  for (int j=0; j<r.get_size(); j++) {
    xx_i(j,j)+=r(j);
  }
  invert_matrix(&xx_i,xx_i);
  Array<Real> xy;
  inner_product(&xy, x,y);
  product(pb,xx_i,xy);
}

void calc_ols_regul(Array<Real> *pb, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight, int nonredundant, const Array<Real> &r) {
  const Array2d<Real> *px=&x;
  const Array<Real> *py=&y;
  const Array<Real> *pr=&r;
  Array2d<Real> wx,nr_x;
  Array<Real> nr_r;
  Array<Real> wy;
  if (weight.get_size()!=0) {
    product_diag(&wx,weight,*px);
    product_diag(&wy,weight,*py);
    px=&wx;
    py=&wy;
  }
  if (nonredundant) {
    Array<int> cols;
    list_nonredundant_columns_regul(&cols,*px,r);
    extract_columns(&nr_x,*px,cols);
    px=&nr_x;
    extract_elements(&nr_r,*pr,cols);
    pr=&nr_r;
    Array<Real> b;
    calc_ols_regul(&b,*px,*py,*pr);
    pb->resize(x.get_size()(1));
    zero_array(pb);
    for (int i=0; i<b.get_size(); i++) {
      (*pb)(cols(i))=b(i);
    }
  }
  else {
    calc_ols_regul(pb,*px,*py,*pr);
  }
}

Real calc_cv_regul(const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &reg) {
  int n=x.get_size()(0);
  int r=x.get_size()(1);
  Array2d<Real> xx_i;
  inner_product(&xx_i, x,x);
  for (int i=0; i<r; i++) {
    xx_i(i,i)+=reg(i);
  }
  invert_matrix(&xx_i,xx_i);
  Array<Real> xy;
  inner_product(&xy, x,y);
  Array<Real> b;
  product(&b,xx_i,xy);
  Array<Real> y_hat;
  product(&y_hat,x,b);
  Real cv=0.;
  for (int i=0; i<n; i++) {
    Real den=1.;
    for (int j=0; j<r; j++) {
      for (int k=0; k<r; k++) {
        den-=x(i,j)*xx_i(j,k)*x(i,k);
      }
    }
    if (den<zero_tolerance) {
      return MAXFLOAT;
    }
    cv+=sqr((y_hat(i)-y(i))/den);
  }
  return sqrt(cv/(Real)n);
}

Real calc_cv_regul(const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight, int nonredundant, const Array<Real> &r) {
  const Array2d<Real> *px=&x;
  const Array<Real> *py=&y;
  const Array<Real> *pr=&r;
  Array2d<Real> wx,nr_x;
  Array<Real> nr_r;
  Array<Real> wy;
  if (weight.get_size()!=0) {
    product_diag(&wx,weight,*px);
    product_diag(&wy,weight,*py);
    px=&wx;
    py=&wy;
  }
  if (nonredundant) {
    Array<int> cols;
    list_nonredundant_columns_regul(&cols,*px,r);
    extract_columns(&nr_x,*px,cols);
    px=&nr_x;
    extract_elements(&nr_r,*pr,cols);
    pr=&nr_r;
  }
  return calc_cv_regul(*px,*py,*pr);
}

void predict_ols_regul(Array<Real> *p_yhat, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight, int nonredundant, const Array<Real> &r) {
  Array<Real> b;
  calc_ols_regul(&b,x,y,weight,nonredundant,r);
  product(p_yhat, x,b);
}

