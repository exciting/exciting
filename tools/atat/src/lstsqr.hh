#ifndef _LSTSQR_H_
#define _LSTSQR_H_

#include "linalg.h"
#include "integer.h"
#include "arraylist.h"

void calc_least_square_matrix(Array2d<Real> *answer, const Array2d<Real> &lhs);
int list_nonredundant_columns(Array<int> *which_col, const Array2d<Real> &x);

extern Array<Real> empty_rArray;

void calc_ols(Array<Real> *pb, const Array2d<Real> &x, const Array<Real> &y);
void calc_ols(Array<Real> *pb, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight, int nonredundant=0);
Real calc_cv(const Array2d<Real> &x, const Array<Real> &y);
Real calc_cv(const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight, int nonredundant=0);
void predict_ols(Array<Real> *p_yhat, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight=empty_rArray, int nonredundant=0);

void calc_ols_var(Array2d<Real> *pvar, const Array2d<Real> &x, const Array<Real> &y);
void calc_gls_var(Array2d<Real> *pvar, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight);
void calc_gls_homo_var(Array2d<Real> *pvar, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight);
int build_basis(Array<ArrayReal> *pbasis, int *pnb_in_basis, const Array<Real> &new_vector);
void gram_schmidt(Array<Array<Real> > *portho_basis, const Array<Array<Real> > &basis);
Real detect_discontinuity(const Array<Real> &y, Real new_y, Real small_value=0.001);

Real find_minimum(const Array<Real> &y, int mesh=100, Real *pos=NULL);

void calc_perfect_fit(Array<Real> *pb, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &w);

#endif
