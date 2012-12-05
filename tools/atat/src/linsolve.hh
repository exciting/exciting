#ifndef _LINSOLVE_H_
#define _LINSOLVE_H_

#include <math.h>
#include <complex.h>
#include "misc.h"
#include "array.h"

void lu_decomposition(Array2d<Real> *pa, Array<int> *ppivot);
void lu_backsub(Array2d<Real> *pa, Array<int> *ppivot, Array<Real> *pb);
void solve_linsys(Array2d<Real> *pa, Array<Real> *pb);
void solve_linsys(Array<Real> *psoln, const Array2d<Real> &mat, const Array<Real> &vect);
void invert_matrix(Array2d<Real> *inv, const Array2d<Real> &a);

/* eigenvalue problems */
void to_tridiag(Array2d<Real> *pa, Array<Real> *pd, Array<Real> *pe, int do_vect);
void tridiag_to_diag(Array<Real> *pd, Array<Real> *pe, Array2d<Real> *pz, int do_vect);
void diagonalize_symmetric_matrix(Array<Real> *p_lambda, Array2d<Real> *p_vect, const Array2d<Real> &mat);
void diagonalize_symmetric_matrix(Array<Real> *p_lambda, Array2d<Complex> *p_vect, const Array2d<Complex> &mat);

#endif
