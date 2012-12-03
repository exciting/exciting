#include "array.h"

int list_nonredundant_columns_regul(Array<int> *which_col, const Array2d<Real> &x, const Array<Real> &r);
void calc_ols_regul(Array<Real> *pb, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &r);
void calc_ols_regul(Array<Real> *pb, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight, int nonredundant, const Array<Real> &r);
Real calc_cv_regul(const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &reg);
Real calc_cv_regul(const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight, int nonredundant, const Array<Real> &r);
void predict_ols_regul(Array<Real> *p_yhat, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight, int nonredundant, const Array<Real> &r);

