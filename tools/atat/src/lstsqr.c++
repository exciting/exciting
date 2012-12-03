#include "lstsqr.h"

void calc_least_square_matrix(Array2d<Real> *answer, const Array2d<Real> &lhs) {
  iVector2d size=lhs.get_size();
  answer->resize(iVector2d(size(1),size(0)));
  Array2d<Real> to_invert;
  to_invert.resize(iVector2d(size(1),size(1)));
  iVector2d dst;
  for (dst(1)=0; dst(1)<size(1); dst(1)++) {
    for (dst(0)=0; dst(0)<size(1); dst(0)++) {
      to_invert(dst)=0.;
      for (int i=0; i<size(0); i++) {
        to_invert(dst)+=lhs(i,dst(0))*lhs(i,dst(1));
      }
    }
  }
  invert_matrix(&to_invert,to_invert);
  for (dst(1)=0; dst(1)<size(0); dst(1)++) {
    for (dst(0)=0; dst(0)<size(1); dst(0)++) {
      (*answer)(dst)=0.;
      for (int i=0; i<size(1); i++) {
        (*answer)(dst)+=to_invert(dst(0),i)*lhs(dst(1),i);
      }
    }
  }
}

int list_nonredundant_columns(Array<int> *which_col, const Array2d<Real> &x) {
  int nb_in_basis=0;
  Array<ArrayReal> basis;
  LinkedList<int> cols;
  for (int i=0; i<x.get_size()(1); i++) {
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
  LinkedList_to_Array(which_col,cols);
  return (which_col->get_size()==x.get_size()(1));
}

void calc_ols(Array<Real> *pb, const Array2d<Real> &x, const Array<Real> &y) {
  Array2d<Real> xx_i;
  inner_product(&xx_i, x,x);
  invert_matrix(&xx_i,xx_i);
  Array<Real> xy;
  inner_product(&xy, x,y);
  product(pb,xx_i,xy);
}

Array<Real> empty_rArray;

void calc_ols(Array<Real> *pb, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight, int nonredundant) {
  const Array2d<Real> *px=&x;
  const Array<Real> *py=&y;
  Array2d<Real> wx,nr_x;
  Array<Real> wy;
  if (weight.get_size()!=0) {
    product_diag(&wx,weight,*px);
    product_diag(&wy,weight,*py);
    px=&wx;
    py=&wy;
  }
  if (nonredundant) {
    Array<int> cols;
    list_nonredundant_columns(&cols,*px);
    extract_columns(&nr_x,*px,cols);
    px=&nr_x;
    Array<Real> b;
    calc_ols(&b,*px,*py);
    pb->resize(x.get_size()(1));
    zero_array(pb);
    for (int i=0; i<b.get_size(); i++) {
      (*pb)(cols(i))=b(i);
    }
  }
  else {
    calc_ols(pb,*px,*py);
  }
}

Real calc_cv(const Array2d<Real> &x, const Array<Real> &y) {
  int n=x.get_size()(0);
  int r=x.get_size()(1);
  Array2d<Real> xx_i;
  inner_product(&xx_i, x,x);
  invert_matrix(&xx_i,xx_i);
  Array<Real> xy;
  inner_product(&xy, x,y);
  Array<Real> b;
  product(&b,xx_i,xy);
  Array<Real> y_hat;
  product(&y_hat,x,b);
  Real cv=0.;
  for (int i=0; i<n; i++) {
    Real den=1;
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

Real calc_cv(const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight, int nonredundant) {
  const Array2d<Real> *px=&x;
  const Array<Real> *py=&y;
  Array2d<Real> wx,nr_x;
  Array<Real> wy;
  if (weight.get_size()!=0) {
    product_diag(&wx,weight,*px);
    product_diag(&wy,weight,*py);
    px=&wx;
    py=&wy;
  }
  if (nonredundant) {
    Array<int> cols;
    list_nonredundant_columns(&cols,*px);
    extract_columns(&nr_x,*px,cols);
    px=&nr_x;
  }
  return calc_cv(*px,*py);
}

void predict_ols(Array<Real> *p_yhat, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight, int nonredundant) {
  Array<Real> b;
  calc_ols(&b,x,y,weight,nonredundant);
  product(p_yhat, x,b);
}

void calc_ols_var(Array2d<Real> *pvar, const Array2d<Real> &x, const Array<Real> &y) {
  Array<Real> y_hat;
  predict_ols(&y_hat,x,y);
  Real sigma2=0.;
  for (int i=0; i<y_hat.get_size(); i++) {
    sigma2+=sqr((y_hat(i)-y(i)));
  }
  sigma2/=(Real)y_hat.get_size();
  inner_product(pvar, x,x);
  invert_matrix(pvar, *pvar);
  product(pvar, *pvar,sigma2);
}

void calc_gls_var(Array2d<Real> *pvar, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight) {
  Array<Real> y_hat;
  predict_ols(&y_hat,x,y,weight);
  Real sigma2=0.;
  for (int i=0; i<y_hat.get_size(); i++) {
    sigma2+=sqr((y_hat(i)-y(i))/weight(i));
  }
  sigma2/=(Real)y_hat.get_size();
  Array2d<Real> w_x;
  product_diag(&w_x, weight,x);
  inner_product(pvar, w_x,w_x);
  invert_matrix(pvar, *pvar);
  product(pvar, *pvar,sigma2);
}

void calc_gls_homo_var(Array2d<Real> *pvar, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &weight) {
  Array<Real> y_hat;
  predict_ols(&y_hat,x,y,weight);
  Real sigma2=0.;
  for (int i=0; i<y_hat.get_size(); i++) {
    sigma2+=sqr(y_hat(i)-y(i));
  }
  sigma2/=(Real)y_hat.get_size();
  Array2d<Real> w_x,ww_x;
  product_diag(&w_x, weight,x);
  product_diag(&ww_x, weight,w_x);
  Array2d<Real> xw2x,xw2xi,xw4x;
  inner_product(&xw2x, w_x,w_x);
  inner_product(&xw4x, ww_x,ww_x);
  invert_matrix(&xw2xi, xw2x);
  Array2d<Real> tmp;
  product(&tmp, xw4x,xw2xi);
  product(pvar, xw2xi, tmp);
  product(pvar, *pvar, sigma2);
}

int build_basis(Array<ArrayReal> *pbasis, int *pnb_in_basis, const Array<Real> &new_vector) {
  Array<ArrayReal> &basis=*pbasis;
  int &nb_in_basis=*pnb_in_basis;
  if (nb_in_basis==0) basis.resize(new_vector.get_size());
  Array2d<Real> linear_sys(nb_in_basis,nb_in_basis);
  for (int i=0; i<nb_in_basis; i++) {
    for (int j=0; j<=i; j++) {
      Real sum=inner_product(basis(i),basis(j));
      linear_sys(i,j)=sum;
      linear_sys(j,i)=sum;
    }
  }
  Array<Real> rhs(nb_in_basis);
  for (int i=0; i<nb_in_basis; i++) {
    rhs(i)=inner_product(basis(i),new_vector);
  }
  solve_linsys(&linear_sys,&rhs);
  Array<Real> projection(new_vector.get_size());
  for (int k=0; k<new_vector.get_size(); k++) {
    projection(k)=0.;
    for (int j=0; j<nb_in_basis; j++) {
      projection(k)+=rhs(j)*basis(j)(k);
    }
  }
  Real delta2=0.;
  for (int k=0; k<new_vector.get_size(); k++) {
    delta2+=sqr(projection(k)-new_vector(k));
  }
  if (sqrt(delta2)<zero_tolerance) {
    return 0;
  }
  else {
    basis(nb_in_basis).resize(new_vector.get_size());
    for (int k=0; k<new_vector.get_size(); k++) {
      basis(nb_in_basis)(k)=new_vector(k);
    }
    nb_in_basis++;
    return 1;
  }
}

Real detect_discontinuity(const Array<Real> &y, Real new_y, Real small_value) {
  if (y.get_size()<3) return 0.;
  int maxsize=20;
  if (y.get_size()>maxsize) {
    Array<Real> short_y(maxsize);
    int org=y.get_size()-maxsize;
    for (int i=0; i<maxsize; i++) {
      short_y(i)=y(org+i);
    }
    return detect_discontinuity(short_y,new_y,small_value);
  }
  Real best_cv=MAXFLOAT;
  int best_p=0;
  for (int p=0; p<y.get_size()-1; p++) {
    Array2d<Real> x(y.get_size(),p);
    for (int i=0; i<y.get_size(); i++) {
      for (int j=0; j<p; j++) {
        x(i,j)=pow((Real)i,(Real)j);
      }
    }
    Real cv=calc_cv(x,y);
    if (cv<best_cv) {
      best_cv=cv;
      best_p=p;
    }
  }
  Array2d<Real> x(y.get_size(),best_p);
  for (int i=0; i<y.get_size(); i++) {
    for (int j=0; j<best_p; j++) {
      x(i,j)=pow((Real)i,(Real)j);
    }
  }
  Array<Real> b;
  Array2d<Real> var;
  calc_ols(&b,x,y);
  calc_ols_var(&var, x,y);
  Array<Real> new_x(best_p);
  for (int j=0; j<best_p; j++) {
    new_x(j)=pow((Real)(y.get_size()),(Real)j);
  }
  Real predicted_new_y=inner_product(b,new_x);
  Array<Real> tmp;
  product(&tmp,var,new_x);
  Real v=inner_product(new_x,tmp);
  Array<Real> py,resid;
  product(&py,x,b);
  diff(&resid,y,py);
  Real sigma2=inner_product(resid,resid)/(Real)(y.get_size()-1);
  //  cerr << "v=" << v << "   sigma2=" << sigma2 << endl;
  Real t=fabs(new_y-predicted_new_y)/sqrt(v+sigma2+sqr(small_value));
  return t;
}

Real find_minimum(const Array<Real> &y, int mesh, Real *pos) {
  Real best_cv=MAXFLOAT;
  int best_p=0;
  for (int p=1; p<y.get_size(); p++) {
    Array2d<Real> x(y.get_size(),p);
    for (int i=0; i<y.get_size(); i++) {
      for (int j=0; j<p; j++) {
        x(i,j)=ipow((Real)i,j);
      }
    }
    Real cv=calc_cv(x,y);
    if (cv<best_cv) {
      best_cv=cv;
      best_p=p;
    }
  }
  if (best_p<3) {best_p=3;}
  Array2d<Real> x(y.get_size(),best_p);
  for (int i=0; i<y.get_size(); i++) {
    for (int j=0; j<best_p; j++) {
      x(i,j)=ipow((Real)i,j);
    }
  }
  Array<Real> b;
  calc_ols(&b,x,y);
  Real best_t;
  Real min_v=MAXFLOAT;
  if (best_p==3) {
    if (b(2)>0.) {
      best_t=-b(1)/(2.*b(2));
      min_v=b(0)+b(1)*best_t+b(2)*sqr(best_t);
    }
    else {
      cerr << "Warning: instability detected." << endl;
    }
  }
  if (min_v==MAXFLOAT) {
    for (int i=0; i<mesh; i++) {
      Real t=(Real)i*(y.get_size()-1)/(Real)mesh;
      Real v=0;
      for (int j=0; j<b.get_size(); j++) {v+=ipow((Real)t,j)*b(j);}
      if (v<min_v) {
	best_t=t;
	min_v=v;
      }
    }
  }
  if (pos) {(*pos)=best_t;}
  return min_v;
}

void gram_schmidt(Array<Array<Real> > *portho_basis, const Array<Array<Real> > &basis) {
  portho_basis->resize(basis.get_size());
  for (int i=0; i<basis.get_size(); i++) {
    if (basis(i).get_size()==0) break;
    for (int j=0; j<i; j++) {
      Real p=inner_product(((*portho_basis)(j)),((*portho_basis)(i)));
      Array<Real> v;
      product(&v, (*portho_basis)(j),p);
      diff(&((*portho_basis)(i)), ((*portho_basis)(i)), v);
    }
    product(&((*portho_basis)(i)), basis(i),1./sqrt(inner_product((*portho_basis)(i),(*portho_basis)(i))));
  }
}

void calc_perfect_fit(Array<Real> *pb, const Array2d<Real> &x, const Array<Real> &y, const Array<Real> &w) {
  Array2d<Real> xw,xwxt;
  product_diag(&xw, x,w);
  outer_product(&xwxt, x,xw);
  invert_matrix(&xwxt,xwxt);
  Array<Real> xwxtiy;
  product(&xwxtiy,xwxt,y);
  inner_product(pb,xw,xwxtiy);
}
