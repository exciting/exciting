#include "linalg.h"

void transpose_matrix(Array2d<Real> *trans,const Array2d<Real> &a) {
  int maxi=a.get_size()(1);
  int maxj=a.get_size()(0);
  trans->resize(maxi,maxj);
  for (int i=0; i<maxi; i++) {
    for (int j=0; j<maxj; j++) {
      (*trans)(i,j)=a(j,i);
    }
  }
}

void diff(Array<Real> *result, const Array<Real> &a, const Array<Real> &b) {
#ifdef DEBUG
  if (a.get_size()!=b.get_size()) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(a.get_size());
  for (int i=0; i<a.get_size(); i++) {
    (*result)(i)=a(i)-b(i);
  }
}

void diff(Array2d<Real> *result, const Array2d<Real> &a, const Array2d<Real> &b) {
#ifdef DEBUG
  if (a.get_size()!=b.get_size()) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(a.get_size());
  for (int i=0; i<a.get_size()(0); i++) {
    for (int j=0; j<a.get_size()(1); j++) {
      (*result)(i,j)=a(i,j)-b(i,j);
    }
  }
}

void sum(Array<Real> *result, const Array<Real> &a, const Array<Real> &b) {
#ifdef DEBUG
  if (a.get_size()!=b.get_size()) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(a.get_size());
  for (int i=0; i<a.get_size(); i++) {
    (*result)(i)=a(i)+b(i);
  }
}

Real mean(const Array<Real> &a) {
    Real m=0.;
    for (int i=0; i<a.get_size(); i++) {
	m+=a(i);
    }
    return (m/(Real)(a.get_size()));
}

void sum(Array2d<Real> *result, const Array2d<Real> &a, const Array2d<Real> &b) {
#ifdef DEBUG
  if (a.get_size()!=b.get_size()) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(a.get_size());
  for (int i=0; i<a.get_size()(0); i++) {
    for (int j=0; j<a.get_size()(1); j++) {
      (*result)(i,j)=a(i,j)+b(i,j);
    }
  }
}

void sum(Array2d<Complex> *result, const Array2d<Complex> &a, const Array2d<Complex> &b) {
#ifdef DEBUG
  if (a.get_size()!=b.get_size()) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(a.get_size());
  for (int i=0; i<a.get_size()(0); i++) {
    for (int j=0; j<a.get_size()(1); j++) {
      (*result)(i,j)=a(i,j)+b(i,j);
    }
  }
}

void product(Array2d<Real> *result, const Array2d<Real> &a, const Array2d<Real> &b) {
#ifdef DEBUG
  if (a.get_size()(1)!=b.get_size()(0)) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(iVector2d(a.get_size()(0),b.get_size()(1)));
  zero_array(result);
  for (int i=0; i<a.get_size()(0); i++) {
    for (int j=0; j<b.get_size()(1); j++) {
      for (int k=0; k<a.get_size()(1); k++) {
        (*result)(i,j)+=a(i,k)*b(k,j);
      }
    }
  }
}

void product_diag(Array2d<Real> *result, const Array<Real> &diag, const Array2d<Real> &b) {
#ifdef DEBUG
  if (diag.get_size()!=b.get_size()(0)) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(b.get_size());
  for (int i=0; i<b.get_size()(0); i++) {
    for (int j=0; j<b.get_size()(1); j++) {
      (*result)(i,j)=diag(i)*b(i,j);
    }
  }
}

void product_diag(Array2d<Real> *result, const Array2d<Real> &b, const Array<Real> &diag) {
#ifdef DEBUG
  if (diag.get_size()!=b.get_size()(1)) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(b.get_size());
  for (int i=0; i<b.get_size()(0); i++) {
    for (int j=0; j<b.get_size()(1); j++) {
      (*result)(i,j)=b(i,j)*diag(j);
    }
  }
}

void product_diag(Array<Real> *result, const Array<Real> &diag, const Array<Real> &b) {
#ifdef DEBUG
  if (diag.get_size()!=b.get_size()) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(b.get_size());
  for (int i=0; i<b.get_size(); i++) {
    (*result)(i)=diag(i)*b(i);
  }
}

void product(Array<Real> *result, const Array2d<Real> &a, const Array<Real> &b) {
#ifdef DEBUG
  if (a.get_size()(1)!=b.get_size()) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(a.get_size()(0));
  zero_array(result);
  for (int i=0; i<a.get_size()(0); i++) {
    for (int k=0; k<a.get_size()(1); k++) {
      (*result)(i)+=a(i,k)*b(k);
    }
  }
}

void product(Array<Real> *result, const Array<Real> &a, const Array2d<Real> &b) {
#ifdef DEBUG
  if (a.get_size()!=b.get_size()(0)) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(b.get_size()(1));
  zero_array(result);
  for (int i=0; i<b.get_size()(1); i++) {
    for (int k=0; k<a.get_size(); k++) {
      (*result)(i)+=a(k)*b(k,i);
    }
  }
}

void product(Array2d<Real> *result, const Array2d<Real> &a, Real b) {
  result->resize(a.get_size());
  for (int i=0; i<a.get_size()(0); i++) {
    for (int j=0; j<a.get_size()(1); j++) {
      (*result)(i,j)=a(i,j)*b;
    }
  }
}

void product(Array<Real> *result, const Array<Real> &a, Real b) {
  result->resize(a.get_size());
  for (int i=0; i<a.get_size(); i++) {
    (*result)(i)=a(i)*b;
  }
}

Real inner_product(const Array<Real> &a, const Array<Real> &b) {
  Real sum=0;
#ifdef DEBUG
  if (a.get_size()!=b.get_size()) ERRORQUITDUMP("Matrices not conformable");
#endif
  for (int i=0; i<a.get_size(); i++) {
    sum+=a(i)*b(i);
  }
  return sum;  
}

void inner_product(Array2d<Real> *result, const Array2d<Real> &a, const Array2d<Real> &b) {
#ifdef DEBUG
  if (a.get_size()(0)!=b.get_size()(0)) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(iVector2d(a.get_size()(1),b.get_size()(1)));
  zero_array(result);
  int nk=a.get_size()(0);
  for (int i=0; i<a.get_size()(1); i++) {
    for (int j=0; j<b.get_size()(1); j++) {
      Real s=0.;
      for (int k=0; k<nk; k++) {
        s+=a(k,i)*b(k,j);
      }
      (*result)(i,j)=s;
    }
  }
}

void inner_product(Array<Real> *result, const Array2d<Real> &a, const Array<Real> &b) {
#ifdef DEBUG
  if (a.get_size()(0)!=b.get_size()) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(a.get_size()(1));
  zero_array(result);
  for (int i=0; i<a.get_size()(1); i++) {
    for (int k=0; k<a.get_size()(0); k++) {
      (*result)(i)+=a(k,i)*b(k);
    }
  }
}

Real normalize(Array<Real> *a) {
  Real l=sqrt(inner_product(*a,*a));
  product(a,*a,1./l);
  return l;
}

void outer_product(Array2d<Real> *result, const Array2d<Real> &a, const Array2d<Real> &b) {
#ifdef DEBUG
  if (a.get_size()(1)!=b.get_size()(1)) ERRORQUITDUMP("Matrices not conformable");
#endif
  result->resize(iVector2d(a.get_size()(0),b.get_size()(0)));
  zero_array(result);
  for (int i=0; i<a.get_size()(0); i++) {
    for (int j=0; j<b.get_size()(0); j++) {
      for (int k=0; k<a.get_size()(1); k++) {
        (*result)(i,j)+=a(i,k)*b(j,k);
      }
    }
  }
}

void outer_product(Array2d<Real> *result, const Array<Real> &a, const Array<Real> &b) {
  result->resize(iVector2d(a.get_size(),b.get_size()));
  for (int i=0; i<a.get_size(); i++) {
    for (int j=0; j<b.get_size(); j++) {
      (*result)(i,j)=a(i)*b(j);
    }
  }
}

void outer_product_diag(Array<Real> *result, const Array2d<Real> &a, const Array2d<Real> &b) {
#ifdef DEBUG
  if (a.get_size()(1)!=b.get_size()(1)) ERRORQUITDUMP("Matrices not conformable");
  if (a.get_size()(0)!=b.get_size()(0)) ERRORQUITDUMP("Matrix is not square");
#endif
  result->resize(a.get_size()(0));
  zero_array(result);
  for (int i=0; i<a.get_size()(0); i++) {
    for (int k=0; k<a.get_size()(1); k++) {
      (*result)(i)+=a(i,k)*b(i,k);
    }
  }
}

Real quadratic_form(const Array2d<Real> &a, const Array<Real> &x) {
#ifdef DEBUG
  if (a.get_size()(0)!=a.get_size()(1)) ERRORQUITDUMP("Matrix is not square");
  if (a.get_size()(0)!=x.get_size())    ERRORQUITDUMP("Matrices not conformable");
#endif
  Real result=0.;
  for (int i=0; i<x.get_size(); i++) {
    for (int j=0; j<x.get_size(); j++) {
      result+=a(i,j)*x(i)*x(j);
    }
  }
  return result;
}

Real quadratic_form(const Array2d<Complex> &a, const Array<Complex> &x) {
#ifdef DEBUG
  if (a.get_size()(0)!=a.get_size()(1)) ERRORQUITDUMP("Matrix is not square");
  if (a.get_size()(0)!=x.get_size())    ERRORQUITDUMP("Matrices not conformable");
#endif
  Real result=0.;
  for (int i=0; i<x.get_size(); i++) {
    for (int j=0; j<x.get_size(); j++) {
      result+=real((a(i,j)*conj(x(i))*x(j)));
    }
  }
  return result;
}

Real trace(const Array2d<Real> &a) {
#ifdef DEBUG
  if (a.get_size()(0)!=a.get_size()(1)) ERRORQUITDUMP("Cannot take trace of nonsquare matrix");
#endif
  Real tr=0;
  for (int i=0; i<a.get_size()(0); i++) {
    tr+=a(i,i);
  }
  return tr;
}

void set_identity(Array2d<Real> *result, int n) {
  result->resize(iVector2d(n,n));
  zero_array(result);
  for (int i=0; i<n; i++) {
    (*result)(i,i)=1.;
  }
}

Real max_norm(const Array<Real> &v) {
  Real m=0.;
  for (int i=0; i<v.get_size(); i++) {
    m=MAX(m,fabs(v(i)));
  }
  return m;
}

void invert_matrix_tall(Array2d<Real> *inv, const Array2d<Real> &a) {
  Array2d<Real> ata,iata;
  inner_product(&ata,a,a);
  invert_matrix(&iata,ata);
  outer_product(inv,iata,a);
}

void invert_matrix_wide(Array2d<Real> *inv, const Array2d<Real> &a) {
  Array2d<Real> aat,iaat;
  outer_product(&aat,a,a);
  invert_matrix(&iaat,aat);
  inner_product(inv,a,iaat);
}
