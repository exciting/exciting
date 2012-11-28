#include "linsolve.h"

void lu_decomposition(Array2d<Real> *pa, Array<int> *ppivot) {
  const Real tiny=1.0e-16;
  Array2d<Real> &a=*pa;
  int n=a.get_size()(0);
  if (ppivot->get_size()==0) {ppivot->resize(n);}
  
  Array<Real> w(n);
  for (int i=0; i<n; i++) {
    Real maxel=0.;
    for (int j=0; j<n; j++) {maxel=MAX(maxel,fabs(a(i,j)));}
    if (maxel==0.0) {
      ERRORQUIT("Singular matrix in lu_decomposition");
    }
    w(i)=1./maxel;
  }
  for (int j=0; j<n; j++) {
    for (int i=0; i<j; i++) {
      Real acc=a(i,j);
      for (int k=0; k<i; k++) {acc-=a(i,k)*a(k,j);}
      a(i,j)=acc;
    }
    Real maxel=0.;
    int imax;
    for (int i=j; i<n; i++) {
      Real acc=a(i,j);
      for (int k=0; k<j; k++) {acc-=a(i,k)*a(k,j);}
      a(i,j)=acc;
      Real dummy=w(i)*fabs(acc);
      if (dummy>=maxel) {
	maxel=dummy;
	imax=i;
      }
    }
    if (j!=imax) {
      for (int k=0; k<n; k++) {swap(&a(imax,k),&a(j,k));}
      w(imax)=w(j);
    }
    (*ppivot)(j)=imax;
    if (a(j,j)==0.) a(j,j)=tiny;
    if (j!=n) {
      for (int i=j+1; i<n; i++) {a(i,j)/=a(j,j);}
    }
  }
}

void lu_backsub(Array2d<Real> *pa, Array<int> *ppivot, Array<Real> *pb) {
  Array2d<Real> &a=*pa;
  int n=a.get_size()(0);
  Array<Real> &b=*pb;
  
  int ii=n;
  for (int i=0; i<n; i++) {
    int ip=(*ppivot)(i);
    Real acc=b(ip);
    b(ip)=b(i);
    for (int j=ii; j<i; j++) {acc-=a(i,j)*b(j);}
    if (ii==n && fabs(acc)>0.) {ii=i;}
    b(i)=acc;
  }
  for (int i=n-1; i>=0; i--) {
    Real acc=b(i);
    for (int j=i+1; j<n; j++) {acc-=a(i,j)*b(j);}
    b(i)=acc/a(i,i);
  }
}

void solve_linsys(Array2d<Real> *pa, Array<Real> *pb) {
  Array<int> pivot;
  lu_decomposition(pa,&pivot);
  lu_backsub(pa,&pivot,pb);
}

void solve_linsys(Array<Real> *psoln, const Array2d<Real> &mat, const Array<Real> &vect) {
  int n=vect.get_size();
  Array2d<Real> a(mat);
  Array<Real> b(vect);
  solve_linsys(&a,&b);
  (*psoln)=b;
}

void invert_matrix(Array2d<Real> *inv, const Array2d<Real> &a) {
  iVector2d size=a.get_size();
#ifdef DEBUG
  if (size(0)!=size(1)) ERRORQUIT("Cannot invert non square matrix!");
#endif
  inv->resize(size);
  int n=size(0);
  Array2d<Real> linear_sys(a);
  Array<int> pivot;
  lu_decomposition(&linear_sys,&pivot);
  Array<Real> rhs(n);
  for (int col=0; col<n; col++) {
    zero_array(&rhs);
    rhs(col)=1.;
    lu_backsub(&linear_sys,&pivot,&rhs);
    for (int i=0; i<n; i++) {
      (*inv)(i,col)=rhs(i);
    }
  }
}

/* eigenvalue problems */

void to_tridiag(Array2d<Real> *pa, Array<Real> *pd, Array<Real> *pe, int do_vect) {
  Array2d<Real> &a=*pa;
  int n=a.get_size()(0);
  Array<Real> &d=*pd;
  d.resize(n);
  Array<Real> &e=*pe;
  e.resize(n);
  
  for (int i=n-1; i>0; i--) {
    int l=i-1;
    Real h=0.;
    if (l==0) {
      e(i)=a(i,l);
    } else {
      Real nrm=0.;
      for (int k=0; k<=l; k++) {
	nrm+=fabs(a(i,k));
      }
      if (nrm==0.) {
	e(i)=a(i,l);
      } else {
	for (int k=0; k<=l; k++) {
	  a(i,k)/=nrm;
	  h+=a(i,k)*a(i,k);
	}
	Real f=a(i,l);
	Real g=(f>0. ? -sqrt(h) : sqrt(h));
	e(i)=nrm*g;
	h-=f*g;
	a(i,l)=f-g;
	f=0.;
	for (int j=0; j<=l; j++) {
	  a(j,i)=a(i,j)/h;
	  g=0.;
	  for (int k=0; k<=j; k++) {
	    g+=a(j,k)*a(i,k);
	  }
	  for (int k=j+1; k<=l; k++) {
	    g+=a(k,j)*a(i,k);
	  }
	  e(j)=g/h;
	  f+=e(j)*a(i,j);
	}
	Real fhh=f/(h+h);
	for (int j=0; j<=l; j++) {
	  f=a(i,j);
	  g=e(j)-fhh*f;
	  e(j)=g;
	  for (int k=0; k<=j; k++) {
	    a(j,k)-=(f*e(k)+g*a(i,k));
	  }
	}
      }
    }
    d(i)=h;
  }
  d(0)=0.;
  e(0)=0.;
  if (do_vect) {
    for (int i=0; i<n; i++) {
      int l=i-1;
      if (d(i)) {
	for (int j=0; j<=l; j++) {
	  Real g=0.;
	  for (int k=0; k<=l; k++) {
	    g+=a(i,k)*a(k,j);
	  }
	  for (int k=0; k<=l; k++) {
	    a(k,j)-=g*a(k,i);
	  }
	}
      }
      d(i)=a(i,i);
      a(i,i)=1.;
      for (int j=0; j<=l; j++) {a(j,i)=0.; a(i,j)=0.;}
    }
  }
  else {
    for (int i=0; i<n; i++) {
      d(i)=a(i,i);
    }
  }
}

void tridiag_to_diag(Array<Real> *pd, Array<Real> *pe, Array2d<Real> *pa, int do_vect) {
  Array<Real> &d=*pd;
  int n=d.get_size();
  Array<Real> &e=*pe;
  Array2d<Real> &a=*pa;
  
  for (int i=1; i<n; i++) {e(i-1)=e(i);}
  e(n-1)=0.;
  int m;
  for (int j=0; j<n; j++) {
    do {
      for (m=j; m<n-1; m++) {
	Real dd=fabs(d(m))+fabs(d(m+1));
	if (fabs(e(m))+dd==dd) break;
      }
      if (m!=j) {
	Real g=(d(j+1)-d(j))/(2.*e(j));
	Real r=sqrt((g*g)+1.);
	g=d[m]-d[j]+e[j]/(g+(g<0. ? -fabs(r) : fabs(r)));
	Real s=1.;
	Real c=1.;
	Real p=0.;
	for (int i=m-1; i>=j; i--) {
	  Real f=s*e(i);
	  Real b=c*e(i);
	  if (fabs(f)>=fabs(g)) {
	    c=g/f;
	    r=sqrt((c*c)+1.);
	    e(i+1)=f*r;
	    s=1./r;
	    c*=s;
	  } else {
	    s=f/g;
	    r=sqrt((s*s)+1.);
	    e(i+1)=g*r;
	    c=1./r;
	    s*=c;
	  }
	  g=d(i+1)-p;
	  r=(d(i)-g)*s+2.*c*b;
	  p=s*r;
	  d(i+1)=g+p;
	  g=c*r-b;
	  if (do_vect) {
	    for (int k=0; k<n; k++) {
	      f=a(k,i+1);
	      a(k,i+1)=s*a(k,i)+c*f;
	      a(k,i)=c*a(k,i)-s*f;
	    }
	  }
	}
	d(j)=d(j)-p;
	e(j)=g;
	e(m)=0.;
      }
    } while (m!=j);
  }
}

void diagonalize_symmetric_matrix(Array<Real> *p_lambda, Array2d<Real> *p_vect, const Array2d<Real> &mat) {
#ifdef DEBUG
  if (mat.get_size()(0)!=mat.get_size()(1)) ERRORQUIT("Matrix not square in diagonalize_symmetric_matrix");
#endif
  int n=mat.get_size()(0);
  int do_vect=(p_vect!=NULL);
  Array2d<Real> a(mat);
  Array<Real> d,e;
  to_tridiag(&a,&d,&e,do_vect);
  tridiag_to_diag(&d,&e,&a,do_vect);
  (*p_lambda)=d;
  if (do_vect) {
    (*p_vect)=a;
  }
}

void diagonalize_symmetric_matrix(Array<Real> *p_lambda, Array2d<Complex> *p_vect, const Array2d<Complex> &mat) {
#ifdef DEBUG
  if (mat.get_size()(0)!=mat.get_size()(1)) ERRORQUIT("Matrix not square in diagonalize_symmetric_matrix");
#endif
  int n=mat.get_size()(0);
  int do_vect=(p_vect!=NULL);
  Array2d<Real> a(2*n,2*n);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      a(i,j)=real(mat(i,j));
      a(n+i,n+j)=real(mat(i,j));
    }
  }
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      a(n+i,j)=imag(mat(i,j));
      a(i,n+j)=imag(mat(j,i));
    }
  }
  Array<Real> d,e;
  to_tridiag(&a,&d,&e,do_vect);
  tridiag_to_diag(&d,&e,&a,do_vect);
  p_lambda->resize(n);
  if (do_vect) {
    p_vect->resize(n,n);
  }
  Array<int> to_omit(2*n);
  zero_array(&to_omit);
  int j=0;
  for (int i=0; i<n; i++) {
    while (to_omit(j)==1) j++;
    (*p_lambda)(i)=d(j);
    if (do_vect) {
      for (int k=0; k<n; k++) {
	(*p_vect)(k,i)=Complex(a(k,j),a(k+n,j));
      }
    }
    to_omit(j)=1;
    int k;
    int best_k=-1;
    Real best_d=MAXFLOAT;
    for (k=j+1; k<2*n; k++) {
      if (to_omit(k)!=1 && fabs(d(k)-d(j))<best_d) {
	best_d=fabs(d(k)-d(j));
	best_k=k;
      }
    }
    to_omit(best_k)=1;
  }
}
