#include "integer.h"

Real zero_tolerance=1e-3;

int least_common_multiple(int a, int b) {
  int test=2;
  int divider=1;
  while (test<=a && test<=b) {
    if ((a/test)*test == a && (b/test)*test == b) divider=test;
    test=test+1;
  }
  return a*b/divider;
}

Real sign(Real x) {
  return (x<0 ? -1. : (x>0 ? 1. : 0.));
}

int integer_ratio(int *pp, int *pq, Real x, Real epsilon) {
  int &p=*pp;
  int &q=*pq;
  int maxit=20;
  int p1,q1,p2,q2,n,g;
  Real r;

  n=(int)rint(x);
  r=x-n;
  g=(int)sign(r);
  p=n;
  p1=p;
  q1=1;
  q=1;
  p2=1;
  q2=0;

  int counter=0;
  while (fabs((Real)p/(Real)q-x) > epsilon) {
    n=(int)rint(g/r);
    r=g/r-n;
    p=n*p1+g*p2;
    q=n*q1+g*q2;
    g=
    p2=p1;
    q2=q1;
    p1=p;
    q1=q;
    counter++;
    if (counter>maxit) return 0;
  }
  return 1;
}

int factorial(int n) {
  int f=1;
  for (int i=2; i<=n; i++) {f*=i;}
  return f;
}
