//#include <strstream.h>
#include "linalg.h"
#include "chull.h"

int main(void) {
  Array<Array<Real> > x(3);
  x(0).resize(3);
  x(0)(0)=0;  x(0)(1)=0;  x(0)(2)=0;
  x(1).resize(3);
  x(1)(0)=0.8;  x(1)(1)=0.8;  x(1)(2)=1;
  x(2).resize(3);
  x(2)(0)=0.4;  x(2)(1)=0.4;  x(2)(2)=3;
  Array<Real> ground(3);
  ground(0)=0; ground(1)=0; ground(2)=-1.;
  LinkedList<PolytopeFace> hull;
  calc_convex_hull_degenerate(&hull, x,ground);
  exit(1);
  {
    Array2d<Real> a(2,3);
    Array2d<Real> ia,tmp;
    a(0,0)=1.; a(0,1)=2.; a(0,2)=3.;
    a(1,0)=1.; a(1,1)=4.; a(1,2)=3.;
    invert_matrix_wide(&ia,a);
    cout << ia;
    product(&tmp,a,ia);
    cout << tmp;
  }
  {
    Array2d<Real> a(3,2);
    Array2d<Real> ia,tmp;
    a(0,0)=1.; a(0,1)=2.;
    a(1,0)=1.; a(1,1)=3.;
    a(2,0)=3.; a(2,1)=4.;
    invert_matrix_tall(&ia,a);
    cout << ia;
    product(&tmp,ia,a);
    cout << tmp;
  }
  /*
  char *str="1.234";
  istrstream buf(str);
  double x;
  buf >> x;
  cout << x << endl;
  */
}

