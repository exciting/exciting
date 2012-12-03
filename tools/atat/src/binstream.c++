#include <iostream.h>
#include <strstream.h>
#include "misc.h"
#include "vectmac.h"
#include "array.h"

#define MAKE_GENERIC_BIN_OSTREAM(T) ostream & bin_ostream(ostream &file, const T &x) { return file.write((char *)&x,sizeof(T)); }

MAKE_GENERIC_BIN_OSTREAM(int)
MAKE_GENERIC_BIN_OSTREAM(Real)

template<class T, int D>
ostream & bin_ostream(ostream &file, FixedVector<T,D> &x) {
  for (int i=0; i<D; i++) {
    bin_ostream(file, (T &)(x(i)));
  }
  return file;
}

template<class T>
ostream & bin_ostream(ostream &file, const Array<T> &x) {
  bin_ostream(file,x.get_size());
  for (int i=0; i<x.get_size(); i++) {
    bin_ostream(file, (T &)(x(i)));
  }
  return file;
}

#define MAKE_GENERIC_BIN_ISTREAM(T) istream & bin_istream(istream &file, T &x) {  return file.read((char *)&x,sizeof(T)); }

MAKE_GENERIC_BIN_ISTREAM(int)
MAKE_GENERIC_BIN_ISTREAM(Real)

template<class T, int D>
istream & bin_istream(istream &file, FixedVector<T,D> &x) {
  for (int i=0; i<D; i++) {
    bin_istream(file, (T &)(x(i)));
  }
  return file;
}

template<class T>
istream & bin_istream(istream &file, Array<T> &x) {
  int n;
  bin_istream(file, n);
  x.resize(n);
  for (int i=0; i<n; i++) {
    bin_istream(file, (T &)(x(i)));
  }
  return file;
}

int main(int argc, char *argv[]) {
  Real x=0.1234;
  rVector3d y(1.1,2.0,3.0);
  Array<Array<Real> > z(2);
  z(0).resize(3);
  z(0)(0)=1.1; z(0)(1)=1.2; z(0)(2)=1.3;
  z(1).resize(2);
  z(1)(0)=2.1; z(1)(1)=2.2;
  ostrstream line;
  bin_ostream(line,x);
  bin_ostream(line,y);
  bin_ostream(line,z);
  strstream buf;
  buf.write(line.str(),line.tellp());
  Real xx;
  rVector3d yy;
  Array<Array<Real> > zz;
  bin_istream(buf,xx);
  bin_istream(buf,yy);
  bin_istream(buf,zz);
  cout << x << " " << xx << endl << y << " " << " " << yy << endl << z << " " << zz << endl;
  return 0;
}
