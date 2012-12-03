#ifndef __BINSTREAM_H__
#define __BINSTREAM_H__

#include <iostream.h>
#include <strstream.h>
#include "misc.h"
#include "vectmac.h"
#include "arraylist.h"

#define MAKE_GENERIC_BIN_OSTREAM(T) inline ostream & bin_ostream(ostream &file, const T &x) { return file.write((char *)&x,sizeof(T)); }

MAKE_GENERIC_BIN_OSTREAM(char)
MAKE_GENERIC_BIN_OSTREAM(int)
MAKE_GENERIC_BIN_OSTREAM(Real)

template<class T, int D>
inline ostream & bin_ostream(ostream &file, const FixedVector<T,D> &x) {
  for (int i=0; i<D; i++) {
    bin_ostream(file, (T &)(x(i)));
  }
  return file;
}

template<class T, int D>
inline ostream & bin_ostream(ostream &file, const FixedMatrix<T,D> &x) {
  for (int i=0; i<D; i++) {
      for (int j=0; j<D; j++) {
	  bin_ostream(file, (T &)(x(i,j)));
      }
  }
  return file;
}

template<class T>
inline ostream & bin_ostream(ostream &file, const Array<T> &x) {
  bin_ostream(file,x.get_size());
  for (int i=0; i<x.get_size(); i++) {
    bin_ostream(file, (T &)(x(i)));
  }
  return file;
}

template<class T>
inline ostream & bin_ostream(ostream &file, const LinkedList<T> &x) {
  bin_ostream(file,x.get_size());
  LinkedListIterator<T> i(x);
  for (; i; i++) {
    bin_ostream(file, (T &)(*i));
  }
  return file;
}

#define MAKE_GENERIC_BIN_ISTREAM(T) inline istream & bin_istream(istream &file, T &x) {  return file.read((char *)&x,sizeof(T)); }

MAKE_GENERIC_BIN_ISTREAM(char)
MAKE_GENERIC_BIN_ISTREAM(int)
MAKE_GENERIC_BIN_ISTREAM(Real)

template<class T, int D>
inline istream & bin_istream(istream &file, FixedVector<T,D> &x) {
  for (int i=0; i<D; i++) {
    bin_istream(file, (T &)(x(i)));
  }
  return file;
}

template<class T, int D>
inline istream & bin_istream(istream &file, FixedMatrix<T,D> &x) {
  for (int i=0; i<D; i++) {
      for (int j=0; j<D; j++) {
	  bin_istream(file, (T &)(x(i,j)));
      }
  }
  return file;
}

template<class T>
inline istream & bin_istream(istream &file, Array<T> &x) {
  int n;
  bin_istream(file, n);
  x.resize(n);
  for (int i=0; i<n; i++) {
    bin_istream(file, (T &)(x(i)));
  }
  return file;
}

template<class T>
inline istream & bin_istream(istream &file, LinkedList<T> &x) {
  x.delete_all();
  int n;
  bin_istream(file, n);
  for (int i=0; i<n; i++) {
    T *obj=new T;
    bin_istream(file, (T &)(*obj));
    x << obj;
  }
  return file;
}

#endif
