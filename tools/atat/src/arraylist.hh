#ifndef _ARRAYLIST_H_
#define _ARRAYLIST_H_

#include "array.h"
#include "linklist.h"

template<class T,class S>
inline LinkedList<T>& operator <<(LinkedList<T>& list, const Array<S> &a) {
  for (int i=0; i<a.get_size(); i++) {
    list << new T(a(i));
  }
  return list;
}

template<class T,class S>
void LinkedList_to_Array(Array<T> *a, const LinkedList<S> &l) {
  if (a) {
    int size=l.get_size();
    a->resize(size);
    LinkedListIterator<S> it(l);
    for (int i=0; i<size; i++) {
      (*a)(i)=(T)(*it);
      it++;
    }
  }
}

template<class T,class S>
void ArrayArray_to_Array2d(Array2d<T> *pmat, const Array<Array<S> > &a) {
  if (a.get_size()==0) {
    pmat->resize(0,0);
  }
  else {
    int h=a(0).get_size();
    pmat->resize(a.get_size(),h);
    for (int i=0; i<a.get_size(); i++) {
      for (int j=0; j<h; j++) {
        (*pmat)(i,j)=(T)(a(i)(j));
      }
    }
  }
}

template<class T,class S>
void LinkedListArray_to_Array2d(Array2d<T> *a, const LinkedList<Array<S> > &l) {
  if (a) {
    LinkedListIterator<S> it(l);
    iVector2d size(l.get_size(),it->get_size());
    a->resize(size);
    for (int i=0; i<size(0); i++) {
      for (int j=0; j<size(1); j++) {
	(*a)(i,j)=(T)((*it)(j));
      }
      it++;
    }
  }
}

template<class T>
class TrivEqual {
 public:
  int operator () (const T& a, const T& b) const {
    return (a==b);
  }
};

template<class T, class C>
inline int add_unique(LinkedList<T> *l, T *pobject, const C &isequal) {
  LinkedListIterator<T> i(*l);
  for ( ; i; i++) {
    if (isequal(*i,*pobject)) break;
  }
  if (!i) {
    (*l) << pobject;
    return 1;
  }
  else {
    return 0;
  }
}

template<class T, class C>
inline int add_unique(LinkedList<T> *l, const T &object, const C &isequal) {
  T *pobject=new T(object);
  if (add_unique(l,pobject,isequal)) {
    return 1;
  }
  else {
    delete pobject;
    return 0;
  }
}

template<class T, class C>
inline void add_sorted(LinkedList<T> *l, T *pobject, const C &less_equal) {
  LinkedListIterator<T> i(*l);
  for ( ; i; i++) {
    if (!less_equal(*i,pobject)) break;
  }
  l->add(pobject,i);
}

template<class T, class C>
inline void add_sorted(LinkedList<T> *l, const T &object, const C &comparator) {
  add_unique_sorted(l,new T(object),comparator);
}

template<class T>
class TrivLessThan {
 public:
  int operator () (const T& a, const T& b) const {
    return (a<b);
  }
};

template<class T, class C>
inline void sort_array(Array<T> *a, const C &comparator) {
  for (int i=1; i<a.get_size(); i++) {
    for (int j=0; j<a.get_size()-i; j++) {
      if (comparator(a(j+1),a(j))) {
	swap(&(a(j+1)),&(a(j)));
      }
    }
  }
}

template<class T>
class MultiDimIterator {
  T current;
  T max;
  T min;
  int valid;
 public:
  MultiDimIterator(): current(), min(), max() {}
  ~MultiDimIterator() {}
  MultiDimIterator(const T &_min, const T &_max): current(), min(), max() {
    init(_min,_max);
  }
  MultiDimIterator(const T &_max): current(), min(), max() {
    init(_max);
  }
  void init(const T &_min, const T &_max) {
    min=_min;
    max=_max;
    valid=1;
    current=_min;
  }
  void init(const T &_max) {
    min=_max;
    max=_max;
    for (int i=0; i<min.get_size(); i++) {
      min(i)=0;
      max(i)=max(i)-1;
    }
    valid=1;
    current=min;
  }
  operator const T& (void) {return current;}
  //  operator void * () {return (void *)valid;}
  operator void * () {return ((int *)NULL+valid);}
  const T& operator++(int);
  void bump_up(int level);  
};

template<class T>
inline const T& MultiDimIterator<T>::operator++(int) {
  int i=0;
  while (i<max.get_size()) {
    if (current(i)<max(i)) break;
    current(i)=min(i);
    i++;
  }
  if (i==max.get_size()) {
    valid=0;
  }
  else {
    current(i)++;
  }
  return current;
}

template<class T>
inline void MultiDimIterator<T>::bump_up(int level) {
  current(level)=min(level);
  int i=level+1;
  while (i<max.get_size() && current(i)>=max(i)) {
    current(i)=min(i);
    i++;
  }
  if (i==max.get_size()) {
    valid=0;
  }
  else {
    current(i)++;
  }
}

class MultipletIterator {
 public:
  enum Style {all,nopermute,norepeat};
 private:
  Array<int> current;
  int min;
  int max;
  int valid;
  Style style;
 public:
  MultipletIterator(): current() {}
  MultipletIterator(int dim, int _min, int _max, Style _style=all): current() {
    init(dim,_min,_max,_style);
  }
  void init(int dim, int _min, int _max, Style _style=all) {
    min=_min;
    max=_max;
    valid=1;
    style=_style;
    current.resize(dim);
    for (int i=0; i<current.get_size(); i++) {
      switch (style) {
        case norepeat:
          current(i)=min+i;
          break;
        default:
          current(i)=min;
          break;
      }
    }
  }
  operator const Array<int>& (void) {return current;}
  operator int () {return valid;}
  const Array<int>& operator++(int);
};

inline const Array<int>&  MultipletIterator::operator++(int) {
  int i;
  switch (style) {
    case MultipletIterator::all:
      i=0;
      while (current(i)>=max && i<current.get_size()) {
        current(i)=min;
        i++;
      }
      if (i==current.get_size()) {valid=0;}
      current(i)++;
      break;
    case MultipletIterator::nopermute:
      i=0;
      while (current(i)==max && i<current.get_size()) i++;
      if (i==current.get_size()) {valid=0; return current;}
      current(i)++;
      { for (int j=i-1; j>=0; j--) {current(j)=current(i);} }
      break;
    case MultipletIterator::norepeat:
      current(i)=min;
      i=0;
      while (current(i)==max && i<current.get_size()-1) i++;
      current(i)++;
      for (int j=i-1; j>=0; j--) {
        current(j)=current(i)+(i-j);
        if (current(j)>max) {valid=0; return current;}
      }
      break;
  }
  return current;
}

template<class T>
void extract_columns(Array2d<T> *pa, const Array2d<T> &b, const LinkedList<int> &cols_list, int nb_cols) {
  Array<int> cols;
  LinkedList_to_Array(&cols,cols_list);
  extract_columns(pa,b,cols,nb_cols);
}

#endif
