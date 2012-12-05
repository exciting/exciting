#ifndef __FIXAGG_H__
#define __FIXAGG_H__

template<class T>
class Aggregate {
public:
  T x;
  Aggregate(void) {}
  Aggregate(const T &r) {x=r;}
  operator T & () {return x;}
};

#endif
