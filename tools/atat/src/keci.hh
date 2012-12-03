#ifndef __KECI_H__
#define __KECI_H__

#include <complex.h>
#include "xtalutil.h"
#include "linklist.h"

class KSpaceECI {
 public:
  virtual void static_init(const Structure &lattice) {}
  virtual void get_k_space_eci(Array2d<Complex> *p_keci, const rVector3d &k, Real x) {}
};

class MultiKSpaceECI: public KSpaceECI, public LinkedList<KSpaceECI> {
  KSpaceECI *theonlyone;
 public:
  MultiKSpaceECI(void): KSpaceECI(), LinkedList<KSpaceECI>() {theonlyone=NULL;}
  void get_k_space_eci(Array2d<Complex> *pft_eci, const rVector3d &k, Real x);
  void check_only_one(void);
};

#endif
