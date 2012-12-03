#ifndef _EQUIL_H_
#define _EQUIL_H_

#include "arraylist.h"
#include "linalg.h"

extern Array<Real> equildummyarray;

class GenericAccumulator {
 public:
  virtual int new_data(const Array<Real> &data) {return 0;}
  virtual const Array<Real> & get_mean(void) {return equildummyarray;}
  virtual const Array<Real> & get_var(void) {return equildummyarray;}
  virtual void get_step(int *pequil,int *pstep) { }
};

class Equilibrator: public GenericAccumulator {
  int granularity;
  int cur_bin;
  int cur_cnt;
  Array<Real> cur_sum;
  Array<Real> cur_sum2;
  Array<ArrayReal> bin_sum;
  Array<ArrayReal> bin_sum2;
  Array<Real> buf_corr;
  int cur_corr_len;
  Real prec;
  int which_elem;
  Array<Real> good_val;
  Array<Real> good_val2;
  Array<Real> buf_val;
  int equil,step;
 public:
  Equilibrator(void);
  Equilibrator(Real _prec, int _which_elem=0, int init_granularity=200, int nb_bin=16);
  void init(Real _prec, int _which_elem=0, int init_granularity=200, int nb_bin=16);
  int new_data(const Array<Real> &data);
  const Array<Real> & get_mean(void) {return good_val;}
  const Array<Real> & get_var(void) {return good_val2;}
  void get_step(int *pequil,int *pstep) {if (pequil) *pequil=equil; if (pstep) *pstep=step;}
};

class Accumulator: public GenericAccumulator {
  int max_n;
  int cur_n;
  Array<Real> cur_sum;
  Array<Real> cur_sum2;
  Array<Real> good_val;
  Array<Real> good_val2;
  void accum(void);
public:
  Accumulator(int _max_n=0): cur_sum(), cur_sum2(), good_val() {max_n=_max_n; cur_n=0;}
  int new_data(const Array<Real> &data);
  const Array<Real> & get_mean(void) {accum(); return good_val; }
  const Array<Real> & get_var(void) {accum(); return good_val2; }
  void get_step(int *pequil,int *pstep) {if (pstep) *pstep=cur_n;}
};

GenericAccumulator *create_accum(int n_step, Real prec, int which_elem=0);
#endif
