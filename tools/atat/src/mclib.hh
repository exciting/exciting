#ifndef _MCLIB_H_
#define _MCLIB_H_

#include <fstream.h>
#include <time.h>
#include "clus_str.h"
#include "linalg.h"
#include "keci.h"

#define SPIN_TYPE signed char

class MonteCarlo {
 protected:
  Structure lattice;
  iVector3d supercell;
  iVector3d margin;
  iVector3d total_box;
  SPIN_TYPE *spin;
  SPIN_TYPE *spin_changed;
  Real rspin_size;
  int site_in_cell;
  int total_clusters;
  int *nb_clusters;
  int *nb_x_clusters;
  int **cluster_size;
  Real *rcluster_mult_per_atom;
  int ***site_offset;
  Real **eci;
  int **which_cluster;
  Real E_ref;
  int which_is_empty;
  Real cur_energy;
  Array<Real> cur_rho;
  Real *pcur_rho;
  Real cur_conc;
  Real cur_disorder_param;

  Real T;
  Real mu;
 public:
  void calc_from_scratch(void);
  void spin_flip(void);
  void spin_double_flip(void);
  MonteCarlo(const Structure &_lattice, const iVector3d &_supercell, const SpaceGroup &space_group,
             const LinkedList<Cluster> &cluster_list);
  ~MonteCarlo(void);
  void set_eci(const Array<Real> &eci);
  void init_random(Real concentration=0.);
  void init_structure(const Structure &str);
  void set_concentration(Real concentration);
  void init_run(Real _T, Real _mu);
  void run(int mc_passes, int mode=1); // 1: grand-canonical, 2: canonical;
  void view(const Array<Arrayint> &labellookup, const Array<AutoString> &atom_label, ofstream &file, const rMatrix3d &axes);

  const Array<Real> & get_cur_corr(void) const {
    return cur_rho;
  }
  Real get_cur_energy(void) {
    return E_ref+cur_energy+extension_get_energy();
  }
  Real get_cur_concentration(void) const {
    return cur_conc;
  }
  Real get_cur_disorder_param(void) const {
    return cur_disorder_param;
  }
  int get_total_clusters(void) const {
    return total_clusters;
  }
  void get_cluster_mult(Array<Real> *mult) const {
    mult->resize(total_clusters);
    for (int i=0; i<total_clusters; i++) {(*mult)(i)=rcluster_mult_per_atom[i];}
  }
  Real get_T(void) const {return T;}
  Real get_mu(void) const {return mu;}

 protected:
  virtual void extension_calc_from_scratch(void) {}
  virtual void extension_update_spin_flip(int *cell, int incell, int newspin, Real d_recip_energy) {}
  virtual void extension_undo_spin_flip(void) {}
  virtual Real extension_spin_flip_energy(int *cell, int incell, int newspin) {return 0.;}
  virtual Real extension_get_energy(void) {return 0.;}
};

class MCOutputData {
 public:
  Real T,mu;
  Real E,x,lro;
  Real heatcap,suscept,var_lro;
  Array<Real> corr;
  Array<Real> mult;
  int n_equil,n_step;
  MCOutputData(void): corr(),mult() {}
};

class Equilibrator;
class Accumulator;

const int nb_value_accum=3;

void run_mc_until_converged(MonteCarlo *pmc, int mode, Equilibrator *eq, Real prec,  int which_elem=1 );
void run_mc(MonteCarlo *pmc, int mode, int n_step, Accumulator *eq);
void run_mc(MCOutputData *pmcdata, MonteCarlo *pmc, int mode, int n_step, Real prec, int which_elem=1 );

extern Array<Real> mclibdummyarray;

class GenericAccumulator {
 public:
  virtual int new_data(const Array<Real> &data) {return 0;}
  virtual const Array<Real> & get_mean(void) {return mclibdummyarray;}
  virtual const Array<Real> & get_var(void) {return mclibdummyarray;}
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
  Equilibrator(Real _prec, int data_size=1, int _which_elem=0, int init_granularity=500, int nb_bin=16);
  void init(Real _prec, int data_size=1, int _which_elem=0, int init_granularity=500, int nb_bin=16);
  int new_data(const Array<Real> &data);
  const Array<Real> & get_mean(void) {return good_val;}
  const Array<Real> & get_var(void) {return good_val2;}
  void get_step(int *pequil,int *pstep) {if (pequil) *pequil=equil; if (pstep) *pstep=step;}
};

class Accumulator: public GenericAccumulator {
  int cur_n;
  Array<Real> cur_sum;
  Array<Real> cur_sum2;
  Array<Real> good_val;
  Array<Real> good_val2;
  void accum(void);
public:
  Accumulator(void): cur_sum(), cur_sum2(), good_val() {cur_n=0;}
  int new_data(const Array<Real> &data);
  const Array<Real> & get_mean(void) {accum(); return good_val; }
  const Array<Real> & get_var(void) {accum(); return good_val2; }
};

class FlippedSpin {
 public:
  FlippedSpin(int *_cell, int _incell, SPIN_TYPE _newspin, Real _d_recip_energy) {
	for (int i=0; i<3; i++) {cell[i]=_cell[i];}
	incell=_incell; newspin=_newspin;
	d_recip_energy=_d_recip_energy;
  }
  int cell[3];
  int incell;
  SPIN_TYPE newspin;
  Real d_recip_energy;
};

inline int unrollsym(int i, int j) {
  int mx=MAX(i,j);
  return mx*(mx+1)/2+MIN(i,j);
}

class KSpaceMonteCarlo: public MonteCarlo {
 protected:
  int nsite;
  int size;
  Real rsize;
  Complex **ft_spin;
  Complex **ft_eci;
  Complex **dir_eci;
  Complex **convol;
  LinkedList<FlippedSpin> flipped_spins;
  KSpaceECI *p_kspace_eci;
  Real cur_recip_E;
  clock_t fft_time;
  clock_t flip_time;
  Real ref_x;
  Real threshold_dx;

 public:  
  KSpaceMonteCarlo(const Structure &_lattice, const iVector3d &_supercell, const SpaceGroup &space_group,
		   const LinkedList<Cluster> &cluster_list, KSpaceECI *_p_kspace_eci);
  ~KSpaceMonteCarlo();
 protected:
  void set_k_space_eci(void);
  void extension_calc_from_scratch(void);
  void extension_update_spin_flip(int *cell, int incell, int newspin, Real d_recip_energy);
  void extension_undo_spin_flip(void);
  Real extension_spin_flip_energy(int *cell, int incell, int newspin);
  Real extension_get_energy(void) {return cur_recip_E;}
};

#endif
