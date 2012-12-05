#ifndef _DRAWPD_H_
#define _DRAWPD_H_

#include <fstream.h>
#include "arraylist.h"
#include "calccorr.h"

class PhaseThermoData {
protected:
  Real cur_phi;
  Real cur_x;
  PhaseThermoData *cur_phase;
public:
  PhaseThermoData(void) {}
  virtual void set_T_mu(Real T, Real mu) {}
  virtual void read(istream &file) {}
  virtual void write(ostream &file) {}
  virtual char *object_label(void) {return "";}
  virtual PhaseThermoData *new_object(void) {return this;}
  virtual void init(  const Structure &lattice, const SpaceGroup &space_group,
		      const Structure &str,
		      const LinkedList<ArrayrVector3d> &cluster_list,
		      const LinkedList<Real> &eci_list) {}
  Real phi(void) {return cur_phi;}
  Real x(void) {return cur_x;}
  PhaseThermoData* get_sub_phase(void) {return cur_phase;}
  Real phi(Real T, Real mu) {set_T_mu(T,mu); return cur_phi;}
  Real x(Real T, Real mu) {set_T_mu(T,mu); return cur_x;}
  PhaseThermoData* get_sub_phase(Real T, Real mu) {set_T_mu(T,mu); return cur_phase;}
};

class MF_ThermoData: public PhaseThermoData {
 protected:
  static Real spin_precision;
  Real volume;
  Real max_flip;
  Array<Real> ideal_spin;
  Array<ArrayArrayint> clusters;
  Array<ArrayReal> eci;
  Real E_ref;
 public:
  MF_ThermoData(void);
  MF_ThermoData(  const Structure &lattice, const SpaceGroup &space_group,
                  const Structure &str,
                  const LinkedList<ArrayrVector3d> &cluster_list,
		  const LinkedList<Real> &eci_list);
  void init(  const Structure &lattice, const SpaceGroup &space_group,
	      const Structure &str,
	      const LinkedList<ArrayrVector3d> &cluster_list,
	      const LinkedList<Real> &eci_list);
  void set_T_mu(Real T, Real mu);
  void read(istream &file);
  void write(ostream &file);
  char *object_label(void) {return "MF";}
  PhaseThermoData *new_object(void) {return new MF_ThermoData();}
};

class LTE_ThermoData: public PhaseThermoData {
 public:
  Real E0;
  Real x0;
  Real volume;
  Array<Real> dE;
  Array<Real> dx;
  Real cur_E;
 public:
  LTE_ThermoData(void);
  LTE_ThermoData( const Structure &lattice, const SpaceGroup &space_group,
                  const Structure &str,
                  const LinkedList<ArrayrVector3d> &cluster_list, const LinkedList<Real> &eci_list);
  void init(  const Structure &lattice, const SpaceGroup &space_group,
	      const Structure &str,
	      const LinkedList<ArrayrVector3d> &cluster_list,
	      const LinkedList<Real> &eci_list);
  void set_T_mu(Real T, Real mu);
  Real E(Real T, Real mu);
  Real E(void) {return cur_E;}
  void read(istream &file);
  void write(ostream &file);
  char *object_label(void) {return "LTE";}
  PhaseThermoData *new_object(void) {return new LTE_ThermoData();}
};

class HTE_ThermoData: public PhaseThermoData {
 protected:
  Array<Real> poly_x;
  static Real dx;
  // protected:
 public:
  Real F_of_x(Real T, Real x);
  Real mu_of_x(Real T, Real x);
  Real x_of_mu(Real T, Real mu);
 public:
  HTE_ThermoData(void);
  HTE_ThermoData(const Structure &lattice,
		 const SpaceGroup &space_group,
		 const Structure &str,
		 const LinkedList<ArrayrVector3d> &cluster_list,
		 const LinkedList<Real> &eci_list);
  void init(  const Structure &lattice, const SpaceGroup &space_group,
	      const Structure &str,
	      const LinkedList<ArrayrVector3d> &cluster_list,
	      const LinkedList<Real> &eci_list);
  void simple_init(int atom_per_unit_cell, const SpaceGroup &space_group,
	      const LinkedList<ArrayrVector3d> &cluster_list, const LinkedList<Real> &eci_list);
  void set_T_mu(Real T, Real mu);
  void read(istream &file);
  void write(ostream &file);
  char *object_label(void) {return "HTE";}
  PhaseThermoData *new_object(void) {return new HTE_ThermoData();}
};

class SampledThermoData: public PhaseThermoData {
  Real T0;
  Real dT;
  Array<Real> mu0;
  Real dmu;
  Array<ArrayReal> xs;
  Array<ArrayReal> phis;
 public:
  SampledThermoData(void);
  SampledThermoData(Real _T0, Real _dT, const Array<Real> &_mu0,
                     Real _dmu, const Array<ArrayReal> &_xs,
                     const Array<ArrayReal> &_phis);
  void set_T_mu(Real T, Real mu);
  void read(istream &file);
  void write(ostream &file);
  char *object_label(void) {return "TAB";}
  PhaseThermoData *new_object(void) {return new SampledThermoData();}
  const Array<Real> & get_mu_data(void) {return mu0;}
  const Array<ArrayReal> & get_x_data(void) {return xs;}
  const Array<ArrayReal> & get_phi_data(void) {return phis;}
};

class PhaseListT: public PhaseThermoData {
 protected:
  LinkedList<PhaseThermoData> list;
  LinkedList<Real> T_limit;
 protected:
  PhaseThermoData* find_phase(Real T);
 public:
  PhaseListT(void);
  PhaseThermoData *add_phase(PhaseThermoData *phase, Real Tmax=MAXFLOAT);
  void set_T_mu(Real T, Real mu);
  void read(istream &file);
  void write(ostream &file);
  char *object_label(void) {return "lstT";}
  PhaseThermoData *new_object(void) {return new PhaseListT();}
};

class PhaseList: public PhaseThermoData {
 protected:
  LinkedList<PhaseThermoData> list;
 protected:
  PhaseThermoData* find_phase(Real T, Real mu);
 public:
  PhaseList(void);
  PhaseThermoData *add_phase(PhaseThermoData *phase);
  void set_T_mu(Real T, Real mu);
  void read(istream &file);
  void write(ostream &file);
  char *object_label(void) {return "lst";}
  PhaseThermoData *new_object(void) {return new PhaseList();}
};

class MiscibilityGap: public PhaseThermoData {
 public:
  PhaseThermoData *lo[2];
  PhaseThermoData *hi;
 public:
  MiscibilityGap(void);
  MiscibilityGap(PhaseThermoData *_hi, PhaseThermoData *_lo0, PhaseThermoData *_lo1);
  ~MiscibilityGap();
  void set_T_mu(Real T, Real mu);
  void read(istream &file);
  void write(ostream &file);
  char *object_label(void) {return "MG";}
  PhaseThermoData *new_object(void) {return new MiscibilityGap();}
};

class PhaseObjectList {
 public:
  LinkedList<PhaseThermoData> list;
  PhaseObjectList(void);
};

extern PhaseObjectList phaseobjectlist;

PhaseThermoData *read_phase(istream &file);

class PhaseBoundary {
 public:
  Real x[2];
  PhaseThermoData *phase[2];
};

extern Real tiny_phi;

class TwoPhaseRegion {
 public:
  LinkedList<ArrayReal> x_list;
  LinkedList<Real> T_list;
  PhaseThermoData *phase[2];
};

void recurse_phase_boundaries(LinkedList<PhaseBoundary> *pb_list, const LinkedList<PhaseThermoData> &phases, Real T, Real mu0, Real mu1, Real dx);
void calc_phase_boundaries(LinkedList<PhaseBoundary> *pb_list, const LinkedList<PhaseThermoData> &phases,
                           Real T, Real dx, const Array<Real> &mu_hints);
int same_phases(const LinkedList<PhaseBoundary> &pb_list1, const LinkedList<PhaseBoundary> &pb_list2);
void add_one_T(LinkedList<TwoPhaseRegion> *two_phase_list, Real T, const LinkedList<PhaseBoundary> &pb_list);

void calc_phase_boundaries(LinkedList<TwoPhaseRegion> *two_phase_list, const LinkedList<PhaseThermoData> &phases,
                           Real dT, Real dx, const Array<Real> &mu_hints);
void write_phase_diagram(const LinkedList<TwoPhaseRegion> &two_phase_list, ostream &file);

Real E_from_phi(Real T, Real mu, PhaseThermoData *p);

#endif
