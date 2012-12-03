#include "stringo.h"
#include "arraylist.h"
#include "linalg.h"
//#include "mclib.h"
#include "drawpd.h"

void make_eci_list(LinkedList<Real> *pecilist, const PolyInterpolatorBig<Array<Real> > &teci, Real T);

void read_t_dep_eci_file(PolyInterpolatorBig<Array<Real> > *pteci, PolyInterpolatorBig<Array<Real> > *peeci, int ncluster, Real kboltzman);

//void fix_energy(MCOutputData *pmcdata, const PolyInterpolatorBig<Array<Real> > &teeci, Real T);

void fix_energy(Real *pE, const Array<Real> &mult, const Array<Real> &corr, const PolyInterpolatorBig<Array<Real> > &teci, Real T);

class TDependentECIPhase: public PhaseThermoData {
  AutoString labeltocreate;
  Structure lattice;
  SpaceGroup space_group;
  Structure str;
  LinkedList<ArrayrVector3d> cluster_list;
  PolyInterpolatorBig<Array<Real> > teci;
  PhaseThermoData *phase;
 public:
  TDependentECIPhase(char *_label,
		     const Structure &_lattice, const SpaceGroup &_space_group,
		     const Structure &_str,
		     const LinkedList<ArrayrVector3d> &_cluster_list,
		     const PolyInterpolatorBig<Array<Real> > &_teci);
  ~TDependentECIPhase(void) {delete phase;}
  void set_T_mu(Real T, Real mu);
  void write(ostream &file);
};
