#include "stringo.h"
#include "arraylist.h"
#include "linalg.h"

void make_eci_list(LinkedList<Real> *pecilist, const PolyInterpolatorBig<Array<Real> > &teci, Real T);

void read_t_dep_eci_file(PolyInterpolatorBig<Array<Real> > *pteci, PolyInterpolatorBig<Array<Real> > *peeci, int ncluster, Real kboltzman);

void fix_energy(Real *pE, const Array<Real> &mult, const Array<Real> &corr, const PolyInterpolatorBig<Array<Real> > &teci, Real T);
