#include "calccorr.h"
#include "linalg.h"

class GCThermoData {
 public:
  virtual void set_T_mu(Real _T, const Array<Real> &_mu) {}
  virtual Real get_phi(void) {return 0.;}
  virtual Real set_get_phi(Real _T, const Array<Real> &_mu) {
    set_T_mu(_T,_mu);
    return get_phi();
  }
  virtual Real get_lro(void) {return 1.;}
};


class iMultiCluster {
 public:
  int which_eci;
  Array<int> site;
  Array<int> site_type;
  Array<int> func;
  iMultiCluster(void): site(), site_type(), func() {}
};

class CalcMeanField: public GCThermoData {
  Array<Array<Real> > ref_prob;
  Array<Array<iMultiCluster> > clusters;
  const Array<Array<Array<Real> > > &corrfunc;
  Array<Real> eci;
  Real phi;
  Real lro;
  int which_is_empty;
  Real empty_mult;
  int do_LTE;
  Real calc_point_corr(const Array<Real> &siteprob, int site_type, int func);
 public:
  CalcMeanField(const Structure &_lattice, const Array<Array<int> > &_site_type_list, const Structure &_str, 
		const SpaceGroup &space_group, const LinkedList<MultiCluster> &cluster_list, 
		const Array<Array<Array<Real> > > &_corrfunc);
  void set_eci(const Array<Real> &_eci) {eci=_eci;}
  void set_ref_prob(const Array<Array<Real> > &_ref_prob) {ref_prob=_ref_prob;}
  void set_LTE(int _do_LTE) {do_LTE=_do_LTE;}
  void set_T_mu(Real _T, const Array<Real> &_mu);
  Real get_phi(void) {return phi;}
  Real get_lro(void) {return lro;}
};

void calc_conc_from_phi(Array<Real> *pconc, GCThermoData *ptd, Real T, const Array<Real> &mu, Real dmu);
Real calc_E_from_phi(GCThermoData *ptd, Real T, const Array<Real> &mu, Real dT);
Real calc_E_from_phi(CalcMeanField *ptd, const PolyInterpolatorBig<Array<Real> > &teci, Real T, const Array<Real> &mu, Real dT);
