#include <fstream.h>
#include "mteci.h"

void make_eci_list(LinkedList<Real> *pecilist, const PolyInterpolatorBig<Array<Real> > &teci, Real T) {
  pecilist->delete_all();
  Array<Real> eciarray;
  teci.interpol(&eciarray,T);
  (*pecilist) << eciarray;
}

void read_t_dep_eci_file(PolyInterpolatorBig<Array<Real> > *pteci, PolyInterpolatorBig<Array<Real> > *pteeci, int ncluster, Real kboltzman)  {
  ifstream tecifile("teci.out");
  if (tecifile) {
    Real Tmax;
    int nT;
    tecifile >> Tmax >> nT;
    Tmax*=kboltzman;
    Array<Array<Real> > ecis(nT);
    for (int it=0; it<nT; it++) {
      ecis(it).resize(ncluster);
      for (int c=0; c<ncluster; c++) {
	tecifile >> ecis(it)(c);
      }
    }
    pteci->init(0.,Tmax,ecis);
    Array<Array<Real> > eecis(nT-1);
    eecis(0)=ecis(0);
    Real dT=Tmax/(Real)(nT-1);
    for (int it=1; it<nT-1; it++) {
      Real T=dT*(Real)it;
      diff(&eecis(it),ecis(it+1),ecis(it-1));
      product(&eecis(it),eecis(it),-T/(2.*dT));
    }
    pteeci->init(0.,Tmax,eecis);
  } else {
    ifstream ecifile("eci.out");
    if (!ecifile) ERRORQUIT("Unable to open eci.out");
    Array<Array<Real> > ecis(1);
    Array<Array<Real> > eecis(1);
    ecis(0).resize(ncluster);
    eecis(0).resize(ncluster);
    for (int c=0; c<ecis(0).get_size(); c++) {
      ecifile >> ecis(0)(c);
      eecis(0)(c)=0.;
    }
    pteci->init(0.,1.,ecis);
    pteeci->init(0.,1.,eecis);
  }
}

void fix_energy(Real *pE, const Array<Real> &mult, const Array<Real> &corr, const PolyInterpolatorBig<Array<Real> > &teci, Real T) {
  Array<Real> eci1,eci0,eeci,tmp;
  Real dT=teci.max()*1e-3;
  teci.interpol(&eci1,T+dT);
  teci.interpol(&eci0,T-dT);
  diff(&eeci,eci1,eci0);
  product(&eeci,eeci,-T/(2.*dT));
  product_diag(&tmp,mult,corr);
  (*pE)+=inner_product(eeci,tmp);
}
