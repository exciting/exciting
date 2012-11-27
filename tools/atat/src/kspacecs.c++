#include <fstream.h>
#include "xtalutil.h"
#include "keci.h"
#include "linalg.h"
#include "plugin.h"

// Kubic Harmonics
// l=0,1,2,3 gives K_0,K_4,K_6,K_8;
Real kubic_harm(const rVector3d &k, int l) {
  Real r2=norm2(k);
  Real A=sqrt(1./(4.*M_PI));
  switch (l) {
  case 0: return(A);
  case 1: return(A*sqrt(21./4.)*(1.-5.*(ipow(k(0),2)*ipow(k(1),2) + ipow(k(0),2)*ipow(k(2),2) + ipow(k(1),2)*ipow(k(2),2)) / ipow(r2,2)));
  case 2: return(A*sqrt(13./2.)*(1./4.)*(7*(ipow(k(0),6)+ipow(k(1),6)+ipow(k(2),6)+30*ipow(k(0),2)*ipow(k(1),2)*ipow(k(2),2))/ipow(r2,3) - 5));
  case 3: return(A*sqrt(561.)*(1./8.)*(ipow(k(0),8) + ipow(k(1),8) + ipow(k(2),8) - 14*(ipow(k(0),6)*ipow(k(1),2) + ipow(k(0),6)*ipow(k(2),2) + ipow(k(1),6)*ipow(k(2),2) + ipow(k(0),2)*ipow(k(1),6) + ipow(k(0),2)*ipow(k(2),6) + ipow(k(1),2)*ipow(k(2),6)) + 35.*(ipow(k(0),4)*ipow(k(1),4) + ipow(k(0),4)*ipow(k(2),4) + ipow(k(1),4)*ipow(k(2),4))) / ipow(r2,4));
  default: ERRORQUIT("Kubic harmonics l > 8 not implemented yet.");
  }
  return 0.;
}

#ifndef NO_OBJECTS_IN_KSPACE_CS

#ifdef STATICCSCOEF
#define STATICCSCOEFP static
#else
#define STATICCSCOEFP
#endif

class CS_KSpaceECI: public KSpaceECI {
  STATICCSCOEFP Array<LinearInterpolator<Real> > kubic_coef;  // kubic_coef(l)(conc), l=0,1,2,3,... conc in [0,1];
  STATICCSCOEFP Real ikc;                                     // inverse k-vector cutoff for attenuation;
  STATICCSCOEFP int nb_atom_in_cell;                          // nb atom per unit cell of lattice;

 public:
  // Constructor: don't touch;
#ifdef STATICCSCOEF
  CS_KSpaceECI(void): KSpaceECI() {}
#else
  CS_KSpaceECI(void): KSpaceECI(), kubic_coef() {}
#endif
  // read in parameters from the file cs.in;
  void static_init(const Structure &lattice) {
    nb_atom_in_cell=lattice.atom_pos.get_size();
    // open file and error checking;
    ifstream file("cs.in");
    if (!file) {
      cerr << "Unable to open cs.in. Will ignore constituent strain." << endl;
      kubic_coef.resize(0);
    }
    else {
      // if ok, read file in;
      file >> ikc;
      int n_term=0;
      int n_conc=0;
      file >> n_term;
      file >> n_conc;
      kubic_coef.resize(n_term);
      for (int i=0; i<kubic_coef.get_size(); i++) {
	// read coef for kubic harmonic i, at all concentration;
	Array<Real> tmp(n_conc);
	for (int j=0; j<n_conc; j++) {
	  file >> tmp(j);
	}
	// initialize and LinearInterpolator object, see linalg.h;
	// (to interpolate between concentration grid points);
	kubic_coef(i).init(0.,1.,tmp);
      }
    }
  }

  // calculate K space ECI associated with constituent strain energy for structure str;
  void get_k_space_eci(Array2d<Complex> *p_keci, const rVector3d &k, Real x) {
    p_keci->resize(iVector2d(nb_atom_in_cell,nb_atom_in_cell));
    Real e=0.;
    if (norm2(k)>sqr(zero_tolerance)) {
      for (int l=0; l<kubic_coef.get_size(); l++) { // loop over kubic harmonics;
	e+=kubic_coef(l)(x) * kubic_harm(k,l) * exp(-norm2(k)*sqr(ikc));
      }
    }
    for (int i=0; i<nb_atom_in_cell; i++) {
      for (int j=0; j<nb_atom_in_cell; j++) {
	(*p_keci)(i,j)=e;
      }
    }
  }
};

SpecificPlugIn<KSpaceECI,CS_KSpaceECI> CS_KSpaceECIRegister("cs"); // register plug-in under name "cs";

#ifdef STATICCSCOEF
// static storage space of parameters;
Array<LinearInterpolator<Real> > CS_KSpaceECI::kubic_coef;
Real CS_KSpaceECI::ikc;
int CS_KSpaceECI::nb_atom_in_cell;
#endif

#endif
