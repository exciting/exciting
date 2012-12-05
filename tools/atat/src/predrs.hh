#include "refine.h"
#include "keci.h"
#include "linalg.h"

class RSPredictorBase: public EnergyPredictor {
  static Structure lattice;                            // the parent lattice;
  virtual KSpaceECI *get_kspace_eci_object(void) {return NULL;}    // return pointer object returning kspace eci;

  Array<rVector3d> kpts; // k-point list;
  Array<Array<Complex> > str_fact; // structure factors;
  Real conc; // concentration between 0 and 1;
  int vol; // in multiple of parent lattice unit cell volume;
  int done; // flag: done with finding k-points;
  
 public:
  // Constructor: don't touch;
  RSPredictorBase(void);

  // init static members of class (those shared among all objects of type RSPredictor;
  void static_init(const ClusterExpansion &ce);

  // calculate constituent strain energy for structure str;
  Real get_energy(const StructureInfo &str);
};

template<class T>
class RSPredictor: public RSPredictorBase {
  static T kspace_eci_object;
  KSpaceECI *get_kspace_eci_object(void) {return &kspace_eci_object;}
 public:
  RSPredictor(void): RSPredictorBase() {}
};
