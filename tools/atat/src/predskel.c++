#include "refine.h"

// replace MyOwn by the name of your new predictor
// look for ?? and edit;

class MyOwnPredictor: public EnergyPredictor {
  static Real ??; 
 public:
  MyOwnPredictor(void): EnergyPredictor() {}
  void static_init(const ClusterExpansion &ce) {
    ??=??;
  }
  Real get_energy(const StructureInfo &str) {
    return (??);
  }
};

SpecificPlugIn<EnergyPredictor,MyOwnPredictor> MyOwnPredictorObject("MyOwn");

static Real MyOwnPredictor::??; 

