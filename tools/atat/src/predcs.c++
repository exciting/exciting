#include "predrs.h"

#define STATICCSCOEF static

#include "kspacecs.cc"

SpecificPlugIn<EnergyPredictor, RSPredictor<CS_KSpaceECI> > CSPredictorRegister("cs"); // register plug-in under name "cs";

// static storage space of parameters;
CS_KSpaceECI RSPredictor<CS_KSpaceECI>::kspace_eci_object=CS_KSpaceECI();     // object returning kspace eci;
