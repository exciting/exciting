#include "lapw.h"
#include "lapw_band.h"

void lapw_density(bloch_states_k *ks, mdarray<double,6>& densmt, mdarray<double,3>& densir);

extern "C" void FORTRAN(lapw_execute)(int32_t *ikloc_, complex16 *apwalm_, double *evalsv_, 
                                      double *occsv_, double *densmt_, double *densir_, int32_t *flags_)
{
    timer t("lapw_execute");
    
    int flags = *flags_;
    unsigned int ikloc = *ikloc_ - 1;
    
    bloch_states_k *ks = lapw_runtime.bloch_states[ikloc];
    
    ks->copy_apwalm(apwalm_);
    
    if (flags & 1) 
    {
        lapw_band<lapw_impl,lapw_diag>(ks);
        
        memcpy(evalsv_, &ks->evalsv[0], lapw_global.nstsv * sizeof(double));
    }
    
    if (flags & 2)
    {
        mdarray<double,6> densmt(densmt_, lapw_global.nrfmtmax, lapw_global.nrfmtmax, lapw_global.lmmaxvr, lapw_global.atoms.size(), lapw_global.nspinor, lapw_global.nspinor); 
        mdarray<double,3> densir(densir_, lapw_global.ngrtot, lapw_global.nspinor, lapw_global.nspinor);
        ks->occsv.resize(lapw_global.nstsv);
        memcpy(&ks->occsv[0], occsv_, lapw_global.nstsv * sizeof(double));
         
        lapw_density(ks, densmt, densir);
    }
}



