#include "lapw.h"
#include "lapw_band.h"

extern "C" void FORTRAN(lapw_execute)(int32_t *ikloc_, complex16 *apwalm_, double *evalsv_, 
                                      double *occsv_, double *densmt_, double *densir_, int32_t *flags_)
{
    timer t("lapw_execute");
    
    int flags = *flags_;
    unsigned int ikloc = *ikloc_ - 1;
    
    lapw_eigen_states eigen_states(p.kpoints[ikloc]);
    eigen_states.pack_apwalm(apwalm_);
    
    if (flags & 1) 
    {
        lapw_band<lapw_impl,lapw_diag>(eigen_states);
        
        memcpy(evalsv_, &eigen_states.evalsv[0], p.nstsv * sizeof(double));
    }
    
    if (flags & 2)
    {
        mdarray<double,6> densmt(densmt_, p.nrfmtmax, p.nrfmtmax, p.lmmaxvr, p.natmtot, p.nspinor, p.nspinor); 
        mdarray<double,3> densir(densir_, p.ngrtot, p.nspinor, p.nspinor);
        eigen_states.occsv.resize(p.nstsv);
        memcpy(&eigen_states.occsv[0], occsv_, p.nstsv * sizeof(double));
         
        lapw_density(eigen_states, densmt, densir);
    }
}



