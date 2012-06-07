#include "lapw.h"

extern "C" void FORTRAN(lapw_get_evec)(int32_t *ikloc_, int32_t *nmatmax_, complex16 *evecfv_, complex16 *evecsv_)
{
    int ikloc = *ikloc_ - 1;
    int nmatmax = *nmatmax_;
    
    mdarray<complex16,2> evecfv(evecfv_, nmatmax, lapw_global.nstfv); 
    mdarray<complex16,2> evecsv(evecsv_, lapw_global.nstsv, lapw_global.nstsv);
    
    for (int i = 0; i < lapw_global.nstfv; i++)
    {
        memcpy(&evecfv(0, i), &lapw_runtime.bloch_states[ikloc]->evecfv(0, i), lapw_runtime.bloch_states[ikloc]->lapw_basis_size * sizeof(complex16));
    }

    for (int i = 0; i < lapw_global.nstsv; i++)
    {
        memcpy(&evecsv(0, i), &lapw_runtime.bloch_states[ikloc]->evecsv(0, i), lapw_global.nstsv * sizeof(complex16));
    }
}
