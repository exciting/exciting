#include "lapw.h"

extern "C" void FORTRAN(lapw_load_species)(int *is_, int *nlorb, int *lorbl, int *apword,
                                           double* rmt, int *nrmt)
{
    int is = *is_ - 1;
    
    lapw_global.species[is]->rmt = *rmt;
    lapw_global.species[is]->nrmt = *nrmt;
    
    for (int i = 0; i < *nlorb; i++) 
        lapw_global.species[is]->lo_descriptors.push_back(radial_l_descriptor(lorbl[i]));
    
    for (unsigned int l = 0; l <= lapw_global.lmaxapw; l++) 
    {
        radial_l_descriptor lch(l);
        for (int io = 0; io < apword[l]; io++)
            lch.radial_solution_descriptors.push_back(radial_solution_descriptor());
        
        lapw_global.species[is]->apw_descriptors.push_back(lch);
    }
}


