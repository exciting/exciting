#include "lapw.h"

extern "C" void FORTRAN(lapw_load_species)(int *is_, int *nlorb, int *lorbl, int *apword,
                                           double* rmt, int *nrmt)
{
    int is = *is_ - 1;
    
    geometry.species[is]->rmt = *rmt;
    geometry.species[is]->nrmt = *nrmt;
    
    for (int i = 0; i < *nlorb; i++) 
        geometry.species[is]->lo_descriptors.push_back(radial_l_channel_descriptor(lorbl[i]));
    
    for (unsigned int l = 0; l <= p.lmaxapw; l++) 
    {
        radial_l_channel_descriptor lch(l);
        for (int io = 0; io < apword[l]; io++)
        {
            radial_solution_descriptor rs;
            lch.radial_solution_descriptors.push_back(rs);
        }
        geometry.species[is]->apw_descriptors.push_back(lch);
    }
}


