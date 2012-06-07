#include "lapw.h"

/*! \file lapw_init.cpp
    \brief initialize internal variables
*/

extern "C" void FORTRAN(lapw_init)()
{
    lapw_global.max_mt_index_size = 0;

    for (int is = 0; is < (int)lapw_global.species.size(); is++)
    {
        for (int l = 0; l <= lapw_global.lmaxapw; l++)
        {
            for (int io = 0; io < (int)lapw_global.species[is]->apw_descriptors[l].radial_solution_descriptors.size(); io++)
                lapw_global.species[is]->radial_index.add(l);
        }
        
        for (int ilo = 0; ilo < (int)lapw_global.species[is]->lo_descriptors.size(); ilo++)
        {
            int l = lapw_global.species[is]->lo_descriptors[ilo].l;
            
            lapw_global.species[is]->radial_index.add(l, ilo);

        }
            
        lapw_global.species[is]->radial_index.init();
        lapw_global.species[is]->index.init(lapw_global.species[is]->radial_index);
        lapw_global.max_mt_index_size = std::max(lapw_global.species[is]->index.size(), lapw_global.max_mt_index_size);
    }

    lapw_global.size_wfmt_apw = 0;
    lapw_global.size_wfmt_lo = 0;
    lapw_global.size_wfmt = 0;
    
    for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
    {
        Atom *atom = lapw_global.atoms[ias];
        Species *species = atom->species;
        
        atom->offset_apw = lapw_global.size_wfmt_apw;
        atom->offset_lo = lapw_global.size_wfmt_lo;
        atom->offset_wfmt = lapw_global.size_wfmt;

        lapw_global.size_wfmt_apw += species->index.apw_size();
        lapw_global.size_wfmt_lo += species->index.lo_size();
        lapw_global.size_wfmt += species->index.size();
    }
    
    /*lapw_global.size_apwalm = 0;
    for (unsigned int ic = 0; ic < lapw_global.natmcls; ic++)
    {
        int ias = lapw_global.ic2ias[ic];
        Species *species = lapw_global.atoms[ias].species;
        lapw_global.size_apwalm += species->size_ci_apw;
    }*/
    
    lapw_global.l_by_lm.resize(lapw_global.lmmaxapw);
    for (int l = 0; l <= lapw_global.lmaxapw; l++)
        for (int m = -l; m <= l; m++)
            lapw_global.l_by_lm[mt_index::idxlm(l, m)] = l;

    assert(lapw_global.size_wfmt == (lapw_global.size_wfmt_apw + lapw_global.size_wfmt_lo));
    
/*    for (unsigned int is = 0; is < lapw_global.species.size(); is++)
    {
        std::cout << "species : " << is << std::endl 
                  << "  size_ci : " << lapw_global.species[is].ci.size() << std::endl
                  << "  size_ci_apw : " << lapw_global.species[is].size_ci_apw << std::endl 
                  << "  size_ci_lo : " << lapw_global.species[is].size_ci_lo << std::endl;
    }
    
    for (unsigned ias = 0; ias < lapw_global.atoms.size(); ias++)
    {
        Atom *atom = &lapw_global.atoms[ias];
        std::cout << "atom : " << ias << std::endl
                  << "  offset_wfmt : " << atom->offset_wfmt << std::endl
                  << "  offset_apw : " << atom->offset_apw << std::endl
                  << "  offset_lo : " << atom->offset_lo << std::endl;
    }
*/    
}

