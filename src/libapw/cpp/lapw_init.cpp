#include "lapw.h"

/*! \file lapw_init.cpp
    \brief initialize internal variables
*/

extern "C" void FORTRAN(lapw_init)()
{
    for (unsigned int is = 0; is < lapw_global.species.size(); is++)
    {
        lapw_global.species[is]->rfmt_order.resize(lapw_global.lmaxapw + 1, 0);
        
        for (unsigned int l = 0; l <= lapw_global.lmaxapw; l++)
        {
            for (unsigned int io = 0; io < lapw_global.species[is]->apw_descriptors[l].radial_solution_descriptors.size(); io++)
            {
                lapw_global.species[is]->ci_by_idxrf.push_back(lapw_global.species[is]->ci.size());
                lapw_global.species[is]->l_by_idxrf.push_back(l);
                
                for (int m = -l; m <= (int)l; m++)
                    lapw_global.species[is]->ci.push_back(mtci(l, m, lapw_global.species[is]->rfmt_order[l], lapw_global.species[is]->nrfmt));
                
                lapw_global.species[is]->rfmt_order[l]++;
                lapw_global.species[is]->nrfmt++;
            }
        }
        lapw_global.species[is]->size_ci_apw = lapw_global.species[is]->ci.size();
        
        for (unsigned int ilo = 0; ilo < lapw_global.species[is]->lo_descriptors.size(); ilo++)
        {
            int l = lapw_global.species[is]->lo_descriptors[ilo].l;
            
            lapw_global.species[is]->ci_by_idxrf.push_back(lapw_global.species[is]->ci.size());
            lapw_global.species[is]->l_by_idxrf.push_back(l);

            for (int m = -l; m <= l; m++)
                lapw_global.species[is]->ci.push_back(mtci(l, m, lapw_global.species[is]->rfmt_order[l], lapw_global.species[is]->nrfmt, ilo));
            
            lapw_global.species[is]->rfmt_order[l]++;
            lapw_global.species[is]->nrfmt++;
        }
        lapw_global.species[is]->size_ci_lo = lapw_global.species[is]->ci.size() - lapw_global.species[is]->size_ci_apw;
        lapw_global.species[is]->ci_lo = &lapw_global.species[is]->ci[lapw_global.species[is]->size_ci_apw];

        int maxorder = 0;
        for (unsigned int l = 0; l <= lapw_global.lmaxapw; l++)
            maxorder = std::max(maxorder, lapw_global.species[is]->rfmt_order[l]);
        
        lapw_global.species[is]->ci_by_lmo.set_dimensions(lapw_global.lmmaxapw, maxorder);
        lapw_global.species[is]->ci_by_lmo.allocate();

        for (unsigned int i = 0; i < lapw_global.species[is]->ci.size(); i++)
            lapw_global.species[is]->ci_by_lmo(lapw_global.species[is]->ci[i].lm, lapw_global.species[is]->ci[i].order) = i;
    }

    lapw_global.size_wfmt_apw = 0;
    lapw_global.size_wfmt_lo = 0;
    lapw_global.size_wfmt = 0;
    
    for (unsigned int ias = 0; ias < lapw_global.atoms.size(); ias++)
    {
        Atom *atom = lapw_global.atoms[ias];
        Species *species = atom->species;
        
        atom->offset_apw = lapw_global.size_wfmt_apw;
        atom->offset_lo = lapw_global.size_wfmt_lo;
        atom->offset_wfmt = lapw_global.size_wfmt;

        lapw_global.size_wfmt_apw += species->size_ci_apw;
        lapw_global.size_wfmt_lo += species->size_ci_lo;
        lapw_global.size_wfmt += species->ci.size();
    }
    
    /*lapw_global.size_apwalm = 0;
    for (unsigned int ic = 0; ic < lapw_global.natmcls; ic++)
    {
        int ias = lapw_global.ic2ias[ic];
        Species *species = lapw_global.atoms[ias].species;
        lapw_global.size_apwalm += species->size_ci_apw;
    }*/
    
    lapw_global.l_by_lm.resize(lapw_global.lmmaxapw);
    for (unsigned int l = 0; l <= lapw_global.lmaxapw; l++)
        for (int m = -l; m <= (int)l; m++)
            lapw_global.l_by_lm[idxlm(l, m)] = l;

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

