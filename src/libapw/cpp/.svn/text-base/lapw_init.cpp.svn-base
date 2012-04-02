#include "lapw.h"

extern "C" void FORTRAN(lapw_init)()
{
    for (unsigned int is = 0; is < geometry.species.size(); is++)
    {
        geometry.species[is]->rfmt_order.resize(p.lmaxapw + 1, 0);
        
        for (unsigned int l = 0; l <= p.lmaxapw; l++)
        {
            for (unsigned int io = 0; io < geometry.species[is]->apw_descriptors[l].radial_solution_descriptors.size(); io++)
            {
                geometry.species[is]->ci_by_idxrf.push_back(geometry.species[is]->ci.size());
                geometry.species[is]->l_by_idxrf.push_back(l);
                
                for (int m = -l; m <= (int)l; m++)
                    geometry.species[is]->ci.push_back(mtci(l, m, geometry.species[is]->rfmt_order[l], geometry.species[is]->nrfmt));
                
                geometry.species[is]->rfmt_order[l]++;
                geometry.species[is]->nrfmt++;
            }
        }
        geometry.species[is]->size_ci_apw = geometry.species[is]->ci.size();
        
        for (unsigned int ilo = 0; ilo < geometry.species[is]->lo_descriptors.size(); ilo++)
        {
            int l = geometry.species[is]->lo_descriptors[ilo].l;
            
            geometry.species[is]->ci_by_idxrf.push_back(geometry.species[is]->ci.size());
            geometry.species[is]->l_by_idxrf.push_back(l);

            for (int m = -l; m <= l; m++)
                geometry.species[is]->ci.push_back(mtci(l, m, geometry.species[is]->rfmt_order[l], geometry.species[is]->nrfmt, ilo));
            
            geometry.species[is]->rfmt_order[l]++;
            geometry.species[is]->nrfmt++;
        }
        geometry.species[is]->size_ci_lo = geometry.species[is]->ci.size() - geometry.species[is]->size_ci_apw;
        geometry.species[is]->ci_lo = &geometry.species[is]->ci[geometry.species[is]->size_ci_apw];

        int maxorder = 0;
        for (unsigned int l = 0; l <= p.lmaxapw; l++)
            maxorder = std::max(maxorder, geometry.species[is]->rfmt_order[l]);
        
        geometry.species[is]->ci_by_lmo.set_dimensions(p.lmmaxapw, maxorder);
        geometry.species[is]->ci_by_lmo.allocate();

        for (unsigned int i = 0; i < geometry.species[is]->ci.size(); i++)
            geometry.species[is]->ci_by_lmo(geometry.species[is]->ci[i].lm, geometry.species[is]->ci[i].order) = i;
    }

    p.size_wfmt_apw = 0;
    p.size_wfmt_lo = 0;
    p.size_wfmt = 0;
    
    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        Species *species = atom->species;
        
        atom->offset_apw = p.size_wfmt_apw;
        atom->offset_lo = p.size_wfmt_lo;
        atom->offset_wfmt = p.size_wfmt;

        p.size_wfmt_apw += species->size_ci_apw;
        p.size_wfmt_lo += species->size_ci_lo;
        p.size_wfmt += species->ci.size();
    }
    
    /*p.size_apwalm = 0;
    for (unsigned int ic = 0; ic < p.natmcls; ic++)
    {
        int ias = p.ic2ias[ic];
        Species *species = geometry.atoms[ias].species;
        p.size_apwalm += species->size_ci_apw;
    }*/
    
    p.l_by_lm.resize(p.lmmaxapw);
    for (unsigned int l = 0; l <= p.lmaxapw; l++)
        for (int m = -l; m <= (int)l; m++)
            p.l_by_lm[idxlm(l, m)] = l;

    assert(p.size_wfmt == (p.size_wfmt_apw + p.size_wfmt_lo));
    
/*    for (unsigned int is = 0; is < geometry.species.size(); is++)
    {
        std::cout << "species : " << is << std::endl 
                  << "  size_ci : " << geometry.species[is].ci.size() << std::endl
                  << "  size_ci_apw : " << geometry.species[is].size_ci_apw << std::endl 
                  << "  size_ci_lo : " << geometry.species[is].size_ci_lo << std::endl;
    }
    
    for (unsigned ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        std::cout << "atom : " << ias << std::endl
                  << "  offset_wfmt : " << atom->offset_wfmt << std::endl
                  << "  offset_apw : " << atom->offset_apw << std::endl
                  << "  offset_lo : " << atom->offset_lo << std::endl;
    }
*/    
}

