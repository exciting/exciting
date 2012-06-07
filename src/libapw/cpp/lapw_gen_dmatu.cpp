#include "lapw.h"

/*! \file lapw_gen_dmatu.cpp 
    \brief add k-contribution to the density matrix of the LDA+U method
*/

/*! \ingroup functions
    \brief adds k-contribution to the density (otherwise called orbital occuapncy) matrix used in the LDA+U method

     Density matrix has the following expression:
     \f[
        n_{\ell,mm'}^{\sigma \sigma'} = \sum_{i {\bf k}}^{occ} \int_{0}^{R_{MT}} r^2 dr 
          \Psi_{\ell m}^{i{\bf k}\sigma *}({\bf r}) \Psi_{\ell m'}^{i{\bf k}\sigma'}({\bf r})
     \f]
*/
void lapw_gen_dmatu(const bloch_states_k* const ks, mdarray<complex16,5>& zdens)
{
    for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
    {
        Species *sp = lapw_global.atoms[ias]->species;
        int l = sp->lu;

        if (l >= 0)
        {
            int nrf = sp->radial_index.nrf(l);

            for (int ispn2 = 0; ispn2 < lapw_global.nspinor; ispn2++)
            {
                for (int ispn1 = 0; ispn1 < lapw_global.nspinor; ispn1++)
                {
                    if (use_spin_block(ispn1, ispn2))
                    {
                        for (int io2 = 0; io2 < nrf; io2++)
                        {
                            for (int lm2 = mt_index::idxlm(l, -l); lm2 <= mt_index::idxlm(l, l); lm2++)
                            {
                                for (int io1 = 0; io1 < nrf; io1++)
                                {
                                    for (int lm1 = mt_index::idxlm(l, -l); lm1 <= mt_index::idxlm(l, l); lm1++)
                                    {
                                        lapw_runtime.dmatu(lm1, lm2, ispn1, ispn2, ias) += 
                                            zdens(sp->index(lm1, io1), sp->index(lm2, io2), ispn1, ispn2, ias) * lapw_runtime.ovlprad(l, io1, io2, ias);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
