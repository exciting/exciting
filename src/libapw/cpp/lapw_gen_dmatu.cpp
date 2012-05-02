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
    for (int ias = 0; ias < lapw_global.atoms.size(); ias++)
    {
        int l = lapw_global.atoms[ias]->species->lu;

        if (l >= 0)
        {
            mdarray<int,2> *ci_by_lmo = &lapw_global.atoms[ias]->species->ci_by_lmo;
            int ordmax = lapw_global.atoms[ias]->species->rfmt_order[l];

            for (int ispn2 = 0; ispn2 < lapw_global.nspinor; ispn2++)
            {
                for (int ispn1 = 0; ispn1 < lapw_global.nspinor; ispn1++)
                {
                    for (int io2 = 0; io2 < ordmax; io2++)
                    {
                        for (int m2 = -l; m2 <= l; m2++)
                        {
                            int lm2 = idxlm(l, m2);
                            for (int io1 = 0; io1 < ordmax; io1++)
                            {
                                for (int m1 = -l; m1 <= l; m1++)
                                {
                                    int lm1 = idxlm(l, m1);
                                    lapw_runtime.dmatu(lm1, lm2, ispn1, ispn2, ias) += 
                                        zdens((*ci_by_lmo)(lm1, io1), (*ci_by_lmo)(lm2, io2), ispn1, ispn2, ias) * lapw_runtime.ovlprad(l, io1, io2, ias);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
