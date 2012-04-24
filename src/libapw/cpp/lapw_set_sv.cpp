#include "lapw.h"

/// compute product of magnetic field with first-variational wave-functions
void b_dot_wf(bloch_states_k *ks, mdarray<complex16,3>& hwf)
{
    timer t("b_dot_wf");
    
    size_t szmax = 0;
    for (unsigned int is = 0; is < lapw_global.species.size(); is++)
      szmax = std::max(lapw_global.species[is]->ci.size(), szmax);
    
    mdarray<complex16,3> zm(NULL, szmax, szmax, lapw_global.ndmag);
    zm.allocate();
            
    for (unsigned int ias = 0; ias < lapw_global.natmtot; ias++)
    {
        int offset = lapw_global.atoms[ias]->offset_wfmt;
        int sz = lapw_global.atoms[ias]->species->ci.size();
        
        zm.zero();

        for (int j2 = 0; j2 < sz; j2++)
        {
            int lm2 = lapw_global.atoms[ias]->species->ci[j2].lm;
            int idxrf2 = lapw_global.atoms[ias]->species->ci[j2].idxrf;
            
            for (unsigned int i = 0; i < lapw_global.ndmag; i++)
            {
                for (int j1 = 0; j1 <= j2; j1++)
                {
                    int lm1 = lapw_global.atoms[ias]->species->ci[j1].lm;
                    int idxrf1 = lapw_global.atoms[ias]->species->ci[j1].idxrf;
                
                    L3_sum_gntyry(lm1, lm2, &lapw_runtime.beffrad(0, idxrf1, idxrf2, ias, i), zm(j1, j2, i));
                }
            }
        }
        // compute hwf = hwf + B_z*|wf_j>
        zhemm<cpu>(0, 0, sz, lapw_global.nstfv, zone, &zm(0, 0, 0), zm.size(0), 
            &ks->scalar_wave_functions(offset, 0), ks->scalar_wave_functions.size(0), zone, &hwf(offset, 0, 0), hwf.size(0));
        
        // compute hwf = hwf + (B_x - iB_y)|wf_j>
        if (lapw_global.ndmag == 3)
        {
            for (int j2 = 0; j2 < sz; j2++)
            {
                for (int j1 = 0; j1 <= j2; j1++)
                    zm(j1, j2, 0) = zm(j1, j2, 1) - zi * zm(j1, j2, 2);
                
                for (int j1 = j2 + 1; j1 < sz; j1++)
                    zm(j1, j2, 0) = conj(zm(j2, j1, 1)) - zi * conj(zm(j2, j1, 2));
            }
              
            zgemm<cpu>(0, 0, sz, lapw_global.nstfv, sz, zone, &zm(0, 0, 0), zm.size(0), 
                &ks->scalar_wave_functions(offset, 0), ks->scalar_wave_functions.size(0), zone, &hwf(offset, 0, 2), hwf.size(0));
        }
    }

    timer *t1 = new timer("b_dot_wf_it");
#pragma omp parallel default(shared)
{        
    std::vector<complex16> wfr(lapw_global.ngrtot);
    std::vector<complex16> zfft(lapw_global.ngrtot);
#pragma omp for
    for (unsigned int i = 0; i < lapw_global.nstfv; i++)
    {
        memset(&wfr[0], 0, lapw_global.ngrtot * sizeof(complex16));
        for (unsigned int ig = 0; ig < ks->ngk; ig++) 
            wfr[ks->idxgfft[ig]] = ks->scalar_wave_functions(lapw_global.size_wfmt + ig, i);
                                    
        lapw_fft(1, &wfr[0]);
                   
        for (unsigned int ir = 0; ir < lapw_global.ngrtot; ir++)
            zfft[ir] = wfr[ir] * lapw_runtime.beffir(ir, 0) * lapw_global.cfunir[ir];
                                                           
        lapw_fft(-1, &zfft[0]);
        
        for (unsigned int ig = 0; ig < ks->ngk; ig++) 
            hwf(lapw_global.size_wfmt + ig, i, 0) += zfft[ks->idxgfft[ig]];

        if (lapw_global.ndmag == 3)
        {
            for (unsigned int ir = 0; ir < lapw_global.ngrtot; ir++)
                zfft[ir] = wfr[ir] * (lapw_runtime.beffir(ir, 1) - zi * lapw_runtime.beffir(ir, 2)) * lapw_global.cfunir[ir];
                                                               
            lapw_fft(-1, &zfft[0]);
            
            for (unsigned int ig = 0; ig < ks->ngk; ig++) 
                hwf(lapw_global.size_wfmt + ig, i, 2) += zfft[ks->idxgfft[ig]];
        }
    }
}
    delete t1;

    // copy -B_z|wf>
    for (unsigned int i = 0; i < lapw_global.nstfv; i++)
        for (unsigned int j = 0; j < lapw_global.size_wfmt + ks->ngk; j++)
            hwf(j, i, 1) = -hwf(j, i, 0);
}

void lapw_set_sv(bloch_states_k *ks)
{
    timer t("lapw_set_sv");
    
    //unsigned int wf_size = eigen_states.scalar_wf.size(0);

    int nhwf;
    if (lapw_global.ndmag == 0) nhwf = 1; // have only one block, nonmagnetic
    if (lapw_global.ndmag == 1) nhwf = 2; // have up-up and dn-dn blocks, collinear
    if (lapw_global.ndmag == 3) nhwf = 3; // have up-up, dn-dn and up-dn blocks, general case

    // product of the second-variational hamiltonian and a wave-function
    mdarray<complex16,3> hwf(NULL, ks->wave_function_size, lapw_global.nstfv, nhwf);
    hwf.allocate();
    hwf.zero();
    
    // compute product of magnetic field and wave-function 
    if (lapw_global.ndmag > 0)
        b_dot_wf(ks, hwf);
    
    // compute <wf_i | (h * wf_j)> for up-up block
    zgemm<cpu>(2, 0, lapw_global.nstfv, lapw_global.nstfv, ks->wave_function_size, zone, &ks->scalar_wave_functions(0, 0), 
        ks->scalar_wave_functions.size(0), &hwf(0, 0, 0), hwf.size(0), zzero, &ks->evecsv(0, 0), ks->evecsv.size(0));
        
    // compute <wf_i | (h * wf_j)> for dn-dn block
    if (lapw_global.ndmag != 0)
        zgemm<cpu>(2, 0, lapw_global.nstfv, lapw_global.nstfv, ks->wave_function_size, zone, &ks->scalar_wave_functions(0, 0), 
        ks->scalar_wave_functions.size(0), &hwf(0, 0, 1), hwf.size(0), zzero, &ks->evecsv(lapw_global.nstfv, lapw_global.nstfv), 
        ks->evecsv.size(0));

    // compute <wf_i | (h * wf_j)> for up-dn block
    if (lapw_global.ndmag == 3)
        zgemm<cpu>(2, 0, lapw_global.nstfv, lapw_global.nstfv, ks->wave_function_size, zone, &ks->scalar_wave_functions(0, 0), 
        ks->scalar_wave_functions.size(0), &hwf(0, 0, 2), hwf.size(0), zzero, &ks->evecsv(0, lapw_global.nstfv), ks->evecsv.size(0));

    //unsigned int nspn = (lapw_global.ndmag == 0) ? 1 : 2;
    for (unsigned int ispn = 0, i = 0; ispn < lapw_global.nspinor; ispn++)
        for (unsigned int ist = 0; ist < lapw_global.nstfv; ist++, i++)
            ks->evecsv(i, i) += ks->evalfv[ist];
}


