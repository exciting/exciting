#include "lapw.h"

/// compute product of magnetic field with first-variational wave-functions
void apply_magnetic_field(bloch_states_k* ks, mdarray<complex16,3>& hwf)
{
    timer t("lapw_set_sv:apply_magnetic_field");
    
    mdarray<complex16,3> zm(NULL, lapw_global.max_mt_index_size, lapw_global.max_mt_index_size, lapw_global.ndmag);
    zm.allocate();
            
    for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
    {
        int offset = lapw_global.atoms[ias]->offset_wfmt;
        int sz = lapw_global.atoms[ias]->species->index.size();
        
        zm.zero();

        for (int j2 = 0; j2 < sz; j2++)
        {
            int lm2 = lapw_global.atoms[ias]->species->index[j2].lm;
            int idxrf2 = lapw_global.atoms[ias]->species->index[j2].idxrf;
            
            for (int i = 0; i < lapw_global.ndmag; i++)
            {
                for (int j1 = 0; j1 <= j2; j1++)
                {
                    int lm1 = lapw_global.atoms[ias]->species->index[j1].lm;
                    int idxrf1 = lapw_global.atoms[ias]->species->index[j1].idxrf;
                
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

    timer *t1 = new timer("lapw_set_sv:apply_magnetic_field:it");
    #pragma omp parallel default(shared)
    {        
        std::vector<complex16> wfr(lapw_global.ngrtot);
        std::vector<complex16> zfft(lapw_global.ngrtot);
        
        #pragma omp for
        for (int i = 0; i < lapw_global.nstfv; i++)
        {
            memset(&wfr[0], 0, lapw_global.ngrtot * sizeof(complex16));
            for (int ig = 0; ig < ks->ngk; ig++) 
                wfr[ks->idxgfft[ig]] = ks->scalar_wave_functions(lapw_global.size_wfmt + ig, i);
                                        
            lapw_fft(1, &wfr[0]);
                       
            for (int ir = 0; ir < lapw_global.ngrtot; ir++)
                zfft[ir] = wfr[ir] * lapw_runtime.beffir(ir, 0) * lapw_global.cfunir[ir];
                                                               
            lapw_fft(-1, &zfft[0]);
            
            for (int ig = 0; ig < ks->ngk; ig++) 
                hwf(lapw_global.size_wfmt + ig, i, 0) += zfft[ks->idxgfft[ig]];
    
            if (lapw_global.ndmag == 3)
            {
                for (int ir = 0; ir < lapw_global.ngrtot; ir++)
                    zfft[ir] = wfr[ir] * (lapw_runtime.beffir(ir, 1) - zi * lapw_runtime.beffir(ir, 2)) * lapw_global.cfunir[ir];
                                                                   
                lapw_fft(-1, &zfft[0]);
                
                for (int ig = 0; ig < ks->ngk; ig++) 
                    hwf(lapw_global.size_wfmt + ig, i, 2) += zfft[ks->idxgfft[ig]];
            }
        }
    }
    delete t1;

    // copy -B_z|wf> TODO: this implementation assumes that hwf was zero on input!!!
    for (int i = 0; i < lapw_global.nstfv; i++)
        for (int j = 0; j < ks->wave_function_size; j++)
            hwf(j, i, 1) = -hwf(j, i, 0);
}

template<spin_block sblock> 
void apply_u_correction(bloch_states_k* ks, mdarray<complex16,3>& hwf)
{
    timer t("lapw_set_sv:apply_u_correction");

    for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
    {
        int offset = lapw_global.atoms[ias]->offset_wfmt;
        int l = lapw_global.atoms[ias]->species->lu;
        Species *sp = lapw_global.atoms[ias]->species;
        
        if (l >= 0)
        {
            int nrf = sp->radial_index.nrf(l);

            for (int ist = 0; ist < lapw_global.nstfv; ist++)
            {
                for (int io2 = 0; io2 < nrf; io2++)
                {
                    for (int lm2 = mt_index::idxlm(l, -l); lm2 <= mt_index::idxlm(l, l); lm2++)
                    {
                        int idx2 = sp->index(lm2, io2);

                        for (int io1 = 0; io1 < nrf; io1++)
                        {
                            for (int lm1 = mt_index::idxlm(l, -l); lm1 <= mt_index::idxlm(l, l); lm1++)
                            {
                                int idx1 = sp->index(lm1, io1);
                                complex16 zt = lapw_runtime.ovlprad(l, io2, io1, ias) * ks->scalar_wave_functions(offset + idx1, ist);

                                if (sblock == uu)
                                    hwf(offset + idx2, ist, 0) += lapw_runtime.vmatu(lm2, lm1, 0, 0, ias) * zt;

                                if (sblock == dd)
                                    hwf(offset + idx2, ist, 1) += lapw_runtime.vmatu(lm2, lm1, 1, 1, ias) * zt; 

                                if (sblock == ud)
                                    hwf(offset + idx2, ist, 2) += lapw_runtime.vmatu(lm2, lm1, 0, 1, ias) * zt; 

                            }
                        }
                    }
                }
            }
        } // l >= 0
    }
}

void apply_so_correction(bloch_states_k* ks, mdarray<complex16,3>& hwf)
{
    timer t("lapw_set_sv:apply_so_correction");
    
    for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
    {
        int offset = lapw_global.atoms[ias]->offset_wfmt;
        Species *sp = lapw_global.atoms[ias]->species;
        
        for (int l = 0; l <= lapw_global.lmaxapw; l++)
        {
            int nrf = sp->radial_index.nrf(l);
            for (int io1 = 0; io1 < nrf; io1++)
            {
                for (int io2 = 0; io2 < nrf; io2++)
                {
                    for (int m = -l; m <= l; m++)
                    {
                        int idx1 = sp->index(l, m, io1);
                        int idx2 = sp->index(l, m, io2);
                        int idx3 = 0;
                        if (m + l) idx3 = sp->index(l, m - 1, io2);
                        
                        for (int ist = 0; ist < lapw_global.nstfv; ist++)
                        {
                            hwf(offset + idx1, ist, 0) += ks->scalar_wave_functions(offset + idx2, ist) * double(m) * lapw_runtime.socrad(l, io1, io2, ias);
                            hwf(offset + idx1, ist, 1) -= ks->scalar_wave_functions(offset + idx2, ist) * double(m) * lapw_runtime.socrad(l, io1, io2, ias);
                            if (m + l) hwf(offset + idx1, ist, 2) += ks->scalar_wave_functions(offset + idx3, ist) * sqrt(double((l + m) * (l - m + 1))) * 
                                lapw_runtime.socrad(l, io1, io2, ias);
                        }
                    }
                }
            }
        }
    }
}


void lapw_set_sv(bloch_states_k *ks)
{
    timer t("lapw_set_sv");
    
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
        apply_magnetic_field(ks, hwf);

    if (lapw_global.ldapu)
    {
        apply_u_correction<uu>(ks, hwf);
        if (lapw_global.ndmag != 0) apply_u_correction<dd>(ks, hwf);
        if (lapw_global.ndmag == 3) apply_u_correction<ud>(ks, hwf);
    }

    if (lapw_global.spinorb)
    {
       apply_so_correction(ks, hwf);
    }
    
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

    for (int ispn = 0, i = 0; ispn < lapw_global.nspinor; ispn++)
        for (int ist = 0; ist < lapw_global.nstfv; ist++, i++)
            ks->evecsv(i, i) += ks->evalfv[ist];
}


