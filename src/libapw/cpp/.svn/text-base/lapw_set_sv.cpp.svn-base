#include "lapw.h"

void b_dot_wf(lapw_eigen_states& eigen_states, mdarray<complex16,3>& hwf)
{
    timer t("b_dot_wf");
    
    unsigned int ngk = eigen_states.kp->ngk;
    
    size_t szmax = 0;
    for (unsigned int is = 0; is < geometry.species.size(); is++)
      szmax = std::max(geometry.species[is]->ci.size(), szmax);
    
    mdarray<complex16,3> zm(NULL, szmax, szmax, p.ndmag);
    zm.allocate();
            
    for (unsigned int ias = 0; ias < p.natmtot; ias++)
    {
        int offset = geometry.atoms[ias].offset_wfmt;
        int sz = geometry.atoms[ias].species->ci.size();
        
        zm.zero();

        for (int j2 = 0; j2 < sz; j2++)
        {
            int lm2 = geometry.atoms[ias].species->ci[j2].lm;
            int idxrf2 = geometry.atoms[ias].species->ci[j2].idxrf;
            
            for (unsigned int i = 0; i < p.ndmag; i++)
            {
                for (int j1 = 0; j1 <= j2; j1++)
                {
                    int lm1 = geometry.atoms[ias].species->ci[j1].lm;
                    int idxrf1 = geometry.atoms[ias].species->ci[j1].idxrf;
                
                    L3_sum_gntyry(lm1, lm2, &p.beffrad(0, idxrf1, idxrf2, ias, i), zm(j1, j2, i));
                }
            }
        }
        // compute hwf = hwf + B_z*|wf_j>
        zhemm<cpu>(0, 0, sz, p.nstfv, zone, &zm(0, 0, 0), zm.size(0), 
            &eigen_states.scalar_wf(offset, 0), eigen_states.scalar_wf.size(0), zone, &hwf(offset, 0, 0), hwf.size(0));
        
        // compute hwf = hwf + (B_x - iB_y)|wf_j>
        if (p.ndmag == 3)
        {
            for (int j2 = 0; j2 < sz; j2++)
            {
                for (int j1 = 0; j1 <= j2; j1++)
                    zm(j1, j2, 0) = zm(j1, j2, 1) - zi * zm(j1, j2, 2);
                
                for (int j1 = j2 + 1; j1 < sz; j1++)
                    zm(j1, j2, 0) = conj(zm(j2, j1, 1)) - zi * conj(zm(j2, j1, 2));
            }
              
            zgemm<cpu>(0, 0, sz, p.nstfv, sz, zone, &zm(0, 0, 0), zm.size(0), 
                &eigen_states.scalar_wf(offset, 0), eigen_states.scalar_wf.size(0), zone, &hwf(offset, 0, 2), hwf.size(0));
        }
    }

    timer *t1 = new timer("b_dot_wf_it");
#pragma omp parallel default(shared)
{        
    std::vector<complex16> wfr(p.ngrtot);
    std::vector<complex16> zfft(p.ngrtot);
#pragma omp for
    for (unsigned int i = 0; i < p.nstfv; i++)
    {
        memset(&wfr[0], 0, p.ngrtot * sizeof(complex16));
        for (unsigned int ig = 0; ig < ngk; ig++) 
            wfr[eigen_states.kp->idxgfft[ig]] = eigen_states.scalar_wf(p.size_wfmt + ig, i);
                                    
        lapw_fft(1, &wfr[0]);
                   
        for (unsigned int ir = 0; ir < p.ngrtot; ir++)
            zfft[ir] = wfr[ir] * p.beffir(ir, 0) * p.cfunir[ir];
                                                           
        lapw_fft(-1, &zfft[0]);
        
        for (unsigned int ig = 0; ig < ngk; ig++) 
            hwf(p.size_wfmt + ig, i, 0) += zfft[eigen_states.kp->idxgfft[ig]];

        if (p.ndmag == 3)
        {
            for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                zfft[ir] = wfr[ir] * (p.beffir(ir, 1) - zi * p.beffir(ir, 2)) * p.cfunir[ir];
                                                               
            lapw_fft(-1, &zfft[0]);
            
            for (unsigned int ig = 0; ig < ngk; ig++) 
                hwf(p.size_wfmt + ig, i, 2) += zfft[eigen_states.kp->idxgfft[ig]];
        }
    }
}
    delete t1;

    // copy -B_z|wf>
    for (unsigned int i = 0; i < p.nstfv; i++)
        for (unsigned int j = 0; j < p.size_wfmt + ngk; j++)
            hwf(j, i, 1) = -hwf(j, i, 0);
}

void lapw_set_sv(lapw_eigen_states& eigen_states)
{
    timer t("lapw_set_sv");
    
    unsigned int wf_size = eigen_states.scalar_wf.size(0);

    int nhwf;
    if (p.ndmag == 0) nhwf = 1; // have only one block, nonmagnetic
    if (p.ndmag == 1) nhwf = 2; // have up-up and dn-dn blocks, collinear
    if (p.ndmag == 3) nhwf = 3; // have up-up, dn-dn and up-dn blocks, general case

    // product of the second-variational hamiltonian and a wave-function
    mdarray<complex16,3> hwf(NULL, wf_size, p.nstfv, nhwf);
    hwf.allocate();
    hwf.zero();
    
    // compute product of magnetic field and wave-function 
    if (p.ndmag > 0)
        b_dot_wf(eigen_states, hwf);
    
    // compute <wf_i | (h * wf_j)> for up-up block
    zgemm<cpu>(2, 0, p.nstfv, p.nstfv, wf_size, zone, &eigen_states.scalar_wf(0, 0), wf_size, 
        &hwf(0, 0, 0), wf_size, zzero, &eigen_states.evecsv(0, 0), p.nstsv);
        
    // compute <wf_i | (h * wf_j)> for dn-dn block
    if (p.ndmag != 0)
        zgemm<cpu>(2, 0, p.nstfv, p.nstfv, wf_size, zone, &eigen_states.scalar_wf(0, 0), wf_size, 
            &hwf(0, 0, 1), wf_size, zzero, &eigen_states.evecsv(p.nstfv, p.nstfv), p.nstsv);

    // compute <wf_i | (h * wf_j)> for up-dn block
    if (p.ndmag == 3)
        zgemm<cpu>(2, 0, p.nstfv, p.nstfv, wf_size, zone, &eigen_states.scalar_wf(0, 0), wf_size,
            &hwf(0, 0, 2), wf_size, zzero, &eigen_states.evecsv(0, p.nstfv), p.nstsv);

    unsigned int nspn = (p.ndmag == 0) ? 1 : 2;
    for (unsigned int ispn = 0, i = 0; ispn < nspn; ispn++)
        for (unsigned int ist = 0; ist < p.nstfv; ist++, i++)
            eigen_states.evecsv(i, i) += eigen_states.evalfv[ist];
}


