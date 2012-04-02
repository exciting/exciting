#include "lapw.h"

template <int N> void sum_zdens(int ias, int lm3, mdarray<complex16,5>& zdens, mdarray<complex16,2>& gtmp, 
                                mdarray<double,6>& densmt)
{
    int sz = geometry.atoms[ias].species->ci.size();
    
    int l3 = p.l_by_lm[lm3];
    
    for (int j2 = 0; j2 < sz; j2++)
    {
        int l2 = geometry.atoms[ias].species->ci[j2].l;
        int lm2 = geometry.atoms[ias].species->ci[j2].lm;
        unsigned int idxrf2 = geometry.atoms[ias].species->ci[j2].idxrf;

        int j1 = 0;
        
        for (unsigned int idxrf1 = 0; idxrf1 <=idxrf2; idxrf1++)
        {
            int l1 = geometry.atoms[ias].species->l_by_idxrf[idxrf1];
            
            if ((l1 + l2 + l3) % 2 == 0)
            {
                for (int lm1 = l1 * l1; lm1 < (l1 + 1) * (l1 + 1); lm1++, j1++) 
                {
                
                    if (N == 0) 
                        densmt(idxrf1, idxrf2, lm3, ias, 0, 0) += real(zdens(j1, j2, 0, 0, ias) * gtmp(lm1, lm2));
                
                    if (N == 1) 
                        densmt(idxrf1, idxrf2, lm3, ias, 1, 1) += real(zdens(j1, j2, 1, 1, ias) * gtmp(lm1, lm2));
                
                    if (N == 2) 
                    { 
                        densmt(idxrf1, idxrf2, lm3, ias, 0, 1) += 2.0 * real(zdens(j1, j2, 0, 1, ias) * gtmp(lm1, lm2));
                        densmt(idxrf1, idxrf2, lm3, ias, 1, 0) -= 2.0 * imag(zdens(j1, j2, 0, 1, ias) * gtmp(lm1, lm2));
                    }
                }
            } 
            else
            {
                j1 += (2 * l1 + 1);
            }
        }
    }
}

void lapw_density(lapw_eigen_states& eigen_states, mdarray<double,6>& densmt, mdarray<double,3>& densir)
{
    timer t("lapw_density");
    
    size_t szmax = 0;
    for (unsigned int is = 0; is < geometry.species.size(); is++)
      szmax = std::max(geometry.species[is]->ci.size(), szmax);

    std::vector<int> idxocc;
    std::vector<double> woccsv;
    for (unsigned int j = 0; j < p.nstsv; j++)
    {
        double wo = eigen_states.occsv[j] * eigen_states.kp->weight;

        if (wo > 1e-14)
        {
            idxocc.push_back(j);
            woccsv.push_back(wo);
        }
    }
    
    mdarray<complex16,5> zdens(NULL, szmax, szmax, p.nspinor, p.nspinor, p.natmtot);
    zdens.allocate();
    mdarray<complex16,3> wf1(NULL, szmax, idxocc.size(), p.nspinor);
    wf1.allocate();
    mdarray<complex16,3> wf2(NULL, szmax, idxocc.size(), p.nspinor);
    wf2.allocate();

    timer *t1 = new timer("lapw_density:zdens");    
    for (unsigned int ias = 0; ias < p.natmtot; ias++)
    {
        int offset = geometry.atoms[ias].offset_wfmt;
        int sz = geometry.atoms[ias].species->ci.size();

        for (unsigned i = 0; i < idxocc.size(); i++)
            for (unsigned int ispn = 0; ispn < p.nspinor; ispn++)
            {
                memcpy(&wf1(0, i, ispn), &eigen_states.spinor_wf(offset, ispn, idxocc[i]), sz * sizeof(complex16));
                for (int k = 0; k < sz; k++) wf2(k, i, ispn) = wf1(k, i, ispn) * woccsv[i];
            }
        
        for (unsigned int ispn2 = 0; ispn2 < p.nspinor; ispn2++)
            for (unsigned int ispn1 = 0; ispn1 < p.nspinor; ispn1++)
                if ((p.ndmag == 1 && ispn1 == ispn2) || (p.ndmag != 1))
                    zgemm<cpu>(0, 2, sz, sz, idxocc.size(), zone, &wf1(0, 0, ispn1), wf1.size(0), 
                        &wf2(0, 0, ispn2), wf2.size(0), zzero, &zdens(0, 0, ispn1, ispn2, ias), zdens.size(0));
    }
    delete t1;
    
    t1 = new timer("lapw_density:densmt");
    #pragma omp parallel default(shared)
    {
        mdarray<complex16,2> gtmp(NULL, p.lmmaxapw, p.lmmaxapw);
        gtmp.allocate();        
        #pragma omp for
        for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++)
        {
            zcopy(gtmp.size(), &p.gntyry(lm3, 0, 0), p.lmmaxvr, &gtmp(0, 0), 1);

            for (unsigned int ias = 0; ias < p.natmtot; ias++)
            {
                sum_zdens<0>(ias, lm3, zdens, gtmp, densmt);
                
                if (p.ndmag > 0)
                    sum_zdens<1>(ias, lm3, zdens, gtmp, densmt);
                
                if (p.ndmag == 3)
                    sum_zdens<2>(ias, lm3, zdens, gtmp, densmt);
            }
        }
    }    
    delete t1;
    
    t1 = new timer("lapw_density:densir"); 
    #pragma omp parallel default(shared)
    {
        mdarray<complex16,2> zfft(NULL, p.ngrtot, p.nspinor);
        zfft.allocate();
        mdarray<double,3> densir_tmp(NULL, p.ngrtot, p.nspinor, p.nspinor);
        densir_tmp.allocate();
        densir_tmp.zero();
        #pragma omp for
        for (unsigned int i = 0; i < idxocc.size(); i++)
        {
            zfft.zero();
            for (unsigned int ispn = 0; ispn < p.nspinor; ispn++)
            {
                for (unsigned int ig = 0; ig < eigen_states.kp->ngk; ig++)
                    zfft(eigen_states.kp->idxgfft[ig], ispn) = eigen_states.spinor_wf(p.size_wfmt + ig, ispn, idxocc[i]);
                lapw_fft(1, &zfft(0, ispn));
            }
            
            for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                densir_tmp(ir, 0, 0) += real(zfft(ir, 0) * conj(zfft(ir, 0))) * woccsv[i] / geometry.omega;
            
            if (p.ndmag > 0)
                for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                    densir_tmp(ir, 1, 1) += real(zfft(ir, 1) * conj(zfft(ir, 1))) * woccsv[i] / geometry.omega;
            
            if (p.ndmag == 3)
            {
                for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                    densir_tmp(ir, 0, 1) += 2.0 * real(zfft(ir, 0) * conj(zfft(ir, 1))) * woccsv[i] / geometry.omega;
                
                for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                    densir_tmp(ir, 1, 0) -= 2.0 * imag(zfft(ir, 0) * conj(zfft(ir, 1))) * woccsv[i] / geometry.omega;
            }
        }
        #pragma omp critical
        {
            for (unsigned int ispn1 = 0; ispn1 < p.nspinor; ispn1++)
                for (unsigned int ispn2 = 0; ispn2 < p.nspinor; ispn2++)
                    for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                        densir(ir, ispn1, ispn2) += densir_tmp(ir, ispn1, ispn2);
        }
    }
    delete t1;
}
