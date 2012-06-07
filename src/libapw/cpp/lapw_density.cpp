#include "lapw.h"

void lapw_gen_dmatu(const bloch_states_k* const ks, mdarray<complex16,5>& zdens);

template <int N> void sum_zdens(int ias, int lm3, mdarray<complex16,5>& zdens, mdarray<complex16,2>& gtmp, 
                                mdarray<double,6>& densmt)
{
    int sz = lapw_global.atoms[ias]->species->index.size();
    
    int l3 = lapw_global.l_by_lm[lm3];
    
    for (int j2 = 0; j2 < sz; j2++)
    {
        int l2 = lapw_global.atoms[ias]->species->index[j2].l;
        int lm2 = lapw_global.atoms[ias]->species->index[j2].lm;
        int idxrf2 = lapw_global.atoms[ias]->species->index[j2].idxrf;

        int j1 = 0;
        
        for (int idxrf1 = 0; idxrf1 <=idxrf2; idxrf1++)
        {
            int l1 = lapw_global.atoms[ias]->species->radial_index[idxrf1].l;
            
            if ((l1 + l2 + l3) % 2 == 0)
            {
                for (int lm1 = mt_index::idxlm(l1, -l1); lm1 <= mt_index::idxlm(l1, l1); lm1++, j1++) 
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

void lapw_density(bloch_states_k *ks, mdarray<double,6>& densmt, mdarray<double,3>& densir)
{
    timer t("lapw_density");
    
    std::vector<int> idxocc;
    std::vector<double> woccsv;
    for (int j = 0; j < lapw_global.nstsv; j++)
    {
        double wo = ks->occsv[j] * ks->weight;

        if (wo > 1e-14)
        {
            idxocc.push_back(j);
            woccsv.push_back(wo);
        }
    }
    
    mdarray<complex16,5> zdens(NULL, lapw_global.max_mt_index_size, lapw_global.max_mt_index_size, lapw_global.nspinor, lapw_global.nspinor, lapw_global.atoms.size());
    zdens.allocate();
    zdens.zero();
    mdarray<complex16,3> wf1(NULL, lapw_global.max_mt_index_size, idxocc.size(), lapw_global.nspinor);
    wf1.allocate();
    mdarray<complex16,3> wf2(NULL, lapw_global.max_mt_index_size, idxocc.size(), lapw_global.nspinor);
    wf2.allocate();

    timer *t1 = new timer("lapw_density:zdens");    
    for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
    {
        int offset = lapw_global.atoms[ias]->offset_wfmt;
        int sz = lapw_global.atoms[ias]->species->index.size();

        for (int i = 0; i < (int)idxocc.size(); i++)
            for (int ispn = 0; ispn < lapw_global.nspinor; ispn++)
            {
                memcpy(&wf1(0, i, ispn), &ks->spinor_wave_functions(offset, ispn, idxocc[i]), sz * sizeof(complex16));
                for (int k = 0; k < sz; k++) wf2(k, i, ispn) = wf1(k, i, ispn) * woccsv[i];
            }
        
        for (int ispn2 = 0; ispn2 < lapw_global.nspinor; ispn2++)
            for (int ispn1 = 0; ispn1 < lapw_global.nspinor; ispn1++)
                if (use_spin_block(ispn1, ispn2))
                {
                    zgemm<cpu>(0, 2, sz, sz, idxocc.size(), zone, &wf1(0, 0, ispn1), wf1.size(0), 
                        &wf2(0, 0, ispn2), wf2.size(0), zzero, &zdens(0, 0, ispn1, ispn2, ias), zdens.size(0));
                }
    } 
    delete t1;

    t1 = new timer("lapw_density:densmt");
    #pragma omp parallel default(shared)
    {
        mdarray<complex16,2> gtmp(NULL, lapw_global.lmmaxapw, lapw_global.lmmaxapw);
        gtmp.allocate();        
        #pragma omp for
        for (int lm3 = 0; lm3 < lapw_global.lmmaxvr; lm3++)
        {
            zcopy(gtmp.size(), &lapw_global.gntyry(lm3, 0, 0), lapw_global.lmmaxvr, &gtmp(0, 0), 1);

            for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
            {
                sum_zdens<0>(ias, lm3, zdens, gtmp, densmt);
                
                if (lapw_global.ndmag > 0)
                    sum_zdens<1>(ias, lm3, zdens, gtmp, densmt);
                
                if (lapw_global.ndmag == 3)
                    sum_zdens<2>(ias, lm3, zdens, gtmp, densmt);
            }
        }
    }    
    delete t1;
    
    t1 = new timer("lapw_density:densir"); 
    #pragma omp parallel default(shared)
    {
        mdarray<complex16,2> zfft(NULL, lapw_global.ngrtot, lapw_global.nspinor);
        zfft.allocate();
        mdarray<double,3> densir_tmp(NULL, lapw_global.ngrtot, lapw_global.nspinor, lapw_global.nspinor);
        densir_tmp.allocate();
        densir_tmp.zero();
        #pragma omp for
        for (int i = 0; i < (int)idxocc.size(); i++)
        {
            zfft.zero();
            for (int ispn = 0; ispn < lapw_global.nspinor; ispn++)
            {
                for (int ig = 0; ig < ks->ngk; ig++)
                    zfft(ks->idxgfft[ig], ispn) = ks->spinor_wave_functions(lapw_global.size_wfmt + ig, ispn, idxocc[i]);
                lapw_fft(1, &zfft(0, ispn));
            }
            
            for (int ir = 0; ir < lapw_global.ngrtot; ir++)
                densir_tmp(ir, 0, 0) += real(zfft(ir, 0) * conj(zfft(ir, 0))) * woccsv[i] / lapw_global.omega;
            
            if (lapw_global.ndmag > 0)
                for (int ir = 0; ir < lapw_global.ngrtot; ir++)
                    densir_tmp(ir, 1, 1) += real(zfft(ir, 1) * conj(zfft(ir, 1))) * woccsv[i] / lapw_global.omega;
            
            if (lapw_global.ndmag == 3)
            {
                for (int ir = 0; ir < lapw_global.ngrtot; ir++)
                    densir_tmp(ir, 0, 1) += 2.0 * real(zfft(ir, 0) * conj(zfft(ir, 1))) * woccsv[i] / lapw_global.omega;
                
                for (int ir = 0; ir < lapw_global.ngrtot; ir++)
                    densir_tmp(ir, 1, 0) -= 2.0 * imag(zfft(ir, 0) * conj(zfft(ir, 1))) * woccsv[i] / lapw_global.omega;
            }
        }
        #pragma omp critical
        {
            for (int ispn1 = 0; ispn1 < lapw_global.nspinor; ispn1++)
                for (int ispn2 = 0; ispn2 < lapw_global.nspinor; ispn2++)
                    if (use_spin_block(ispn1, ispn2))
                    {
                        for (int ir = 0; ir < lapw_global.ngrtot; ir++)
                            densir(ir, ispn1, ispn2) += densir_tmp(ir, ispn1, ispn2);
                    }
        }
    }
    delete t1;
    
    if (lapw_global.ldapu)
    {
        t1 = new timer("lapw_density:dmatu");
        
        lapw_gen_dmatu(ks, zdens);
        
        delete t1;
    }
 }
