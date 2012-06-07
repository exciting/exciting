#ifndef __LAPW_BAND_H__
#define __LAPW_BAND_H__

#include "lapw.h"

void lapw_set_sv(bloch_states_k *ks);

/*! \ingroup functions
    \brief applies the muffin-tin part of the first-variational Hamiltonian to the
           APW basis function
    
    The following vector is computed:
    \f[
      b_{L_2 \nu_2}^{\alpha}({\bf G'}) = \sum_{L_1 \nu_1} \sum_{L_3} 
        a_{L_1\nu_1}^{\alpha*}({\bf G'}) 
        \langle u_{\ell_1\nu_1}^{\alpha} | h_{L3}^{\alpha} |  u_{\ell_2\nu_2}^{\alpha}  
        \rangle  \langle Y_{L_1} | R_{L_3} | Y_{L_2} \rangle +  
        \frac{1}{2} \sum_{\nu_1} a_{L_2\nu_1}^{\alpha *}({\bf G'})
        u_{\ell_2\nu_1}^{\alpha}(R_{\alpha})
        u_{\ell_2\nu_2}^{'\alpha}(R_{\alpha})R_{\alpha}^{2}
    \f] 
*/
template <spin_block sblock> 
void apply_hfvmt_to_apw(bloch_states_k* const ks, mdarray<complex16,2>& hapw)
{
    timer t("apply_hfvmt_to_apw");
   
    #pragma omp parallel default(shared)
    {
        std::vector<complex16> zv(ks->ngk);
        std::vector<double> v1(lapw_global.lmmaxvr);
        std::vector<complex16> v2(lapw_global.lmmaxvr);
        #pragma omp for
        for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
        {
            Atom *atom = lapw_global.atoms[ias];
            Species *species = atom->species;
        
            // precompute apw block
            for (int j2 = 0; j2 < (int)species->index.apw_size(); j2++)
            {
                memset(&zv[0], 0, ks->ngk * sizeof(complex16));
                
                int lm2 = species->index[j2].lm;
                int idxrf2 = species->index[j2].idxrf;
                
                for (int j1 = 0; j1 < (int)species->index.apw_size(); j1++)
                {
                    int lm1 = species->index[j1].lm;
                    int idxrf1 = species->index[j1].idxrf;
                    
                    complex16 zsum(0, 0);
                    
                    if (sblock == nm)
                    {
                        L3_sum_gntyry(lm1, lm2, &lapw_runtime.hmltrad(0, idxrf1, idxrf2, ias), zsum);
                    }

                    if (sblock == uu)
                    {
                        for (int lm3 = 0; lm3 < lapw_global.lmmaxvr; lm3++) 
                            v1[lm3] = lapw_runtime.hmltrad(lm3, idxrf1, idxrf2, ias) + lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 0);
                        L3_sum_gntyry(lm1, lm2, &v1[0], zsum);
                    }
                    
                    if (sblock == dd)
                    {
                        for (int lm3 = 0; lm3 < lapw_global.lmmaxvr; lm3++) 
                            v1[lm3] = lapw_runtime.hmltrad(lm3, idxrf1, idxrf2, ias) - lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 0);
                        L3_sum_gntyry(lm1, lm2, &v1[0], zsum);
                    }
                    
                    if (sblock == ud)
                    {
                        for (int lm3 = 0; lm3 < lapw_global.lmmaxvr; lm3++) 
                            v2[lm3] = complex16(lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 1), -lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 2));
                        L3_sum_gntyry(lm1, lm2, &v2[0], zsum);
                    }
        
                    if (abs(zsum) > 1e-14) 
                        for (int ig = 0; ig < ks->ngk; ig++) 
                            zv[ig] += zsum * ks->apwalm(ig, atom->offset_apw + j1); 
                }
                
                // surface term
                if (sblock != ud)
                {
                    int l2 = species->index[j2].l;
                    int io2 = species->index[j2].order;
                    
                    for (int io1 = 0; io1 < (int)species->apw_descriptors[l2].radial_solution_descriptors.size(); io1++)
                    {
                        double t1 = 0.5 * pow(species->rmt, 2) * lapw_runtime.apwfr(species->nrmt - 1, 0, io1, l2, ias) * lapw_runtime.apwdfr(io2, l2, ias); 
                        for (int ig = 0; ig < ks->ngk; ig++) 
                            zv[ig] += t1 * ks->apwalm(ig, atom->offset_apw + species->index(lm2, io1));
                    }
                }

                memcpy(&hapw(0, atom->offset_apw + j2), &zv[0], ks->ngk * sizeof(complex16));
            }
        } 
    }
}

/*! \ingroup functions
    \brief sets up a Hamiltonian matrix
    
    Hamiltonian matrix has the following structure:
    \f[
       H_{\mu' \mu}=\langle \varphi_{\mu' } | \hat H | \varphi_{\mu } \rangle  = 
       \left( \begin{array}{cc} 
         H_{\bf G'G} & H_{{\bf G'}j} \\
         H_{j'{\bf G}} & H_{j'j}
       \end{array} \right)
    \f]
    with APW-APW, APW-lo, lo-APW and lo-lo blocks. Two Hamiltonian 
    representations are used to evaluate muffin-tin and interstitial contribution 
    to the matrix elements: inside muffin-tin spheres it's a spherical expansion
    \f[
      H({\bf r}) = \sum_{L_3} h_{L_3}^{\alpha}(r)  R_{L_3}({\hat {\bf r}}) 
    \f]
    where 
    \f[
      h_{00}^{\alpha}(r) = \frac{1}{R_{00}}\Big(-\frac{1}{2} \nabla^{2} + v_{00}^{\alpha}(r) \Big)
    \f]
    and in the interstitial it's a plane-wave expansion
    \f[
      H({\bf r})=-\frac{1}{2}\nabla^2 + \sum_{{\bf G}} e^{i{\bf Gr}}V({\bf G})
    \f]

      
*/
template <implementation impl, spin_block sblock> 
void lapw_set_h(bloch_states_k* const ks, mdarray<complex16,2>& h)
{
    timer t("lapw_set_h");
   
    mdarray<complex16,2> zm(NULL, ks->ngk, lapw_global.size_wfmt_apw);
    zm.allocate();
    
    apply_hfvmt_to_apw<sblock>(ks, zm);

    if (impl == cpu)
    {
        zgemm<cpu>(0, 2, ks->ngk, ks->ngk, lapw_global.size_wfmt_apw, zone, &zm(0, 0), zm.size(0), 
                   &ks->apwalm(0, 0), ks->apwalm.size(0), zzero, &h(0, 0), h.size(0));
    }
    if (impl == gpu)
    {
        ks->apwalm.allocate_on_device();
        ks->apwalm.copy_to_device();
        zm.allocate_on_device();
        zm.copy_to_device();
        h.allocate_on_device();
        h.zero_on_device();
        zgemm<gpu>(0, 2, ks->ngk, ks->ngk, lapw_global.size_wfmt_apw, zone, zm.get_ptr_device(), zm.size(0),
                   ks->apwalm.get_ptr_device(), ks->apwalm.size(0), zzero, h.get_ptr_device(), h.size(0));
        h.copy_to_host();
        h.deallocate_on_device();
        zm.deallocate_on_device();
    }
 
    #pragma omp parallel default(shared)
    {
        std::vector<double> v1(lapw_global.lmmaxvr);
        std::vector<complex16> v2(lapw_global.lmmaxvr);
        #pragma omp for
        for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
        {
            Atom *atom = lapw_global.atoms[ias];
            Species *species = atom->species;
    
            int lo_index_offset = species->index.apw_size();
            
            for (int j2 = 0; j2 < species->index.lo_size(); j2++) // loop over columns (local-orbital block) 
            {
                int lm2 = species->index[lo_index_offset + j2].lm;
                int idxrf2 = species->index[lo_index_offset + j2].idxrf;
                
                // apw-lo block
                for (int j1 = 0; j1 < species->index.apw_size(); j1++) // loop over rows
                {
                    int lm1 = species->index[j1].lm;
                    int idxrf1 = species->index[j1].idxrf;
                    
                    complex16 zsum(0, 0);
                    
                    if (sblock == nm)
                    {
                        L3_sum_gntyry(lm1, lm2, &lapw_runtime.hmltrad(0, idxrf2, idxrf1, ias), zsum);
                    }

                    if (sblock == uu)
                    {
                        for (int lm3 = 0; lm3 < lapw_global.lmmaxvr; lm3++) 
                            v1[lm3] = lapw_runtime.hmltrad(lm3, idxrf2, idxrf1, ias) + lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 0);
                        L3_sum_gntyry(lm1, lm2, &v1[0], zsum);
                    }
                    
                    if (sblock == dd)
                    {
                        for (int lm3 = 0; lm3 < lapw_global.lmmaxvr; lm3++) 
                            v1[lm3] = lapw_runtime.hmltrad(lm3, idxrf2, idxrf1, ias) - lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 0);
                        L3_sum_gntyry(lm1, lm2, &v1[0], zsum);
                    }
                    
                    if (sblock == ud)
                    {
                        for (int lm3 = 0; lm3 < lapw_global.lmmaxvr; lm3++) 
                            v2[lm3] = complex16(lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 1), -lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 2));
                        L3_sum_gntyry(lm1, lm2, &v2[0], zsum);
                    }
        
                    if (abs(zsum) > 1e-14)
                        for (int ig = 0; ig < ks->ngk; ig++)
                            h(ig, ks->ngk + atom->offset_lo + j2) += zsum * ks->apwalm(ig, atom->offset_apw + j1);
                }
                
                int j1_last = j2;
                if (sblock == ud) j1_last = species->index.lo_size() - 1;
                
                // lo-lo block 
                for (int j1 = 0; j1 <= j1_last; j1++)
                {
                    int lm1 = species->index[lo_index_offset + j1].lm;
                    int idxrf1 = species->index[lo_index_offset + j1].idxrf;
                    
                    complex16 zsum(0, 0);

                    if (sblock == nm)
                    {
                        L3_sum_gntyry(lm1, lm2, &lapw_runtime.hmltrad(0, idxrf1, idxrf2, ias), zsum);
                    }

                    if (sblock == uu)
                    {
                        for (int lm3 = 0; lm3 < lapw_global.lmmaxvr; lm3++) 
                            v1[lm3] = lapw_runtime.hmltrad(lm3, idxrf1, idxrf2, ias) + lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 0);
                        L3_sum_gntyry(lm1, lm2, &v1[0], zsum);
                    }
                    
                    if (sblock == dd)
                    {
                        for (int lm3 = 0; lm3 < lapw_global.lmmaxvr; lm3++) 
                            v1[lm3] = lapw_runtime.hmltrad(lm3, idxrf1, idxrf2, ias) - lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 0);
                        L3_sum_gntyry(lm1, lm2, &v1[0], zsum);
                    }
                    
                    if (sblock == ud)
                    {
                        for (int lm3 = 0; lm3 < lapw_global.lmmaxvr; lm3++) 
                            v2[lm3] = complex16(lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 1), -lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 2));
                        L3_sum_gntyry(lm1, lm2, &v2[0], zsum);
                    }

                    h(ks->ngk + atom->offset_lo + j1, ks->ngk + atom->offset_lo + j2) += zsum;
                }
            }

            if (sblock == ud)
            {
                for (int j2 = 0; j2 < species->index.apw_size(); j2++)
                {
                    int lm2 = species->index[j2].lm;
                    int idxrf2 = species->index[j2].idxrf;
                    
                    for (int j1 = 0; j1 < species->index.lo_size(); j1++)
                    {
                        int lm1 = species->index[lo_index_offset + j1].lm;
                        int idxrf1 = species->index[lo_index_offset + j1].idxrf;
                        
                        complex16 zsum(0, 0);
                        
                        for (int lm3 = 0; lm3 < lapw_global.lmmaxvr; lm3++) 
                            v2[lm3] = complex16(lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 1), -lapw_runtime.beffrad(lm3, idxrf1, idxrf2, ias, 2));
                        L3_sum_gntyry(lm1, lm2, &v2[0], zsum);
                        
                        if (abs(zsum) > 1e-14)
                            for (int ig = 0; ig < ks->ngk; ig++)
                                h(ks->ngk + atom->offset_lo + j1, ig) += zsum * conj(ks->apwalm(ig, atom->offset_apw + j2));
                    }
                }
            }
        } 
    }
    
    for (int j2 = 0; j2 < ks->ngk; j2++) // loop over columns
    {
        int j1_last = j2;
        if (sblock == ud) j1_last = ks->ngk - 1;
   
        double v2[3];
        for (int k = 0; k < 3; k++) v2[k] = ks->vgkc(k, j2);
        for (int j1 = 0; j1 <= j1_last; j1++) // for each column loop over rows
        {
            int ig = idxG12(ks->idxg[j1], ks->idxg[j2]);
            double t1 = 0.5 * (ks->vgkc(0, j1) * v2[0] + 
                               ks->vgkc(1, j1) * v2[1] + 
                               ks->vgkc(2, j1) * v2[2]);
            if (sblock == nm)
                h(j1, j2) += (lapw_runtime.veffig(ig) + t1 * lapw_global.cfunig[ig]);
            
            if (sblock == uu)
                h(j1, j2) += (lapw_runtime.veffig(ig) + t1 * lapw_global.cfunig[ig] + lapw_runtime.beffig(ig, 0));
            
            if (sblock == dd)
                h(j1, j2) += (lapw_runtime.veffig(ig) + t1 * lapw_global.cfunig[ig] - lapw_runtime.beffig(ig, 0));
            
            if (sblock == ud)
                h(j1, j2) += (lapw_runtime.beffig(ig, 1) - zi * lapw_runtime.beffig(ig, 2));
        }
    }
}

template <implementation impl> 
void lapw_set_o(bloch_states_k* const ks, mdarray<complex16,2>& o)
{
    timer t("lapw_set_o");

    if (impl == cpu)
    {
        zgemm<cpu>(0, 2, ks->ngk, ks->ngk, lapw_global.size_wfmt_apw, zone, &ks->apwalm(0, 0), ks->apwalm.size(0), 
            &ks->apwalm(0, 0), ks->apwalm.size(0), zzero, &o(0, 0), o.size(0));
    }
    if (impl == gpu)
    {
        o.allocate_on_device();
        o.zero_on_device();
        zgemm<gpu>(0, 2, ks->ngk, ks->ngk, lapw_global.size_wfmt_apw, zone, ks->apwalm.get_ptr_device(), ks->apwalm.size(0),
                   ks->apwalm.get_ptr_device(), ks->apwalm.size(0), zzero, o.get_ptr_device(), o.size(0));
        o.copy_to_host();
        o.deallocate_on_device();
        ks->apwalm.deallocate_on_device();
    }
    
    for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
    {
        Atom *atom = lapw_global.atoms[ias];
        Species *species = atom->species;
        
        int lo_index_offset = species->index.apw_size();

        for (int j2 = 0; j2 < species->index.lo_size(); j2++) // loop over columns (local-orbital block) 
        {
            int l2 = species->index[lo_index_offset + j2].l;
            int lm2 = species->index[lo_index_offset + j2].lm;
            int order2 = species->index[lo_index_offset + j2].order;
            
            // apw-lo block 
            for (int io1 = 0; io1 < (int)species->apw_descriptors[l2].radial_solution_descriptors.size(); io1++)
                for (int ig = 0; ig < ks->ngk; ig++)
                    o(ig, ks->ngk + atom->offset_lo + j2) += lapw_runtime.ovlprad(l2, io1, order2, ias) * ks->apwalm(ig, atom->offset_apw + species->index(lm2, io1)); 

            // lo-lo block
            for (int j1 = 0; j1 <= j2; j1++)
            {
                int lm1 = species->index[lo_index_offset + j1].lm;
                int order1 = species->index[lo_index_offset + j1].order;
                if (lm1 == lm2) 
                    o(ks->ngk + atom->offset_lo + j1, ks->ngk + atom->offset_lo + j2) += lapw_runtime.ovlprad(l2, order1, order2, ias);
            }
        }
    }
    
    for (int j2 = 0; j2 < ks->ngk; j2++) // loop over columns
        for (int j1 = 0; j1 <= j2; j1++) // for each column loop over rows
            o(j1, j2) += lapw_global.cfunig[idxG12(ks->idxg[j1], ks->idxg[j2])];
}

template <implementation impl, diagonalization mode> 
void lapw_band(bloch_states_k* const ks)
{
    timer t("lapw_band");

    if (mode == second_variational)
    {
        mdarray<complex16,2> h(NULL, ks->lapw_basis_size, ks->lapw_basis_size);
        h.allocate();
        h.zero();
        mdarray<complex16,2> o(NULL, ks->lapw_basis_size, ks->lapw_basis_size);
        o.allocate();
        o.zero();
  
        lapw_set_h<impl,nm>(ks, h);
        lapw_set_o<impl>(ks, o);
        
        ks->evecfv.set_dimensions(ks->lapw_basis_size, lapw_global.nstfv);
        ks->evecfv.allocate();
        ks->evalfv.resize(lapw_global.nstfv);
       
        timer *t1 = new timer("lapw_band:zhegv<impl>");
        zhegv<impl>(ks->lapw_basis_size, lapw_global.nstfv, lapw_global.evaltol, &h(0, 0), &o(0, 0), &ks->evalfv[0], 
            &ks->evecfv(0, 0), ks->evecfv.size(0));
        delete t1;

        ks->generate_scalar_wave_functions();

        if (check_scalar_wf)
            for (int i = 0; i < 3; i++) ks->test_scalar_wave_functions(i);
        
        ks->evecsv.set_dimensions(lapw_global.nstsv, lapw_global.nstsv);
        ks->evecsv.allocate();
        ks->evecsv.zero();
        ks->evalsv.resize(lapw_global.nstsv);
 
        if (lapw_global.ndmag == 0)
        {
            for (int i = 0; i < lapw_global.nstfv; i++)
            {
                ks->evecsv(i, i) = zone;
                ks->evalsv[i] = ks->evalfv[i];
            }
        } 
        else
        {
            lapw_set_sv(ks);
            t1 = new timer("lapw_band:zheev<cpu>"); 
            if (lapw_global.ndmag == 1)
            {
                zheev<cpu>(lapw_global.nstfv, &ks->evecsv(0, 0), ks->evecsv.size(0), &ks->evalsv[0]);
                zheev<cpu>(lapw_global.nstfv, &ks->evecsv(lapw_global.nstfv, lapw_global.nstfv), ks->evecsv.size(0), &ks->evalsv[lapw_global.nstfv]);
            } 
            if (lapw_global.ndmag == 3)
            {
                zheev<cpu>(lapw_global.nstsv, &ks->evecsv(0, 0), ks->evecsv.size(0), &ks->evalsv[0]);
            } 
            delete t1;
        }
    }

    if (mode == full)
    {
        ks->evecfd.set_dimensions(ks->lapw_basis_size * lapw_global.nspinor, lapw_global.nstsv);
        ks->evecfd.allocate();
        ks->evalsv.resize(lapw_global.nstsv);
        
        mdarray<complex16,2> o(NULL, ks->lapw_basis_size, ks->lapw_basis_size);
        o.allocate();
        o.zero();
        lapw_set_o<impl>(ks, o);

        if (lapw_global.ndmag == 0)
        {
            mdarray<complex16,2> h(NULL, ks->lapw_basis_size, ks->lapw_basis_size);
            h.allocate();
            h.zero();
            lapw_set_h<impl,nm>(ks, h);
            
            timer *t1 = new timer("lapw_band:zhegv<impl>");
            zhegv<impl>(ks->lapw_basis_size, lapw_global.nstfv, lapw_global.evaltol, &h(0, 0), &o(0, 0), &ks->evalsv[0],
                        &ks->evecfd(0, 0), ks->evecfd.size(0));
            delete t1;
        }

        if (lapw_global.ndmag == 1)
        {
            mdarray<complex16,2> o1(NULL, ks->lapw_basis_size, ks->lapw_basis_size);
            o1.allocate();
            memcpy(&o1(0, 0), &o(0, 0), o.size() * sizeof(complex16));
            
            mdarray<complex16,2> h(NULL, ks->lapw_basis_size, ks->lapw_basis_size);
            h.allocate();
            h.zero();
            lapw_set_h<impl,uu>(ks, h);
            
            timer *t1 = new timer("lapw_band:zhegv<impl>");
            zhegv<impl>(ks->lapw_basis_size, lapw_global.nstfv, lapw_global.evaltol, &h(0, 0), &o(0, 0), &ks->evalsv[0],
                        &ks->evecfd(0, 0), ks->evecfd.size(0));
            delete t1;
            
            h.zero();
            lapw_set_h<impl,dd>(ks, h);

            t1 = new timer("lapw_band:zhegv<impl>");
            zhegv<impl>(ks->lapw_basis_size, lapw_global.nstfv, lapw_global.evaltol, &h(0, 0), &o1(0, 0), &ks->evalsv[lapw_global.nstfv],
                        &ks->evecfd(ks->lapw_basis_size, lapw_global.nstfv), ks->evecfd.size(0));
            delete t1;
        }
        
        if (lapw_global.ndmag == 3)
        {
            mdarray<complex16,2> o1(NULL, ks->lapw_basis_size * lapw_global.nspinor, ks->lapw_basis_size * lapw_global.nspinor);
            o1.allocate();
            o1.zero();
            for (int i = 0; i < ks->lapw_basis_size; i++)
            {
                memcpy(&o1(0, i), &o(0, i), ks->lapw_basis_size * sizeof(complex16));
                memcpy(&o1(ks->lapw_basis_size, ks->lapw_basis_size + i), &o(0, i), ks->lapw_basis_size * sizeof(complex16));
            }

            mdarray<complex16,2> h(NULL, ks->lapw_basis_size * lapw_global.nspinor, ks->lapw_basis_size * lapw_global.nspinor);
            h.allocate();
            h.zero();

            mdarray<complex16,2> huu(&h(0, 0), ks->lapw_basis_size * lapw_global.nspinor, ks->lapw_basis_size);
            mdarray<complex16,2> hdd(&h(ks->lapw_basis_size, ks->lapw_basis_size), ks->lapw_basis_size * lapw_global.nspinor, ks->lapw_basis_size);
            mdarray<complex16,2> hud(&h(0, ks->lapw_basis_size), ks->lapw_basis_size * lapw_global.nspinor, ks->lapw_basis_size);

            lapw_set_h<impl,uu>(ks, huu);
            lapw_set_h<impl,ud>(ks, hud);
            lapw_set_h<impl,dd>(ks, hdd);
            
            timer *t1 = new timer("lapw_band:zhegv<impl>");
            zhegv<impl>(ks->lapw_basis_size * lapw_global.nspinor, lapw_global.nstsv, lapw_global.evaltol, &h(0, 0), &o1(0, 0), &ks->evalsv[0],
                        &ks->evecfd(0, 0), ks->evecfd.size(0));
            delete t1;
        }
    }
    
    ks->generate_spinor_wave_functions(mode);
    
    if (check_spinor_wf)
        for (int i = 0; i < 3; i++) ks->test_spinor_wave_functions(i);
}

#endif // __LAPW_BAND_H__
