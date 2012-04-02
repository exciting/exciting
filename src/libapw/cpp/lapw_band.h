#ifndef _LAPW_GEV_H_
#define _LAPW_GEV_H_

#include "lapw.h"

template <spin_block sblock> 
void hmt_dot_apw(lapw_eigen_states& eigen_states, mdarray<complex16,2>& hapw)
{
    timer t("hmt_dot_apw");
   
    unsigned int ngk = eigen_states.kp->ngk;

    #pragma omp parallel default(shared)
    {
        std::vector<complex16> zv(ngk);
        std::vector<double> v1(p.lmmaxvr);
        std::vector<complex16> v2(p.lmmaxvr);
        #pragma omp for
        for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
        {
            Atom *atom = &geometry.atoms[ias];
            Species *species = atom->species;
        
            // precompute apw block
            for (unsigned int j2 = 0; j2 < species->size_ci_apw; j2++)
            {
                memset(&zv[0], 0, ngk * sizeof(complex16));
                
                int lm2 = species->ci[j2].lm;
                int idxrf2 = species->ci[j2].idxrf;
                
                for (unsigned int j1 = 0; j1 < species->size_ci_apw; j1++)
                {
                    int lm1 = species->ci[j1].lm;
                    int idxrf1 = species->ci[j1].idxrf;
                    
                    complex16 zsum(0, 0);
                    
                    if (sblock == nm)
                    {
                        L3_sum_gntyry(lm1, lm2, &p.hmltrad(0, idxrf1, idxrf2, ias), zsum);
                    }

                    if (sblock == uu)
                    {
                        for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++) 
                            v1[lm3] = p.hmltrad(lm3, idxrf1, idxrf2, ias) + p.beffrad(lm3, idxrf1, idxrf2, ias, 0);
                        L3_sum_gntyry(lm1, lm2, &v1[0], zsum);
                    }
                    
                    if (sblock == dd)
                    {
                        for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++) 
                            v1[lm3] = p.hmltrad(lm3, idxrf1, idxrf2, ias) - p.beffrad(lm3, idxrf1, idxrf2, ias, 0);
                        L3_sum_gntyry(lm1, lm2, &v1[0], zsum);
                    }
                    
                    if (sblock == ud)
                    {
                        for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++) 
                            v2[lm3] = complex16(p.beffrad(lm3, idxrf1, idxrf2, ias, 1), -p.beffrad(lm3, idxrf1, idxrf2, ias, 2));
                        L3_sum_gntyry(lm1, lm2, &v2[0], zsum);
                    }
        
                    if (abs(zsum) > 1e-14) 
                        for (unsigned int ig = 0; ig < ngk; ig++) 
                            zv[ig] += zsum * eigen_states.apwalm(ig, atom->offset_apw + j1); 
                }
                
                // surface term
                if (sblock != ud)
                {
                    int l2 = species->ci[j2].l;
                    int io2 = species->ci[j2].order;
                    
                    for (unsigned int io1 = 0; io1 < species->apw_descriptors[l2].radial_solution_descriptors.size(); io1++)
                    {
                        double t1 = 0.5 * pow(species->rmt, 2) * p.apwfr(species->nrmt - 1, 0, io1, l2, ias) * p.apwdfr(io2, l2, ias); 
                        for (unsigned int ig = 0; ig < ngk; ig++) 
                            zv[ig] += t1 * eigen_states.apwalm(ig, atom->offset_apw + species->ci_by_lmo(lm2, io1));
                    }
                }

                memcpy(&hapw(0, atom->offset_apw + j2), &zv[0], ngk * sizeof(complex16));
            }
        } 
    }
}


template <implementation impl, spin_block sblock> 
void lapw_set_h(lapw_eigen_states& eigen_states, mdarray<complex16,2>& h)
{
    timer t("lapw_set_h");
   
    unsigned int ngk = eigen_states.kp->ngk;
    
    mdarray<complex16,2> zm(NULL, ngk, p.size_wfmt_apw);
    zm.allocate();
    
    hmt_dot_apw<sblock>(eigen_states, zm);

    if (impl == cpu)
    {
        zgemm<cpu>(0, 2, ngk, ngk, p.size_wfmt_apw, zone, &zm(0, 0), zm.size(0), 
                   &eigen_states.apwalm(0, 0), eigen_states.apwalm.size(0), zzero, 
                   &h(0, 0), h.size(0));
    }
    if (impl == gpu)
    {
        eigen_states.apwalm.allocate_on_device();
        eigen_states.apwalm.copy_to_device();
        zm.allocate_on_device();
        zm.copy_to_device();
        h.allocate_on_device();
        h.zero_on_device();
        zgemm<gpu>(0, 2, ngk, ngk, p.size_wfmt_apw, zone, zm.get_ptr_device(), zm.size(0),
                   eigen_states.apwalm.get_ptr_device(), eigen_states.apwalm.size(0), zzero, 
                   h.get_ptr_device(), h.size(0));
        h.copy_to_host();
        h.deallocate_on_device();
        zm.deallocate_on_device();
    }
 
    #pragma omp parallel default(shared)
    {
        std::vector<double> v1(p.lmmaxvr);
        std::vector<complex16> v2(p.lmmaxvr);
        #pragma omp for
        for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
        {
            Atom *atom = &geometry.atoms[ias];
            Species *species = atom->species;
    
            for (unsigned int j2 = 0; j2 < species->size_ci_lo; j2++) // loop over columns (local-orbital block) 
            {
                int lm2 = species->ci_lo[j2].lm;
                int idxrf2 = species->ci_lo[j2].idxrf;
                
                // apw-lo block
                for (unsigned int j1 = 0; j1 < species->size_ci_apw; j1++) // loop over rows
                {
                    int lm1 = species->ci[j1].lm;
                    int idxrf1 = species->ci[j1].idxrf;
                    
                    complex16 zsum(0, 0);
                    
                    if (sblock == nm)
                    {
                        L3_sum_gntyry(lm1, lm2, &p.hmltrad(0, idxrf2, idxrf1, ias), zsum);
                    }

                    if (sblock == uu)
                    {
                        for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++) 
                            v1[lm3] = p.hmltrad(lm3, idxrf2, idxrf1, ias) + p.beffrad(lm3, idxrf1, idxrf2, ias, 0);
                        L3_sum_gntyry(lm1, lm2, &v1[0], zsum);
                    }
                    
                    if (sblock == dd)
                    {
                        for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++) 
                            v1[lm3] = p.hmltrad(lm3, idxrf2, idxrf1, ias) - p.beffrad(lm3, idxrf1, idxrf2, ias, 0);
                        L3_sum_gntyry(lm1, lm2, &v1[0], zsum);
                    }
                    
                    if (sblock == ud)
                    {
                        for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++) 
                            v2[lm3] = complex16(p.beffrad(lm3, idxrf1, idxrf2, ias, 1), -p.beffrad(lm3, idxrf1, idxrf2, ias, 2));
                        L3_sum_gntyry(lm1, lm2, &v2[0], zsum);
                    }
        
                    if (abs(zsum) > 1e-14)
                        for (unsigned int ig = 0; ig < ngk; ig++)
                            h(ig, ngk + atom->offset_lo + j2) += zsum * eigen_states.apwalm(ig, atom->offset_apw + j1);
                }
                
                unsigned int j1_end;
                if (sblock == ud) j1_end = species->size_ci_lo - 1;
                else j1_end = j2;
    
                // lo-lo block 
                for (unsigned int j1 = 0; j1 <= j1_end; j1++)
                {
                    int lm1 = species->ci_lo[j1].lm;
                    int idxrf1 = species->ci_lo[j1].idxrf;
                    
                    complex16 zsum(0, 0);

                    if (sblock == nm)
                    {
                        L3_sum_gntyry(lm1, lm2, &p.hmltrad(0, idxrf1, idxrf2, ias), zsum);
                    }

                    if (sblock == uu)
                    {
                        for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++) 
                            v1[lm3] = p.hmltrad(lm3, idxrf1, idxrf2, ias) + p.beffrad(lm3, idxrf1, idxrf2, ias, 0);
                        L3_sum_gntyry(lm1, lm2, &v1[0], zsum);
                    }
                    
                    if (sblock == dd)
                    {
                        for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++) 
                            v1[lm3] = p.hmltrad(lm3, idxrf1, idxrf2, ias) - p.beffrad(lm3, idxrf1, idxrf2, ias, 0);
                        L3_sum_gntyry(lm1, lm2, &v1[0], zsum);
                    }
                    
                    if (sblock == ud)
                    {
                        for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++) 
                            v2[lm3] = complex16(p.beffrad(lm3, idxrf1, idxrf2, ias, 1), -p.beffrad(lm3, idxrf1, idxrf2, ias, 2));
                        L3_sum_gntyry(lm1, lm2, &v2[0], zsum);
                    }

                    h(ngk + atom->offset_lo + j1, ngk + atom->offset_lo + j2) += zsum;
                }
            }

            if (sblock == ud)
            {
                for (unsigned int j2 = 0; j2 < species->size_ci_apw; j2++)
                {
                    int lm2 = species->ci[j2].lm;
                    int idxrf2 = species->ci[j2].idxrf;
                    
                    for (unsigned int j1 = 0; j1 < species->size_ci_lo; j1++)
                    {
                        int lm1 = species->ci_lo[j1].lm;
                        int idxrf1 = species->ci_lo[j1].idxrf;
                        
                        complex16 zsum(0, 0);
                        
                        for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++) 
                            v2[lm3] = complex16(p.beffrad(lm3, idxrf1, idxrf2, ias, 1), -p.beffrad(lm3, idxrf1, idxrf2, ias, 2));
                        L3_sum_gntyry(lm1, lm2, &v2[0], zsum);
                        
                        if (abs(zsum) > 1e-14)
                            for (unsigned int ig = 0; ig < ngk; ig++)
                                h(ngk + atom->offset_lo + j1, ig) += zsum * conj(eigen_states.apwalm(ig, atom->offset_apw + j2));
                    }
                }
            }
        } 
    }
    
    for (unsigned int j2 = 0; j2 < ngk; j2++) // loop over columns
    {
        unsigned int j1_end;
        if (sblock == ud) j1_end = ngk - 1;
        else j1_end = j2;
   
        double v2[3];
        for (int k = 0; k < 3; k++) v2[k] = eigen_states.kp->vgkc(k, j2);
        for (unsigned int j1 = 0; j1 <= j1_end; j1++) // for each column loop over rows
        {
            int ig = idxG12(eigen_states.kp, j1, j2);
            double t1 = 0.5 * (eigen_states.kp->vgkc(0, j1) * v2[0] + 
                               eigen_states.kp->vgkc(1, j1) * v2[1] + 
                               eigen_states.kp->vgkc(2, j1) * v2[2]);
            if (sblock == nm)
                h(j1, j2) += (p.veffig(ig) + t1 * p.cfunig[ig]);
            
            if (sblock == uu)
                h(j1, j2) += (p.veffig(ig) + t1 * p.cfunig[ig] + p.beffig(ig, 0));
            
            if (sblock == dd)
                h(j1, j2) += (p.veffig(ig) + t1 * p.cfunig[ig] - p.beffig(ig, 0));
            
            if (sblock == ud)
                h(j1, j2) += (p.beffig(ig, 1) - zi * p.beffig(ig, 2));
        }
    }
}

template <implementation impl> 
void lapw_set_o(lapw_eigen_states& eigen_states, mdarray<complex16,2>& o)
{
    timer t("lapw_set_o");

    unsigned int ngk = eigen_states.kp->ngk;
    
    if (impl == cpu)
    {
        zgemm<cpu>(0, 2, ngk, ngk, p.size_wfmt_apw, zone, &eigen_states.apwalm(0, 0), eigen_states.apwalm.size(0), 
            &eigen_states.apwalm(0, 0), eigen_states.apwalm.size(0), zzero, &o(0, 0), o.size(0));
    }
    if (impl == gpu)
    {
        o.allocate_on_device();
        o.zero_on_device();
        zgemm<gpu>(0, 2, ngk, ngk, p.size_wfmt_apw, zone, eigen_states.apwalm.get_ptr_device(), eigen_states.apwalm.size(0),
                   eigen_states.apwalm.get_ptr_device(), eigen_states.apwalm.size(0), zzero, o.get_ptr_device(), o.size(0));
        o.copy_to_host();
        o.deallocate_on_device();
        eigen_states.apwalm.deallocate_on_device();
    }
    
    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        Species *species = atom->species;

        for (unsigned int j2 = 0; j2 < species->size_ci_lo; j2++) // loop over columns (local-orbital block) 
        {
            int l2 = species->ci_lo[j2].l;
            int lm2 = species->ci_lo[j2].lm;
            int order2 = species->ci_lo[j2].order;
            
            // apw-lo block 
            for (unsigned int io1 = 0; io1 < species->apw_descriptors[l2].radial_solution_descriptors.size(); io1++)
                for (unsigned int ig = 0; ig < ngk; ig++)
                    o(ig, ngk + atom->offset_lo + j2) += p.ovlprad(l2, io1, order2, ias) * eigen_states.apwalm(ig, atom->offset_apw + species->ci_by_lmo(lm2, io1)); 

            // lo-lo block
            for (unsigned int j1 = 0; j1 <= j2; j1++)
            {
                int lm1 = species->ci_lo[j1].lm;
                int order1 = species->ci_lo[j1].order;
                if (lm1 == lm2) 
                    o(ngk + atom->offset_lo + j1, ngk + atom->offset_lo + j2) += p.ovlprad(l2, order1, order2, ias);
            }
        }
    }
    
    for (unsigned int j2 = 0; j2 < ngk; j2++) // loop over columns
        for (unsigned int j1 = 0; j1 <= j2; j1++) // for each column loop over rows
            o(j1, j2) += p.cfunig[idxG12(eigen_states.kp, j1, j2)];
}

template <implementation impl, diagonalization mode> 
void lapw_band(lapw_eigen_states& eigen_states)
{
    timer t("lapw_band");

    unsigned int msize = eigen_states.kp->ngk + p.size_wfmt_lo;
    
    if (mode == second_variational)
    {
        mdarray<complex16,2> h(NULL, msize, msize);
        h.allocate();
        h.zero();
        mdarray<complex16,2> o(NULL, msize, msize);
        o.allocate();
        o.zero();
  
        lapw_set_h<impl,nm>(eigen_states, h);
        lapw_set_o<impl>(eigen_states, o);
        
        eigen_states.evecfv.set_dimensions(msize, p.nstfv);
        eigen_states.evecfv.allocate();
        eigen_states.evalfv.resize(p.nstfv);
       
        timer *t1 = new timer("lapw_band:zhegv<impl>");
        zhegv<impl>(msize, p.nstfv, p.evaltol, &h(0, 0), &o(0, 0), &eigen_states.evalfv[0], 
            &eigen_states.evecfv(0, 0), eigen_states.evecfv.size(0));
        delete t1;

        eigen_states.generate_scalar_wf();

        if (check_scalar_wf)
            for (int i = 0; i < 3; i++) eigen_states.test_scalar_wf(i);
        
        eigen_states.evecsv.set_dimensions(p.nstsv, p.nstsv);
        eigen_states.evecsv.allocate();
        eigen_states.evecsv.zero();
        eigen_states.evalsv.resize(p.nstsv);
 
        if (p.ndmag == 0)
        {
            for (unsigned int i = 0; i < p.nstfv; i++)
            {
                eigen_states.evecsv(i, i) = zone;
                eigen_states.evalsv[i] = eigen_states.evalfv[i];
            }
        } 
        else
        {
            lapw_set_sv(eigen_states);
            t1 = new timer("lapw_band:zheev<cpu>"); 
            if (p.ndmag == 1)
            {
                zheev<cpu>(p.nstfv, &eigen_states.evecsv(0, 0), eigen_states.evecsv.size(0), 
                    &eigen_states.evalsv[0]);
                zheev<cpu>(p.nstfv, &eigen_states.evecsv(p.nstfv, p.nstfv), eigen_states.evecsv.size(0), 
                    &eigen_states.evalsv[p.nstfv]);
            } 
            if (p.ndmag == 3)
            {
                zheev<cpu>(p.nstsv, &eigen_states.evecsv(0, 0), eigen_states.evecsv.size(0), 
                    &eigen_states.evalsv[0]);
            } 
            delete t1;
        }
    }

    if (mode == full)
    {
        eigen_states.evecfd.set_dimensions(msize * p.nspinor, p.nstsv);
        eigen_states.evecfd.allocate();
        eigen_states.evalsv.resize(p.nstsv);
        
        mdarray<complex16,2> o(NULL, msize, msize);
        o.allocate();
        o.zero();
        lapw_set_o<impl>(eigen_states, o);

        if (p.ndmag == 0)
        {
            mdarray<complex16,2> h(NULL, msize, msize);
            h.allocate();
            h.zero();
            lapw_set_h<impl,nm>(eigen_states, h);
            
            timer *t1 = new timer("lapw_band:zhegv<impl>");
            zhegv<impl>(msize, p.nstfv, p.evaltol, &h(0, 0), &o(0, 0), &eigen_states.evalsv[0],
                        &eigen_states.evecfd(0, 0), eigen_states.evecfd.size(0));
            delete t1;
        }

        if (p.ndmag == 1)
        {
            mdarray<complex16,2> o1(NULL, msize, msize);
            o1.allocate();
            memcpy(&o1(0, 0), &o(0, 0), o.size() * sizeof(complex16));
            
            mdarray<complex16,2> h(NULL, msize, msize);
            h.allocate();
            h.zero();
            lapw_set_h<impl,uu>(eigen_states, h);
            
            timer *t1 = new timer("lapw_band:zhegv<impl>");
            zhegv<impl>(msize, p.nstfv, p.evaltol, &h(0, 0), &o(0, 0), &eigen_states.evalsv[0],
                        &eigen_states.evecfd(0, 0), eigen_states.evecfd.size(0));
            delete t1;
            
            h.zero();
            lapw_set_h<impl,dd>(eigen_states, h);

            t1 = new timer("lapw_band:zhegv<impl>");
            zhegv<impl>(msize, p.nstfv, p.evaltol, &h(0, 0), &o1(0, 0), &eigen_states.evalsv[p.nstfv],
                        &eigen_states.evecfd(msize, p.nstfv), eigen_states.evecfd.size(0));
            delete t1;
        }
        
        if (p.ndmag == 3)
        {
            mdarray<complex16,2> o1(NULL, msize * p.nspinor, msize * p.nspinor);
            o1.allocate();
            o1.zero();
            for (unsigned int i = 0; i < msize; i++)
            {
                memcpy(&o1(0, i), &o(0, i), msize * sizeof(complex16));
                memcpy(&o1(msize, msize + i), &o(0, i), msize * sizeof(complex16));
            }

            mdarray<complex16,2> h(NULL, msize * p.nspinor, msize * p.nspinor);
            h.allocate();
            h.zero();

            mdarray<complex16,2> huu(&h(0, 0), msize * p.nspinor, msize);
            mdarray<complex16,2> hdd(&h(msize, msize), msize * p.nspinor, msize);
            mdarray<complex16,2> hud(&h(0, msize), msize * p.nspinor, msize);

            lapw_set_h<impl,uu>(eigen_states, huu);
            lapw_set_h<impl,ud>(eigen_states, hud);
            lapw_set_h<impl,dd>(eigen_states, hdd);
            
            timer *t1 = new timer("lapw_band:zhegv<impl>");
            zhegv<impl>(msize * p.nspinor, p.nstsv, p.evaltol, &h(0, 0), &o1(0, 0), &eigen_states.evalsv[0],
                        &eigen_states.evecfd(0, 0), eigen_states.evecfd.size(0));
            delete t1;
        }
    }
    
    eigen_states.generate_spinor_wf(mode);
    
    if (check_spinor_wf)
        for (int i = 0; i < 3; i++) eigen_states.test_spinor_wf(i);
}

#endif // _LAPW_GEV_H_
