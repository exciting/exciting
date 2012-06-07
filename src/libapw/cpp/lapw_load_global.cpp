#include "lapw.h"

extern "C" void FORTRAN(lapw_load_global)(int *lmaxvr_,
                                          int *lmaxapw_,
                                          int *apwordmax_,
                                          int *nrmtmax_,
                                          int *ngkmax_,
                                          int *ngvec_,
                                          int *ngrtot_,
                                          int *intgv_,
                                          int *ivg_,
                                          int *ivgig_,
                                          int *ngrid_,
                                          int *igfft_,
                                          double *cfunir_,
                                          complex16 *cfunig_,
                                          complex16 *gntyry_,
                                          int *nstfv_,
                                          int *nrfmtmax_,
                                          int *ordrfmtmax_,
                                          double *evaltol_,
                                          int *spinpol_,
                                          int *spinorb_,
                                          int *ndmag_,
                                          double *omega_,
                                          int *natmcls_,
                                          int *ic2ias_,
                                          int *natoms_in_class_,
                                          int *ldapu_)
{
    lapw_global.lmaxvr = *lmaxvr_;
    lapw_global.lmmaxvr = pow(lapw_global.lmaxvr + 1, 2);
    lapw_global.lmaxapw = *lmaxapw_;
    lapw_global.lmmaxapw = pow(lapw_global.lmaxapw + 1, 2);
    lapw_global.apwordmax = *apwordmax_;
    lapw_global.nrmtmax = *nrmtmax_;
    lapw_global.ngkmax = *ngkmax_;
    lapw_global.ngvec = *ngvec_;
    lapw_global.ngrtot = *ngrtot_;
    lapw_global.nstfv = *nstfv_;
    lapw_global.nrfmtmax = *nrfmtmax_;
    lapw_global.ordrfmtmax = *ordrfmtmax_;
    lapw_global.evaltol = *evaltol_;
    
    lapw_global.lmaxlu = 3;
    lapw_global.lmmaxlu = pow(lapw_global.lmaxlu + 1, 2);
    
    lapw_global.intgv.set_dimensions(3, 2);
    lapw_global.intgv.set_ptr(intgv_);
    
    lapw_global.ivgig.set_dimensions(dimension(lapw_global.intgv(0, 0), lapw_global.intgv(0, 1)),
                                     dimension(lapw_global.intgv(1, 0), lapw_global.intgv(1, 1)),
                                     dimension(lapw_global.intgv(2, 0), lapw_global.intgv(2, 1)));
    lapw_global.ivgig.allocate();  
    mdarray<int,3> ivgig_tmp(ivgig_, dimension(lapw_global.intgv(0, 0), lapw_global.intgv(0, 1)),
                                     dimension(lapw_global.intgv(1, 0), lapw_global.intgv(1, 1)),
                                     dimension(lapw_global.intgv(2, 0), lapw_global.intgv(2, 1)));
    
    for (int i = lapw_global.intgv(0, 0); i <= lapw_global.intgv(0, 1); i++)
        for (int j = lapw_global.intgv(1, 0); j <= lapw_global.intgv(1, 1); j++)
            for (int k = lapw_global.intgv(2, 0); k <= lapw_global.intgv(2, 1); k++)
                lapw_global.ivgig(i, j, k) = ivgig_tmp(i, j, k) - 1;
    
    lapw_global.ivg.set_dimensions(3, lapw_global.ngrtot);
    lapw_global.ivg.set_ptr(ivg_);
    
    lapw_global.igfft.resize(lapw_global.ngrtot);
    lapw_global.cfunir.resize(lapw_global.ngrtot);
    lapw_global.cfunig.resize(lapw_global.ngrtot);
    for (int i = 0; i < lapw_global.ngrtot; i++)
    {
        lapw_global.cfunig[i] = cfunig_[i];
        lapw_global.cfunir[i] = cfunir_[i];
        lapw_global.igfft[i] = igfft_[i] - 1;
    }
    lapw_global.ngrid[0] = ngrid_[0];
    lapw_global.ngrid[1] = ngrid_[1];
    lapw_global.ngrid[2] = ngrid_[2];
    
    lapw_global.spinpol = (*spinpol_ != 0);
    lapw_global.ndmag = *ndmag_;
    lapw_global.nspinor = (lapw_global.spinpol) ? 2 : 1;
    lapw_global.nstsv = lapw_global.nstfv * lapw_global.nspinor;
    
    lapw_global.spinorb = (*spinorb_ != 0);
    lapw_global.ldapu = (*ldapu_ != 0);

    lapw_global.natmcls = *natmcls_;
    lapw_global.ic2ias.resize(lapw_global.natmcls);
    lapw_global.natoms_in_class.resize(lapw_global.natmcls);
    for (int ic = 0; ic < (int)lapw_global.ic2ias.size(); ic++)
    {
        lapw_global.ic2ias[ic] = ic2ias_[ic] - 1;
        lapw_global.natoms_in_class[ic] = natoms_in_class_[ic];
    }
    
    lapw_global.gntyry.set_dimensions(lapw_global.lmmaxvr, lapw_global.lmmaxapw, lapw_global.lmmaxapw);
    lapw_global.gntyry.set_ptr(gntyry_);
    
    lapw_global.L3_gntyry.set_dimensions(lapw_global.lmmaxapw, lapw_global.lmmaxapw);
    lapw_global.L3_gntyry.allocate();
    
    for (int lm1 = 0; lm1 < lapw_global.lmmaxapw; lm1++)
        for (int lm2 = 0; lm2 < lapw_global.lmmaxapw; lm2++)
            for (int lm3 = 0; lm3 < lapw_global.lmmaxvr; lm3++) 
                if (abs(lapw_global.gntyry(lm3, lm1, lm2)) > 1e-14)
                {
                    lapw_global.L3_gntyry(lm1, lm2).push_back(lm3);
                }
    
    
    lapw_global.omega = *omega_;
    
    for (int i = 0; i < (int)lapw_global.species.size(); i++)
        delete lapw_global.species[i];
    lapw_global.species.clear();

    for (int i = 0; i < (int)lapw_global.atoms.size(); i++)
        delete lapw_global.atoms[i];
    lapw_global.atoms.clear();

    for (int i = 0; i < (int)lapw_runtime.bloch_states.size(); i++)
        delete lapw_runtime.bloch_states[i];
    lapw_runtime.bloch_states.clear();
}


