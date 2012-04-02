#include "lapw.h"

extern "C" void FORTRAN(lapw_load_global)(int *natmtot_,
                                          int *nspecies_,
                                          int *lmaxvr_,
                                          int *lmaxapw_,
                                          int *apwordmax_,
                                          int *nrmtmax_,
                                          int *ngkmax_,
                                          int *ngvec_,
                                          int *ngrtot_,
                                          int *nlomax_,
                                          int *ias2is_,
                                          int *intgv_,
                                          int *ivg_,
                                          int *ivgig_,
                                          int *ngrid_,
                                          int *igfft_,
                                          double *cfunir_,
                                          std::complex<double> *cfunig_,
                                          std::complex<double> *gntyry_,
                                          int *nstfv_,
                                          int *nstsv_,
                                          int *nmatmax_,
                                          int *nrfmtmax_,
                                          int *ordrfmtmax_,
                                          double *evaltol_,
                                          int *spinpol_,
                                          int *ndmag_,
                                          double *omega_,
                                          int *natmcls_,
                                          int *ic2ias_,
                                          int *natoms_in_class_)
{
    p.natmtot = *natmtot_;
    p.nspecies = *nspecies_;
    p.lmaxvr = *lmaxvr_;
    p.lmmaxvr = pow(p.lmaxvr + 1, 2);
    p.lmaxapw = *lmaxapw_;
    p.lmmaxapw = pow(p.lmaxapw + 1, 2);
    p.apwordmax = *apwordmax_;
    p.nrmtmax = *nrmtmax_;
    p.ngkmax = *ngkmax_;
    p.ngvec = *ngvec_;
    p.ngrtot = *ngrtot_;
    p.nlomax = *nlomax_;
    p.nstfv = *nstfv_;
    p.nstsv = *nstsv_;
    p.nmatmax = *nmatmax_;
    p.nrfmtmax = *nrfmtmax_;
    p.ordrfmtmax = *ordrfmtmax_;
    p.evaltol = *evaltol_;
    
    p.intgv.set_dimensions(3, 2);
    p.intgv.set_ptr(intgv_);
    
    p.ivgig.set_dimensions(dimension(p.intgv(0, 0), p.intgv(0, 1)),
                           dimension(p.intgv(1, 0), p.intgv(1, 1)),
                           dimension(p.intgv(2, 0), p.intgv(2, 1)));
    p.ivgig.allocate();  
    mdarray<int,3> ivgig_tmp(ivgig_, dimension(p.intgv(0, 0), p.intgv(0, 1)),
                                     dimension(p.intgv(1, 0), p.intgv(1, 1)),
                                     dimension(p.intgv(2, 0), p.intgv(2, 1)));
    
    p.ivg.set_dimensions(3, p.ngrtot);
    p.ivg.set_ptr(ivg_);
    
    p.igfft.resize(p.ngrtot);
    p.cfunir.resize(p.ngrtot);
    p.cfunig.resize(p.ngrtot);
    for (unsigned int i = 0; i < p.ngrtot; i++)
    {
        p.cfunig[i] = cfunig_[i];
        p.cfunir[i] = cfunir_[i];
        p.igfft[i] = igfft_[i] - 1;
    }
    p.ngrid[0] = ngrid_[0];
    p.ngrid[1] = ngrid_[1];
    p.ngrid[2] = ngrid_[2];
    
    p.spinpol = (*spinpol_ != 0);
    p.ndmag = *ndmag_;
    p.nspinor = (p.spinpol) ? 2 : 1;

    p.natmcls = *natmcls_;
    p.ic2ias.resize(p.natmcls);
    p.natoms_in_class.resize(p.natmcls);
    for (unsigned int ic = 0; ic < p.ic2ias.size(); ic++)
    {
        p.ic2ias[ic] = ic2ias_[ic] - 1;
        p.natoms_in_class[ic] = natoms_in_class_[ic];
    }
    
    p.gntyry.set_dimensions(p.lmmaxvr, p.lmmaxapw, p.lmmaxapw);
    p.gntyry.set_ptr(gntyry_);
    
    p.L3_gntyry.set_dimensions(p.lmmaxapw, p.lmmaxapw);
    p.L3_gntyry.allocate();
    
    p.L3_gntyry_data.set_dimensions(p.lmmaxapw, p.lmmaxapw);
    p.L3_gntyry_data.allocate();
    
    for (unsigned int lm1 = 0; lm1 < p.lmmaxapw; lm1++)
        for (unsigned int lm2 = 0; lm2 < p.lmmaxapw; lm2++)
            for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++) 
                if (abs(p.gntyry(lm3, lm1, lm2)) > 1e-14)
                {
                    p.L3_gntyry(lm1, lm2).push_back(lm3);
                    p.L3_gntyry_data(lm1, lm2).push_back(p.gntyry(lm3, lm1, lm2));
                }
    
    for (int i = p.intgv(0, 0); i <= p.intgv(0, 1); i++)
        for (int j = p.intgv(1, 0); j <= p.intgv(1, 1); j++)
            for (int k = p.intgv(2, 0); k <= p.intgv(2, 1); k++)
                p.ivgig(i, j, k) = ivgig_tmp(i, j, k) - 1;
    
    geometry.omega = *omega_;
    
    for (unsigned int i = 0; i < geometry.species.size(); i++)
        delete geometry.species[i];
    geometry.species.clear();
    for (unsigned int i = 0; i < p.nspecies; i++)
        geometry.species.push_back(new Species());
    
    geometry.atoms.clear();
    for (unsigned int i = 0; i < p.natmtot; i++)
        geometry.atoms.push_back(Atom(geometry.species[ias2is_[i] - 1]));

    for (unsigned int i = 0; i < p.kpoints.size(); i++)
        delete p.kpoints[i];
    p.kpoints.clear();
}


