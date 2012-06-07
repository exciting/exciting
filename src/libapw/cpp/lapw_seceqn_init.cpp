#include <cstdio>
#include "lapw.h"

extern "C" void FORTRAN(lapw_seceqn_init)(double *hmltrad_, double *ovlprad_, double *beffrad_,
                                          double *apwfr_, double *apwdfr_, double *beffir_,
                                          complex16 *veffig_, complex16* vmatu_, double *rfmt_, double *socrfmt_)
                                     
{
    lapw_runtime.hmltrad.set_dimensions(lapw_global.lmmaxvr, lapw_global.nrfmtmax, lapw_global.nrfmtmax, lapw_global.atoms.size());
    lapw_runtime.hmltrad.set_ptr(hmltrad_);
    
    lapw_runtime.rfmt.set_dimensions(lapw_global.nrmtmax, lapw_global.nrfmtmax, lapw_global.natmcls);
    lapw_runtime.rfmt.set_ptr(rfmt_);

    lapw_runtime.ovlprad.set_dimensions(lapw_global.lmaxapw + 1, lapw_global.ordrfmtmax, lapw_global.ordrfmtmax, lapw_global.atoms.size());
    lapw_runtime.ovlprad.set_ptr(ovlprad_);
    
    /*lapw_runtime.ovlprad.allocate();

    std::vector<double> f(lapw_global.nrmtmax);
    for (int ias = 0; ias < lapw_global.atoms.size(); ias++)
    {
        Species *sp = lapw_global.atoms[ias]->species;
        int ic = lapw_global.atoms[ias]->idxclass;

        for (int l = 0; l <= lapw_global.lmaxapw; l++)
        {
            int nrf = sp->idxmap.getnrf(l);
            for (int io1 = 0; io1 < nrf; io1++)
            {
                for (int io2 = 0; io2 < nrf; io2++)
                {
                    for (int ir = 0; ir < sp->nrmt; ir++)
                        f[ir] = lapw_runtime.rfmt(ir, sp->idxmap.getidxrf(l, io1), ic) * lapw_runtime.rfmt(ir, sp->idxmap.getidxrf(l, io2), ic) * 
                            pow(sp->radial_mesh(ir), 2);
                    lapw_runtime.ovlprad(l, io1, io2, ias) = lapw_spline_integrate(sp->nrmt, &sp->radial_mesh(0), &f[0]); 
                }
            }
        }
    }*/
    if (lapw_global.spinorb)
    {
        lapw_runtime.socrad.set_dimensions(lapw_global.lmaxapw + 1, lapw_global.ordrfmtmax, lapw_global.ordrfmtmax, lapw_global.atoms.size());
        lapw_runtime.socrad.allocate();

        mdarray<double,2> socrfmt(socrfmt_, lapw_global.nrmtmax, lapw_global.atoms.size());
        
        std::vector<double> f(lapw_global.nrmtmax);
        for (int ias = 0; ias < (int)lapw_global.atoms.size(); ias++)
        {
            Species *sp = lapw_global.atoms[ias]->species;
            int ic = lapw_global.atoms[ias]->idxclass;

            for (int l = 0; l <= lapw_global.lmaxapw; l++)
            {
                int nrf = sp->radial_index.nrf(l);
                for (int io1 = 0; io1 < nrf; io1++)
                {
                    for (int io2 = 0; io2 < nrf; io2++)
                    {
                        for (int ir = 0; ir < sp->nrmt; ir++)
                            f[ir] = lapw_runtime.rfmt(ir, sp->radial_index(l, io1), ic) * socrfmt(ir, ias) * 
                                lapw_runtime.rfmt(ir, sp->radial_index(l, io2), ic) * pow(sp->radial_mesh(ir), 2);
                        lapw_runtime.socrad(l, io1, io2, ias) = lapw_spline_integrate(sp->nrmt, &sp->radial_mesh(0), &f[0]); 
                    }
                }
            }
        }
    }
    
    
    
    lapw_runtime.apwfr.set_dimensions(lapw_global.nrmtmax, 2, lapw_global.apwordmax, lapw_global.lmaxapw + 1, lapw_global.atoms.size());
    lapw_runtime.apwfr.set_ptr(apwfr_);

    lapw_runtime.apwdfr.set_dimensions(lapw_global.apwordmax, lapw_global.lmaxapw + 1, lapw_global.atoms.size());
    lapw_runtime.apwdfr.set_ptr(apwdfr_);

    lapw_runtime.veffig.set_dimensions(lapw_global.ngvec);
    lapw_runtime.veffig.set_ptr(veffig_);

    if (lapw_global.ndmag > 0) 
    {
        lapw_runtime.beffrad.set_dimensions(lapw_global.lmmaxvr, lapw_global.nrfmtmax, lapw_global.nrfmtmax, lapw_global.atoms.size(), lapw_global.ndmag);
        lapw_runtime.beffrad.set_ptr(beffrad_);

        lapw_runtime.beffir.set_dimensions(lapw_global.ngrtot, lapw_global.ndmag);
        lapw_runtime.beffir.set_ptr(beffir_);
        
        if (lapw_diag == full)
        {
            lapw_runtime.beffig.set_dimensions(lapw_global.ngvec, lapw_global.ndmag);
            lapw_runtime.beffig.allocate();
            std::vector<complex16> zfft(lapw_global.ngrtot);
            for (int i = 0; i < lapw_global.ndmag; i++)
            {
                for (int ir = 0; ir < lapw_global.ngrtot; ir++) zfft[ir] = zone * lapw_runtime.beffir(ir, i) * lapw_global.cfunir[ir];
                lapw_fft(-1, &zfft[0]);
                for (int ig = 0; ig < lapw_global.ngvec; ig++) lapw_runtime.beffig(ig, i) = zfft[lapw_global.igfft[ig]];
            }
        }
    }
    
    if (lapw_global.ldapu)
    {
        lapw_runtime.dmatu.set_dimensions(lapw_global.lmmaxlu, lapw_global.lmmaxlu, lapw_global.nspinor, lapw_global.nspinor, lapw_global.atoms.size());
        lapw_runtime.dmatu.allocate();
        lapw_runtime.dmatu.zero();
        
        lapw_runtime.vmatu.set_ptr(vmatu_);
        lapw_runtime.vmatu.set_dimensions(lapw_global.lmmaxlu, lapw_global.lmmaxlu, lapw_global.nspinor, lapw_global.nspinor, lapw_global.atoms.size());
    }
    
    /*printf("Overlap integrals : \n");
    for (int ias = 0; ias < lapw_global.atoms.size(); ias++)
    {
        printf("  atom : %i\n", ias);
        
        for (int l = 0; l <= lapw_global.lmaxapw; l++)
        {
            printf("     l : %i\n", l);
            int ordmax = lapw_global.atoms[ias]->species->rfmt_order[l];
            for (int io1 = 0; io1 < ordmax; io1++)
            {
                for (int io2 = 0; io2 < ordmax; io2++)
                {
                    printf("%12.6f", lapw_runtime.ovlprad(l, io1, io2, ias));
                }
                printf("\n");
            }
        }
    }*/
}
