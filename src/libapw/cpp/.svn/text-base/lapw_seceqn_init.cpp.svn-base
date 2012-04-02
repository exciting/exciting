#include "lapw.h"

extern "C" void FORTRAN(lapw_seceqn_init)(double *hmltrad_, double *ovlprad_, double *beffrad_,
                                          double *apwfr_, double *apwdfr_, double *beffir_,
                                          complex16 *veffig_)
                                     
{
    p.hmltrad.set_dimensions(p.lmmaxvr, p.nrfmtmax, p.nrfmtmax, p.natmtot);
    p.hmltrad.set_ptr(hmltrad_);
    
    p.ovlprad.set_dimensions(p.lmaxapw + 1, p.ordrfmtmax, p.ordrfmtmax, p.natmtot);
    p.ovlprad.set_ptr(ovlprad_);
    
    p.apwfr.set_dimensions(p.nrmtmax, 2, p.apwordmax, p.lmaxapw + 1, p.natmtot);
    p.apwfr.set_ptr(apwfr_);

    p.apwdfr.set_dimensions(p.apwordmax, p.lmaxapw + 1, p.natmtot);
    p.apwdfr.set_ptr(apwdfr_);

    p.veffig.set_dimensions(p.ngvec);
    p.veffig.set_ptr(veffig_);

    if (p.ndmag > 0) 
    {
        p.beffrad.set_dimensions(p.lmmaxvr, p.nrfmtmax, p.nrfmtmax, p.natmtot, p.ndmag);
        p.beffrad.set_ptr(beffrad_);

        p.beffir.set_dimensions(p.ngrtot, p.ndmag);
        p.beffir.set_ptr(beffir_);
        
        if (lapw_diag == full)
        {
            p.beffig.set_dimensions(p.ngvec, p.ndmag);
            p.beffig.allocate();
            std::vector<complex16> zfft(p.ngrtot);
            for (unsigned int i = 0; i < p.ndmag; i++)
            {
                for (unsigned int ir = 0; ir < p.ngrtot; ir++) zfft[ir] = zone * p.beffir(ir, i) * p.cfunir[ir];
                lapw_fft(-1, &zfft[0]);
                for (unsigned int ig = 0; ig < p.ngvec; ig++) p.beffig(ig, i) = zfft[p.igfft[ig]];
            }
        }
    }
}
