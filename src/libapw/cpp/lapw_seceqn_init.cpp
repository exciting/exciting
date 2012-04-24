#include "lapw.h"

extern "C" void FORTRAN(lapw_seceqn_init)(double *hmltrad_, double *ovlprad_, double *beffrad_,
                                          double *apwfr_, double *apwdfr_, double *beffir_,
                                          complex16 *veffig_)
                                     
{
    lapw_runtime.hmltrad.set_dimensions(lapw_global.lmmaxvr, lapw_global.nrfmtmax, lapw_global.nrfmtmax, lapw_global.natmtot);
    lapw_runtime.hmltrad.set_ptr(hmltrad_);
    
    lapw_runtime.ovlprad.set_dimensions(lapw_global.lmaxapw + 1, lapw_global.ordrfmtmax, lapw_global.ordrfmtmax, lapw_global.natmtot);
    lapw_runtime.ovlprad.set_ptr(ovlprad_);
    
    lapw_runtime.apwfr.set_dimensions(lapw_global.nrmtmax, 2, lapw_global.apwordmax, lapw_global.lmaxapw + 1, lapw_global.natmtot);
    lapw_runtime.apwfr.set_ptr(apwfr_);

    lapw_runtime.apwdfr.set_dimensions(lapw_global.apwordmax, lapw_global.lmaxapw + 1, lapw_global.natmtot);
    lapw_runtime.apwdfr.set_ptr(apwdfr_);

    lapw_runtime.veffig.set_dimensions(lapw_global.ngvec);
    lapw_runtime.veffig.set_ptr(veffig_);

    if (lapw_global.ndmag > 0) 
    {
        lapw_runtime.beffrad.set_dimensions(lapw_global.lmmaxvr, lapw_global.nrfmtmax, lapw_global.nrfmtmax, lapw_global.natmtot, lapw_global.ndmag);
        lapw_runtime.beffrad.set_ptr(beffrad_);

        lapw_runtime.beffir.set_dimensions(lapw_global.ngrtot, lapw_global.ndmag);
        lapw_runtime.beffir.set_ptr(beffir_);
        
        if (lapw_diag == full)
        {
            lapw_runtime.beffig.set_dimensions(lapw_global.ngvec, lapw_global.ndmag);
            lapw_runtime.beffig.allocate();
            std::vector<complex16> zfft(lapw_global.ngrtot);
            for (unsigned int i = 0; i < lapw_global.ndmag; i++)
            {
                for (unsigned int ir = 0; ir < lapw_global.ngrtot; ir++) zfft[ir] = zone * lapw_runtime.beffir(ir, i) * lapw_global.cfunir[ir];
                lapw_fft(-1, &zfft[0]);
                for (unsigned int ig = 0; ig < lapw_global.ngvec; ig++) lapw_runtime.beffig(ig, i) = zfft[lapw_global.igfft[ig]];
            }
        }
    }
}
