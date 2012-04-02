#include "lapw.h"

extern "C" void FORTRAN(lapw_load_kpoint)(int32_t *ngp_, int32_t *igpig_, double *vgpc_, double *weight_)
{
    int ngp = *ngp_;
    
    kpoint *kp = new kpoint(ngp);
    kp->vgkc.set_ptr(vgpc_);
    kp->vgkc.set_dimensions(3, p.ngkmax);

    p.kpoints.push_back(kp);
    for (int ig = 0; ig < ngp; ig++)
    {
        p.kpoints.back()->idxg[ig] = igpig_[ig] - 1;
        p.kpoints.back()->idxgfft[ig] = p.igfft[igpig_[ig] - 1];
    }
    p.kpoints.back()->weight = *weight_;
}
