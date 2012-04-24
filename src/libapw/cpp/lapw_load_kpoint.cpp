#include "lapw.h"

extern "C" void FORTRAN(lapw_load_kpoint)(int32_t *ngp_, int32_t *igpig_, double *vgpc_, double *weight_)
{
    int ngp = *ngp_;
    
    bloch_states_k *ks = new bloch_states_k(ngp);
    
    ks->vgkc.set_ptr(vgpc_);
    ks->vgkc.set_dimensions(3, lapw_global.ngkmax);
    
    for (int ig = 0; ig < ngp; ig++)
    {
        ks->idxg[ig] = igpig_[ig] - 1;
        ks->idxgfft[ig] = lapw_global.igfft[igpig_[ig] - 1];
    }
    ks->weight = *weight_;
    
    lapw_runtime.bloch_states.push_back(ks);
}
