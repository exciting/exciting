#include "lapw.h"

extern "C" void FORTRAN(lapw_get_dmatu)(complex16* dmatu_)
{
    memcpy(dmatu_, lapw_runtime.dmatu.get_ptr(), lapw_runtime.dmatu.size() * sizeof(complex16)); 
}
 
