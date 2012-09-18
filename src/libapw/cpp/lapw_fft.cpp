#include "lapw.h"

extern "C" void FORTRAN(zfftifc)(int32_t *dim, int32_t *ngrid, int32_t *dir, complex16 *data);

void lapw_fft(int32_t direction, complex16 *data)
{
    int dim = 3;
    FORTRAN(zfftifc)(&dim, &lapw_global.ngrid[0], &direction, data);
}

