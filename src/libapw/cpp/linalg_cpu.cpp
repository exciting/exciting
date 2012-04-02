#include "linalg_cpu.h"

void zcopy(int32_t n, complex16 *zx, int32_t incx, complex16 *zy, int32_t incy)
{
    FORTRAN(zcopy)(&n, zx, &incx, zy, &incy);
}

