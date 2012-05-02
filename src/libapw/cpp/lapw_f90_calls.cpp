#include "lapw.h"

extern "C" void FORTRAN(zfftifc)(int32_t *dim, int32_t *ngrid, int32_t *dir, complex16 *data);

void lapw_fft(int32_t direction, complex16 *data)
{
    int dim = 3;
    FORTRAN(zfftifc)(&dim, &lapw_global.ngrid[0], &direction, data);
}

extern "C" void FORTRAN(fderiv)(int32_t *m, int32_t *n, double *x, double *f, double *g, double *cf);

double lapw_spline_integrate(int32_t n, double *x, double *f)
{
    std::vector<double> g(n);
    std::vector<double> cf(4*n);
    
    int32_t m = -1;
    FORTRAN(fderiv)(&m, &n, x, f, &g[0], &cf[0]);
    
    return g[n-1];
}


