#include "anyfft.h"
#include "fftn.h"

void fftnd(Complex *data, int ndim, int *dims, int isign ) {
  fftn(ndim,dims,(Real *)data,((Real *)data)+1,isign*2,(isign==1 ? 1.0 : -1.0));
  fft_free();
}
