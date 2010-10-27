/*
	decl.h
		Type declarations
		this file is part of Vegas
		last modified 13 Sep 10 th
*/


#include "stddecl.h"

#define MAXGRIDS 10

#define NBINS 128

typedef unsigned char bin_t;
/* Note: bin_t must be wide enough to hold the numbers 0..NBINS */

typedef const bin_t cbin_t;

typedef real Grid[NBINS];

typedef struct {
  real sum, sqsum;
  real weightsum, avgsum;
  real chisum, chisqsum, guess;
  real avg, err, chisq;
} Cumulants;

typedef const Cumulants cCumulants;

typedef int (*Integrand)(ccount *, creal *, ccount *, real *,
  void *, creal *, cint *);

typedef struct _this {
  count ndim, ncomp;
#ifndef MLVERSION
  Integrand integrand;
  void *userdata;
#endif
  real epsrel, epsabs;
  int flags, seed;
  number mineval, maxeval;
  number nstart, nincrease, nbatch;
  int gridno;
  cchar *statefile;
  number neval;
  RNGState rng;
  jmp_buf abort;
} This;

typedef const This cThis;

static Grid *gridptr_[MAXGRIDS];
static count griddim_[MAXGRIDS];

