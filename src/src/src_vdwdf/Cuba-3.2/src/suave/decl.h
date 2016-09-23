/*
	decl.h
		Type declarations
		this file is part of Suave
		last modified 29 Jul 13 th
*/


#include "stddecl.h"

#define MINSAMPLES 10

#define NBINS 64

typedef unsigned char bin_t;
/* Note: bin_t must be wide enough to hold the numbers 0..NBINS */

typedef const bin_t cbin_t;

typedef real Grid[NBINS];

typedef const Grid cGrid;

typedef struct {
  real avg, err, sigsq, chisq;
} Result;

typedef const Result cResult;

typedef struct {
  real lower, upper;
  Grid grid;
} Bounds;

typedef const Bounds cBounds;

typedef int (*Integrand)(ccount *, creal *, ccount *, real *,
  void *, creal *, cint *);

typedef struct _this {
  count ndim, ncomp;
#ifndef MLVERSION
  Integrand integrand;
  void *userdata;
#ifdef HAVE_FORK
  int ncores, *child;
  real *frame;
  SHM_ONLY(int shmid;)
#endif
#endif
  real epsrel, epsabs;
  int flags, seed;
  number mineval, maxeval;
  number nnew;
  real flatness;
  cchar *statefile;
  count nregions;
  number neval;
  RNGState rng;  
  jmp_buf abort;
} This;

#define nframe nnew

typedef const This cThis;

typedef struct region {
  struct region *next;
  size_t size;
  count div, df;
  number n;
  Result result[];
} Region;

#define RegionBounds(r) ((Bounds *)(r->result + t->ncomp))
#define RegionW(r) ((real *)(RegionBounds(r) + t->ndim))

