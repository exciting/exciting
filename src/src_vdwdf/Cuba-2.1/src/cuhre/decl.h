/*
	decl.h
		Type declarations
		this file is part of Cuhre
		last modified 7 Jun 10 th */


#include "stddecl.h"

typedef struct {
  real avg, err;
  count bisectdim;
} Result;

typedef const Result cResult;

typedef struct {
  real avg, err, lastavg, lasterr;
  real weightsum, avgsum;
  real guess, chisum, chisqsum, chisq;
} Totals;

typedef const Totals cTotals;

typedef struct {
  real lower, upper;
} Bounds;

typedef const Bounds cBounds;

typedef struct {
  real *x, *f;
  void *first, *last;
  real errcoeff[3];
  count n;
} Rule;

typedef const Rule cRule;

typedef int (*Integrand)(ccount *, creal *, ccount *, real *, void *);

typedef struct _this {
  count ndim, ncomp;
#ifndef MLVERSION
  Integrand integrand;
  void *userdata;
#endif
  real epsrel, epsabs;
  int flags;
  number mineval, maxeval;
  count key, nregions;
  number neval;
  Rule rule;
  jmp_buf abort;
} This;

typedef const This cThis;

#define TYPEDEFREGION \
  typedef struct region { \
    count div; \
    Result result[NCOMP]; \
    Bounds bounds[NDIM]; \
  } Region

