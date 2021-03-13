/*
	decl.h
		Type declarations
		this file is part of Cuhre
		last modified 26 Jul 13 th
*/


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

enum { nrules = 5 };

typedef struct {
  count n;
  real weight[nrules], scale[nrules], norm[nrules];
  real gen[];
} Set;

#define SetSize (sizeof(Set) + t->ndim*sizeof(real))

typedef struct {
  Set *first, *last;
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
#ifdef HAVE_FORK
  int ncores, *child;
  SHM_ONLY(int shmid;)
#endif
#endif
  real *frame;
  real epsrel, epsabs;
  int flags;
  number mineval, maxeval;
  count key, nregions;
  cchar *statefile;
  number neval;
  Rule rule;
  jmp_buf abort;
} This;

#define nframe rule.n

typedef const This cThis;

typedef struct region {
  count div;
  Bounds bounds[];
} Region;

#define RegionSize (sizeof(Region) + t->ndim*sizeof(Bounds) + t->ncomp*sizeof(Result))

#define RegionResult(r) ((Result *)(r->bounds + t->ndim))

#define RegionPtr(p, n) ((Region *)((char *)p->region + (n)*regionsize))

