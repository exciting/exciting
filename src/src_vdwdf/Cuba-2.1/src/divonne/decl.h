/*
	decl.h
		Type declarations
		this file is part of Divonne
		last modified 8 Jun 10 th
*/


#include "stddecl.h"

#define Tag(x) ((x) | INT_MIN)
#define Untag(x) ((x) & INT_MAX)

typedef struct {
  real lower, upper;
} Bounds;

typedef const Bounds cBounds;

typedef struct {
  real avg, spreadsq;
  real spread, secondspread;
  real nneed, maxerrsq, mindevsq;
  int iregion;
} Totals;

typedef struct {
  void *first, *last;
  real errcoeff[3];
  count n;
} Rule;

typedef const Rule cRule;

typedef struct samples {
  real weight;
  real *x, *f, *avg, *err;
  void (*sampler)(struct _this *t, const struct samples *, cBounds *, creal);
  cRule *rule;
  count coeff;
  number n, neff;
} Samples;

typedef const Samples cSamples;

typedef int (*Integrand)(ccount *, creal *, ccount *, real *, void *, cint *);

typedef void (*PeakFinder)(ccount *, cBounds *, number *, real *);

typedef struct _this {
  count ndim, ncomp;
#ifndef MLVERSION
  Integrand integrand;
  void *userdata;
  PeakFinder peakfinder;
#endif
  real epsrel, epsabs;
  int flags, seed;
  number mineval, maxeval;
  int key1, key2, key3;
  count maxpass;
  Bounds border;
  real maxchisq, mindeviation;
  number ngiven, nextra;
  real *xgiven, *xextra, *fgiven, *fextra;
  count ldxgiven;
  count nregions;
  number neval, neval_opt, neval_cut;
  int phase;
  count selectedcomp, size;
  Samples samples[3];
  Totals *totals;
  Rule rule7, rule9, rule11, rule13;
  RNGState rng;
  void *voidregion;
  jmp_buf abort;
} This;

typedef const This cThis;

#define CHUNKSIZE 4096

#define TYPEDEFREGION \
  typedef struct { \
    real avg, err, spread, chisq; \
    real fmin, fmax; \
    real xmin[NDIM], xmax[NDIM]; \
  } Result; \
  typedef const Result cResult; \
  typedef struct region { \
    count cutcomp, depth, xmajor; \
    real fmajor, fminor, vol; \
    Bounds bounds[NDIM]; \
    Result result[NCOMP]; \
  } Region

#define RegionPtr(n) (&((Region *)t->voidregion)[n])

