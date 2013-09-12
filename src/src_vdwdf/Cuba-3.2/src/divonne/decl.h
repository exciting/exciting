/*
	decl.h
		Type declarations
		this file is part of Divonne
		last modified 26 Jul 13 th
*/


#include "stddecl.h"

#define INIDEPTH 3
#define DEPTH 5
#define POSTDEPTH 15

#define Tag(x) ((x) | INT_MIN)
#define Untag(x) ((x) & INT_MAX)

typedef struct {
  real lower, upper;
} Bounds;

typedef const Bounds cBounds;

typedef struct {
  real avg, err;
} PhaseResult;

typedef struct {
  real avg, spreadsq;
  real spread, secondspread;
  real nneed, maxerrsq, mindevsq;
  real integral, sigsq, chisq;
  PhaseResult phase[2];
  int iregion;
} Totals;

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

typedef struct samples {
  real *x, *f;
  void (*sampler)(struct _this *t, ccount);
  cRule *rule;
  number n, neff;
  count coeff;
} Samples;

typedef const Samples cSamples;

typedef struct {
  real diff, err, spread;
} Errors;

typedef const Errors cErrors;

typedef struct {
  real avg, err, spread, chisq;
  real fmin, fmax;
  real xminmax[];
} Result;

typedef const Result cResult;

#define ResultSize (sizeof(Result) + t->ndim*2*sizeof(real))

typedef struct region {
  int depth, next;
  count isamples, cutcomp, xmajor;
  real fmajor, fminor, vol;
  Bounds bounds[];
} Region;

#define RegionSize (sizeof(Region) + t->ndim*sizeof(Bounds) + t->ncomp*ResultSize)

#define RegionResult(r) ((Result *)(r->bounds + t->ndim))

#define RegionPtr(n) ((Region *)((char *)t->region + (n)*regionsize))


typedef int (*Integrand)(ccount *, creal *, ccount *, real *, void *, cint *);

typedef void (*PeakFinder)(ccount *, cBounds *, number *, real *);

typedef struct _this {
  count ndim, ncomp;
#ifndef MLVERSION
  Integrand integrand;
  void *userdata;
  PeakFinder peakfinder;
#ifdef HAVE_FORK
  int ncores, running, *child;
  real *frame;
  number nframe;
  SHM_ONLY(int shmid;)
#endif
#endif
  real epsrel, epsabs;
  int flags, seed;
  number mineval, maxeval;
  int key1, key2, key3;
  count maxpass;
  Bounds border;
  real maxchisq, mindeviation;
  number ngiven, nextra;
  real *xgiven, *fgiven;
  real *xextra, *fextra;
  count ldxgiven;
  count nregions;
  cchar *statefile;
  number neval, neval_opt, neval_cut, nrand;
  count phase;
  count selectedcomp, size;
  Samples samples[3];
  Totals *totals;
  Rule rule7, rule9, rule11, rule13;
  RNGState rng;
  Region *region;
  jmp_buf abort;
} This;

typedef const This cThis;


#define CHUNKSIZE 4096

#define AllocRegions(t) \
  MemAlloc((t)->region, (t)->size*regionsize)

#define EnlargeRegions(t, n) if( (t)->nregions + n > (t)->size ) \
  ReAlloc((t)->region, ((t)->size += CHUNKSIZE)*regionsize)

#define SAMPLERDEFS \
  csize_t regionsize = RegionSize; \
  Region *region = RegionPtr(iregion); \
  cBounds *b = region->bounds; \
  Result *res = RegionResult(region); \
  cSamples *samples = &t->samples[region->isamples]; \
  real *x = samples->x, *f = samples->f; \
  cnumber n = samples->n

