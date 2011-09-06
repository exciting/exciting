/*
	Divonne.c
		Multidimensional integration by partitioning
		originally by J.H. Friedman and M.H. Wright
		(CERNLIB subroutine D151)
		this version by Thomas Hahn
		last modified 16 Jun 10 th
*/

#include "decl.h"

#define Print(s) puts(s); fflush(stdout)

/*********************************************************************/

static inline void DoSample(This *t, number n, ccount ldx, creal *x, real *f)
{
  t->neval += n;
  while( n-- ) {
    if( t->integrand(&t->ndim, x, &t->ncomp, f, t->userdata, &t->phase) == ABORT )
      longjmp(t->abort, 1);
    x += ldx;
    f += t->ncomp;
  }
}

/*********************************************************************/

static inline count SampleExtra(This *t, cBounds *b)
{
  number n = t->nextra;
  t->peakfinder(&t->ndim, b, &n, t->xextra);
  DoSample(t, n, t->ldxgiven, t->xextra, t->fextra);
  return n;
}

/*********************************************************************/

static inline void AllocGiven(This *t, creal *xgiven)
{
  if( t->ngiven | t->nextra ) {
    cnumber nxgiven = t->ngiven*(t->ldxgiven = IMax(t->ldxgiven, t->ndim));
    cnumber nxextra = t->nextra*t->ldxgiven;
    cnumber nfgiven = t->ngiven*t->ncomp;
    cnumber nfextra = t->nextra*t->ncomp;

    Alloc(t->xgiven, nxgiven + nxextra + nfgiven + nfextra);
    t->xextra = t->xgiven + nxgiven;
    t->fgiven = t->xextra + nxextra;
    t->fextra = t->fgiven + nfgiven;

    if( nxgiven ) {
      t->phase = 0;
      Copy(t->xgiven, xgiven, nxgiven);
      DoSample(t, t->ngiven, t->ldxgiven, t->xgiven, t->fgiven);
    }
  }
}

/*********************************************************************/

#include "common.c"

Extern void EXPORT(Divonne)(ccount ndim, ccount ncomp,
  Integrand integrand, void *userdata,
  creal epsrel, creal epsabs,
  cint flags, cint seed,
  cnumber mineval, cnumber maxeval,
  cint key1, cint key2, cint key3, ccount maxpass,
  creal border, creal maxchisq, creal mindeviation,
  cnumber ngiven, ccount ldxgiven, creal *xgiven,
  cnumber nextra, PeakFinder peakfinder,
  int *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  This t;
  t.ndim = ndim;
  t.ncomp = ncomp;
  t.integrand = integrand;
  t.userdata = userdata;
  t.epsrel = epsrel;
  t.epsabs = epsabs;
  t.flags = flags;
  t.seed = seed;
  t.mineval = mineval;
  t.maxeval = maxeval;
  t.key1 = key1;
  t.key2 = key2;
  t.key3 = key3;
  t.maxpass = maxpass;
  t.border.upper = 1 - (t.border.lower = border);
  t.maxchisq = maxchisq;
  t.mindeviation = mindeviation;
  t.ngiven = ngiven;
  t.xgiven = NULL;
  t.ldxgiven = ldxgiven;
  t.nextra = nextra;
  t.peakfinder = peakfinder;
  t.nregions = 0;
  t.neval = 0;

  AllocGiven(&t, xgiven);

  *pfail = Integrate(&t, integral, error, prob);
  *pnregions = t.nregions;
  *pneval = t.neval;

  free(t.xgiven);
}

/*********************************************************************/

Extern void EXPORT(divonne)(ccount *pndim, ccount *pncomp,
  Integrand integrand, void *userdata,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cint *pseed,
  cnumber *pmineval, cnumber *pmaxeval,
  cint *pkey1, cint *pkey2, cint *pkey3, ccount *pmaxpass,
  creal *pborder, creal *pmaxchisq, creal *pmindeviation,
  cnumber *pngiven, ccount *pldxgiven, creal *xgiven,
  cnumber *pnextra, PeakFinder peakfinder,
  int *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  This t;
  t.ndim = *pndim;
  t.ncomp = *pncomp;
  t.integrand = integrand;
  t.userdata = userdata;
  t.epsrel = *pepsrel;
  t.epsabs = *pepsabs;
  t.flags = *pflags;
  t.seed = *pseed;
  t.mineval = *pmineval;
  t.maxeval = *pmaxeval;
  t.key1 = *pkey1;
  t.key2 = *pkey2;
  t.key3 = *pkey3;
  t.maxpass = *pmaxpass;
  t.border.upper = 1 - (t.border.lower = *pborder);
  t.maxchisq = *pmaxchisq;
  t.mindeviation = *pmindeviation;
  t.ngiven = *pngiven;
  t.xgiven = NULL;
  t.ldxgiven = *pldxgiven;
  t.nextra = *pnextra;
  t.peakfinder = peakfinder;
  t.nregions = 0;
  t.neval = 0;

  AllocGiven(&t, xgiven);

  *pfail = Integrate(&t, integral, error, prob);
  *pnregions = t.nregions;
  *pneval = t.neval;

  free(t.xgiven);
}

