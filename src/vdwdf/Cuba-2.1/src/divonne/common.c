/*
	common.c
		includes most of the modules
		this file is part of Divonne
		last modified 8 Jun 10 th
*/


#include "Random.c"
#include "ChiSquare.c"
#include "Rule.c"
#include "Sample.c"
#include "FindMinimum.c"

static bool Explore(This *t, count iregion, cSamples *samples,
  cint depth, cint flags);

static void Split(This *t, count iregion, int depth);

#include "Explore.c"
#include "Split.c"

static inline bool BadDimension(cThis *t, ccount key)
{
  if( t->ndim > NDIM ) return true;
  if( IsSobol(key) ) return
    t->ndim < SOBOL_MINDIM || (t->seed == 0 && t->ndim > SOBOL_MAXDIM);
  if( IsRule(key, t->ndim) ) return t->ndim < 1;
  return t->ndim < KOROBOV_MINDIM || t->ndim > KOROBOV_MAXDIM;
}

static inline bool BadComponent(cThis *t)
{
  if( t->ncomp > NCOMP ) return true;
  return t->ncomp < 1;
}

#include "Integrate.c"

