/*
	common.c
		includes most of the modules
		this file is part of Suave
		last modified 2 Jun 10 th
*/


#define RegionAlloc(p, n, nnew) MemAlloc(p, \
  sizeof(Region) + \
  (n)*(t->ndim + t->ncomp + 1)*sizeof(real) + \
  (nnew)*t->ndim*sizeof(bin_t))

static inline bool BadDimension(cThis *t)
{
  if( t->ndim > NDIM ) return true;
  return t->ndim < SOBOL_MINDIM ||
    (t->seed == 0 && t->ndim > SOBOL_MAXDIM);
}

static inline bool BadComponent(cThis *t)
{
  if( t->ncomp > NCOMP ) return true;
  return t->ncomp < 1;
}

#include "Random.c"
#include "ChiSquare.c"
#include "Grid.c"
#include "Sample.c"
#include "Fluct.c"
#include "Integrate.c"

