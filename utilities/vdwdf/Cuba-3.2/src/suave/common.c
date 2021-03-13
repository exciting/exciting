/*
	common.c
		includes most of the modules
		this file is part of Suave
		last modified 29 Jul 13 th
*/


static inline Region *RegionAlloc(cThis *t, cnumber n, cnumber nnew)
{
  csize_t size = sizeof(Region) +
    t->ncomp*sizeof(Result) +
    t->ndim*sizeof(Bounds) +
    t->ncomp*t->ndim*2*sizeof(real) +
    n*SAMPLESIZE +
    nnew*t->ndim*sizeof(bin_t);
  Region *p;
  MemAlloc(p, size);
  p->size = size;
  return p;
}

static inline bool BadDimension(cThis *t)
{
  if( t->ndim > MAXDIM ) return true;
  return t->ndim < SOBOL_MINDIM ||
    (t->seed == 0 && t->ndim > SOBOL_MAXDIM);
}

static inline bool BadComponent(cThis *t)
{
  if( t->ncomp > MAXCOMP ) return true;
  return t->ncomp < 1;
}

#include "Random.c"
#include "ChiSquare.c"
#include "Grid.c"
#include "Sample.c"
#include "Fluct.c"

