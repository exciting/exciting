/*
	common.c
		Code common to Vegas.c and Vegas.tm
		this file is part of Vegas
		last modified 6 Jun 10 th
*/


#include "Random.c"
#include "ChiSquare.c"
#include "Grid.c"

#define SamplesAlloc(p) MemAlloc(p, \
  t->nbatch*((t->ndim + t->ncomp + 1)*sizeof(real) + \
             t->ndim*sizeof(bin_t)))

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

#include "Integrate.c"

