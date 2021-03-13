/*
	common.c
		includes most of the modules
		this file is part of Cuhre
		last modified 2 Aug 13 11 th
*/


#include "ChiSquare.c"
#include "Rule.c"

static inline bool BadDimension(cThis *t)
{
  if( t->ndim > MAXDIM ) return true;
  return t->ndim < 2;
}

static inline bool BadComponent(cThis *t)
{
  if( t->ncomp > MAXCOMP ) return true;
  return t->ncomp < 1;
}

