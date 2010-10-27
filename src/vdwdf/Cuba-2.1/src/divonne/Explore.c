/*
	Explore.c
		sample region, determine min and max, split if necessary
		this file is part of Divonne
		last modified 8 Jun 10 th
*/


typedef struct {
  real fmin, fmax;
  creal *xmin, *xmax;
} Extrema;

/*********************************************************************/

static bool Explore(This *t, count iregion, cSamples *samples,
  cint depth, cint flags)
{
#define SPLICE (flags & 1)
#define HAVESAMPLES (flags & 2)

  TYPEDEFREGION;

  count n, dim, comp, maxcomp;
  Extrema extrema[NCOMP];
  Result *r;
  creal *x;
  real *f;
  real halfvol, maxerr;
  Region *region;
  Bounds *bounds;
  Result *result;

  /* needed as of gcc 3.3 to make gcc correctly address region #@$&! */
  sizeof(*region);

  if( SPLICE ) {
    if( t->nregions == t->size ) {
      t->size += CHUNKSIZE;
      ReAlloc(t->voidregion, t->size*sizeof(Region));
    }
    VecCopy(RegionPtr(t->nregions)->bounds, RegionPtr(iregion)->bounds);
    iregion = t->nregions++;
  }
  region = RegionPtr(iregion);
  bounds = region->bounds;
  result = region->result;

  for( comp = 0; comp < t->ncomp; ++comp ) {
    Extrema *e = &extrema[comp];
    e->fmin = INFTY;
    e->fmax = -INFTY;
    e->xmin = e->xmax = NULL;
  }

  if( !HAVESAMPLES ) {
    real vol = 1;
    for( dim = 0; dim < t->ndim; ++dim ) {
      cBounds *b = &bounds[dim];
      vol *= b->upper - b->lower;
    }
    region->vol = vol;

    for( comp = 0; comp < t->ncomp; ++comp ) {
      Result *r = &result[comp];
      r->fmin = INFTY;
      r->fmax = -INFTY;
    }

    x = t->xgiven;
    f = t->fgiven;
    n = t->ngiven;
    if( t->nextra ) n += SampleExtra(t, bounds);

    for( ; n; --n ) {
      for( dim = 0; dim < t->ndim; ++dim ) {
        cBounds *b = &bounds[dim];
        if( x[dim] < b->lower || x[dim] > b->upper ) goto skip;
      }
      for( comp = 0; comp < t->ncomp; ++comp ) {
        Extrema *e = &extrema[comp];
        creal y = f[comp];
        if( y < e->fmin ) e->fmin = y, e->xmin = x;
        if( y > e->fmax ) e->fmax = y, e->xmax = x;
      }
skip:
      x += t->ldxgiven;
      f += t->ncomp;
    }

    samples->sampler(t, samples, bounds, vol);
  }

  x = samples->x;
  f = samples->f;
  for( n = samples->n; n; --n ) {
    for( comp = 0; comp < t->ncomp; ++comp ) {
      Extrema *e = &extrema[comp];
      creal y = *f++;
      if( y < e->fmin ) e->fmin = y, e->xmin = x;
      if( y > e->fmax ) e->fmax = y, e->xmax = x;
    }
    x += t->ndim;
  }
  t->neval_opt -= t->neval;

  halfvol = .5*region->vol;
  maxerr = -INFTY;
  maxcomp = -1;

  for( comp = 0; comp < t->ncomp; ++comp ) {
    Extrema *e = &extrema[comp];
    Result *r = &result[comp];
    real xtmp[NDIM], ftmp, err;

    if( e->xmin ) {	/* not all NaNs */
      t->selectedcomp = comp;
      VecCopy(xtmp, e->xmin);
      ftmp = FindMinimum(t, bounds, xtmp, e->fmin);
      if( ftmp < r->fmin ) {
        r->fmin = ftmp;
        VecCopy(r->xmin, xtmp);
      }

      t->selectedcomp = Tag(comp);
      VecCopy(xtmp, e->xmax);
      ftmp = -FindMinimum(t, bounds, xtmp, -e->fmax);
      if( ftmp > r->fmax ) {
        r->fmax = ftmp;
        VecCopy(r->xmax, xtmp);
      }
    }

    r->avg = samples->avg[comp];
    r->err = samples->err[comp];
    r->spread = halfvol*(r->fmax - r->fmin);

    err = r->spread/Max(fabs(r->avg), NOTZERO);
    if( err > maxerr ) {
      maxerr = err;
      maxcomp = comp;
    }
  }

  t->neval_opt += t->neval;

  if( maxcomp == -1 ) { /* all NaNs */
    region->depth = 0;
    return false;
  }

  region->cutcomp = maxcomp;
  r = &region->result[maxcomp];
  if( halfvol*(r->fmin + r->fmax) > r->avg ) {
    region->fminor = r->fmin;
    region->fmajor = r->fmax;
    region->xmajor = r->xmax - (real *)region->result;
  }
  else {
    region->fminor = r->fmax;
    region->fmajor = r->fmin;
    region->xmajor = r->xmin - (real *)region->result;
  }

  region->depth = IDim(depth);

  if( !HAVESAMPLES ) {
    if( samples->weight*r->spread < r->err ||
        r->spread < t->totals[maxcomp].secondspread ) region->depth = 0;
    if( region->depth == 0 )
      for( comp = 0; comp < t->ncomp; ++comp )
        t->totals[comp].secondspread =
          Max(t->totals[comp].secondspread, result[comp].spread);
  }

  if( region->depth ) Split(t, iregion, region->depth);
  return true;
}

