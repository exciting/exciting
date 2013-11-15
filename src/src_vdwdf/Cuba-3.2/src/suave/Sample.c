/*
	Sample.c
		the sampling step of Suave
		this file is part of Suave
		last modified 30 Jul 13 th
*/


typedef struct {
  real sum, sqsum;
  real weight, weightsum, avg, avgsum;
  real guess, chisum, chisqsum;
} Cumulants;

/*********************************************************************/

static void Sample(This *t, cnumber nnew, Region *region,
  real *lastw, real *lastx, real *lastf)
{
  count comp, df;
  number n;
  Vector(Cumulants, cumul, NCOMP);
  Cumulants *C = cumul + t->ncomp, *c;
  Bounds *bounds = RegionBounds(region), *B = bounds + t->ndim, *b;
  Result *res;
  char **ss = NULL, *s = NULL;
  ccount chars = 128*(region->div + 1);

  creal jacobian = 1/ldexp((real)nnew, region->div);
  real *w = lastw, *f = lastx;
  bin_t *bin = (bin_t *)(lastf + nnew*t->ncomp);

  for( n = nnew; n; --n ) {
    real weight = jacobian;

    t->rng.getrandom(t, f);

    for( b = bounds; b < B; ++b ) {
      creal pos = *f*NBINS;
      ccount ipos = (count)pos;
      creal prev = (ipos == 0) ? 0 : b->grid[ipos - 1];
      creal diff = b->grid[ipos] - prev;
      *f++ = b->lower + (prev + (pos - ipos)*diff)*(b->upper - b->lower);
      *bin++ = ipos;
      weight *= diff*NBINS;
    }

    *w++ = weight;
  }

  DoSample(t, nnew, lastx, lastf, lastw, region->div + 1);

  w[-1] = -w[-1];
  lastw = w;
  w = RegionW(region);
  region->n = lastw - w;

  if( VERBOSE > 2 ) {
    char *p0;
    MemAlloc(ss, t->ndim*64 + t->ncomp*(sizeof(char *) + chars));
    s = (char *)(ss + t->ncomp);
    p0 = s + t->ndim*64;
    for( comp = 0; comp < t->ncomp; ++comp ) {
      ss[comp] = p0;
      p0 += chars;
    }
  }

  FClear(cumul);
  df = n = 0;

  while( w < lastw ) {
    cbool final = (*w < 0);
    creal weight = fabs(*w++);
    ++n;

    for( c = cumul, comp = 0; c < C; ++c ) {
      creal wfun = weight*(*f++);
      c->sum += wfun;
      c->sqsum += Sq(wfun);

      if( final ) {
        if( n > 1 ) {
          real w = Weight(c->sum, c->sqsum, n);
          c->weightsum += c->weight = w;
          c->avgsum += c->avg = w*c->sum;

          if( VERBOSE > 2 ) {
            creal sig = sqrt(1/w);
            ss[comp] += (df == 0) ?
              sprintf(ss[comp], "\n[" COUNT "] "
                REAL " +- " REAL " (" NUMBER ")", comp + 1,
                c->sum, sig, n) :
              sprintf(ss[comp], "\n    "
                REAL " +- " REAL " (" NUMBER ")",
                c->sum, sig, n);
          }

          if( df == 0 ) c->guess = c->sum;
          else {
            c->chisum += w *= c->sum - c->guess;
            c->chisqsum += w*c->sum;
          }
        }

        c->sum = c->sqsum = 0;
      }
    }

    if( final ) ++df, n = 0;
  }

  region->df = --df;

  for( c = cumul, res = region->result; c < C; ++c, ++res ) {
    creal sigsq = 1/c->weightsum;
    creal avg = sigsq*c->avgsum;

    if( LAST ) {
      res->sigsq = 1/c->weight;
      res->avg = res->sigsq*c->avg;
    }
    else {
      res->sigsq = sigsq;
      res->avg = avg;
    }
    res->err = sqrt(res->sigsq);

    res->chisq = (sigsq < .9*NOTZERO) ? 0 : c->chisqsum - avg*c->chisum;
      /* This catches the special case where the integrand is constant
         over the entire region.  Unless that constant is zero, only the
         first set of samples will have zero variance, and hence weight
         (n - 1) 1e30 (see above).  All other sets have been sampled
         from a non-constant weight function and therefore inevitably
         show some variance.  This is an artificial effect, brought about
         by the fact that the constancy of the integrand in the region is
         seen only in this subdivision, and can degrade the chi-square
         score quite a bit.  If the constancy was determined from more
         than two samples (hence .9*NOTZERO), the chi-squares from the
         other sets are removed here. */
  }

  if( VERBOSE > 2 ) {
    char *p = s;
    char *p0 = p + t->ndim*64;
    char *msg = "\nRegion (" REALF ") - (" REALF ")";

    for( b = bounds; b < B; ++b ) {
      p += sprintf(p, msg, b->lower, b->upper);
      msg = "\n       (" REALF ") - (" REALF ")";
    }

    for( comp = 0, res = region->result;
         comp < t->ncomp; ++comp, ++res ) {
      p += sprintf(p, "%s  \tchisq " REAL " (" COUNT " df)",
        p0, res->chisq, df);
      p0 += chars;
    }

    Print(s);
    free(ss);
  }
}

