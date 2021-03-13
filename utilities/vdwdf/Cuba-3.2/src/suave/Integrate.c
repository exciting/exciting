/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Suave
		checkpointing by B. Chokoufe
		last modified 5 Aug 13 th
*/


typedef struct {
  signature_t signature;
  count nregions, df;
  number neval;
  Result totals[];
} State;

static int Integrate(This *t, real *integral, real *error, real *prob)
{
  StateDecl;
  csize_t statesize = sizeof(State) + NCOMP*sizeof(Result);
  Sized(State, state, statesize);
  Array(Var, var, NDIM, 2);
  Vector(char, out, 128*NCOMP + 256);

  Region *anchor = NULL, *region = NULL;
  Result *tot, *Tot = state->totals + t->ncomp;
  Result *res, *resL, *resR;
  Bounds *b, *B;
  count dim, comp;
  int fail;

  if( VERBOSE > 1 ) {
    sprintf(out, "Suave input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  epsrel " REAL "\n  epsabs " REAL "\n"
      "  flags %d\n  seed %d\n"
      "  mineval " NUMBER "\n  maxeval " NUMBER "\n"
      "  nnew " NUMBER "\n  flatness " REAL "\n"
      "  statefile \"%s\"\n",
      t->ndim, t->ncomp,
      t->epsrel, t->epsabs,
      t->flags, t->seed,
      t->mineval, t->maxeval,
      t->nnew, t->flatness,
      t->statefile);
    Print(out);
  }

  if( BadComponent(t) ) return -2;
  if( BadDimension(t) ) return -1;

  ShmAlloc(t, ShmRm(t));
  ForkCores(t);

  if( (fail = setjmp(t->abort)) ) goto abort;

  t->epsabs = Max(t->epsabs, NOTZERO);
  IniRandom(t);

  StateSetup(t);

  if( StateReadTest(t) ) {
    size_t *regsize = NULL;
    StateReadOpen(t, fd) {
      count iregion;
      size_t totsize;
      Region **last = &anchor;
      if( read(fd, state, statesize) != statesize ||
          state->signature != StateSignature(t, 2) ) break;
      t->nregions = state->nregions;
      totsize = t->nregions*sizeof *regsize;
      MemAlloc(regsize, totsize);
      StateRead(fd, regsize, totsize);
      for( iregion = 0; iregion < t->nregions; ++iregion )
        totsize += regsize[iregion];
      if( st.st_size != statesize + totsize ) break;
      for( iregion = 0; iregion < t->nregions; ++iregion ) {
        Region *reg;
        MemAlloc(reg, regsize[iregion]);
        StateRead(fd, reg, regsize[iregion]);
        *last = reg;
        last = &reg->next;
      }
      *last = NULL;
    } StateReadClose(t, fd);
    free(regsize);
    t->neval = state->neval;
    t->rng.skiprandom(t, t->neval);
  }

  if( ini ) {
    count bin;
    Bounds *b0;

    anchor = RegionAlloc(t, t->nnew, t->nnew);
    anchor->next = NULL;
    anchor->div = 0;
    t->nregions = 1;

    b0 = RegionBounds(anchor);
    b0->lower = 0;
    b0->upper = 1;
    /* define the initial distribution of bins */
    for( bin = 0; bin < NBINS; ++bin )
      b0->grid[bin] = (bin + 1)/(real)NBINS;

    for( b = b0 + 1, B = b0 + t->ndim; b < B; ++b )
      Copy(b, b0, 1);

    t->neval = 0;
    Sample(t, t->nnew, anchor, RegionW(anchor),
      RegionW(anchor) + t->nnew,
      RegionW(anchor) + t->nnew + t->ndim*t->nnew);
    state->df = anchor->df;
    FCopy(state->totals, anchor->result);
  }

  /* main iteration loop */
  for( ; ; ) {
    Var *vLR;
    real maxratio, maxerr, minfluct, bias, mid;
    Region *regionL, *regionR, *reg, **parent, **par;
    Bounds *bounds, *boundsL, *boundsR;
    count maxcomp, bisectdim;
    number n, nL, nR, nnewL, nnewR;
    real *w, *wL, *wR, *x, *xL, *xR, *f, *fL, *fR, *wlast, *flast;

    if( VERBOSE ) {
      char *oe = out + sprintf(out, "\n"
        "Iteration " COUNT ":  " NUMBER " integrand evaluations so far",
        t->nregions, t->neval);
      for( tot = state->totals, comp = 0; tot < Tot; ++tot )
        oe += sprintf(oe, "\n[" COUNT "] "
          REAL " +- " REAL "  \tchisq " REAL " (" COUNT " df)",
          ++comp, tot->avg, tot->err, tot->chisq, state->df);
      Print(out);
    }

    maxratio = -INFTY;
    maxcomp = 0;
    for( tot = state->totals, comp = 0; tot < Tot; ++tot ) {
      creal ratio = tot->err/MaxErr(tot->avg);
      if( ratio > maxratio ) {
        maxratio = ratio;
        maxcomp = comp;
      }
    }

    if( maxratio <= 1 && t->neval >= t->mineval ) break;

    if( t->neval >= t->maxeval ) {
      fail = 1;
      break;
    }

    maxerr = -INFTY;
    parent = &anchor;
    region = anchor;
    for( par = &anchor; (reg = *par); par = &reg->next ) {
      creal err = reg->result[maxcomp].err;
      if( err > maxerr ) {
        maxerr = err;
        parent = par;
        region = reg;
      }
    }

    bounds = RegionBounds(region);
    w = RegionW(region);

    /* find the bisectdim which minimizes the fluctuations */
    Fluct(t, var[0], bounds, w, region->n, maxcomp,
      region->result[maxcomp].avg, Max(maxerr, t->epsabs));

    bias = (t->epsrel < 1e-50) ? 2 :
      Max(pow(2., -(real)region->div/t->ndim)/t->epsrel, 2.);
    minfluct = INFTY;
    bisectdim = 0;
    for( dim = 0; dim < t->ndim; ++dim ) {
      creal fluct = (var[dim][0].fluct + var[dim][1].fluct)*
        (bias - bounds[dim].upper + bounds[dim].lower);
      if( fluct < minfluct ) {
        minfluct = fluct;
        bisectdim = dim;
      }
    }

    /* apply stratified sampling to distribute points in bisected region */
    vLR = var[bisectdim];
    minfluct = vLR[0].fluct + vLR[1].fluct;
    nnewL = IMax(
      (minfluct == 0) ? t->nnew/2 : (count)(vLR[0].fluct/minfluct*t->nnew),
      MINSAMPLES );
    nL = vLR[0].n + nnewL;
    nnewR = IMax(t->nnew - nnewL, MINSAMPLES);
    nR = vLR[1].n + nnewR;

    regionL = RegionAlloc(t, nL, nnewL);
    regionR = RegionAlloc(t, nR, nnewR);

    *parent = regionL;
    regionL->next = regionR;
    regionR->next = region->next;
    regionL->div = regionR->div = region->div + 1;

    mid = .5*(bounds[bisectdim].lower + bounds[bisectdim].upper);
    n = region->n;
    wlast = w;              x = w + n;     f = flast = x + n*t->ndim;
    wL = RegionW(regionL);  xL = wL + nL;  fL = xL + nL*t->ndim;
    wR = RegionW(regionR);  xR = wR + nR;  fR = xR + nR*t->ndim;
    while( n-- ) {
      cbool final = (*w < 0);
      if( x[bisectdim] < mid ) {
        if( final && wR > RegionW(regionR) ) wR[-1] = -fabs(wR[-1]);
        *wL++ = *w++;
        XCopy(xL, x);
        xL += t->ndim;
        FCopy(fL, f);
        fL += t->ncomp;
      }
      else {
        if( final && wL > RegionW(regionL) ) wL[-1] = -fabs(wL[-1]);
        *wR++ = *w++;
        XCopy(xR, x);
        xR += t->ndim;
        FCopy(fR, f);
        fR += t->ncomp;
      }
      x += t->ndim;
      f += t->ncomp;
      if( n && final ) wlast = w, flast = f;
    }

    Reweight(t, bounds, wlast, flast, f, state->totals);
    boundsL = RegionBounds(regionL);
    XCopy(boundsL, bounds);
    boundsR = RegionBounds(regionR);
    XCopy(boundsR, bounds);

    boundsL[bisectdim].upper = mid;
    boundsR[bisectdim].lower = mid;
    StretchGrid(bounds[bisectdim].grid,
      boundsL[bisectdim].grid, boundsR[bisectdim].grid);

    Sample(t, nnewL, regionL, wL, xL, fL);
    Sample(t, nnewR, regionR, wR, xR, fR);

    state->df += regionL->df + regionR->df - region->df;

    for( res = region->result,
         resL = regionL->result,
         resR = regionR->result,
         tot = state->totals;
         tot < Tot; ++res, ++resL, ++resR, ++tot ) {
      real diff, sigsq;

      tot->avg += diff = resL->avg + resR->avg - res->avg;

      diff = Sq(.25*diff);
      sigsq = resL->sigsq + resR->sigsq;
      if( sigsq > 0 ) {
        creal c = Sq(1 + sqrt(diff/sigsq));
        resL->sigsq *= c;
        resR->sigsq *= c;
      }
      resL->err = sqrt(resL->sigsq += diff);
      resR->err = sqrt(resR->sigsq += diff);

      tot->sigsq += resL->sigsq + resR->sigsq - res->sigsq;
      tot->err = sqrt(tot->sigsq);

      tot->chisq += resL->chisq + resR->chisq - res->chisq;
    }

    free(region);
    region = NULL;
    ++t->nregions;

    if( StateWriteTest(t) ) {
      StateWriteOpen(t, fd) {
        Region *reg;
        size_t totsize, *regsize, *s;
        state->signature = StateSignature(t, 2);
        state->neval = t->neval;
        state->nregions = t->nregions;
        StateWrite(fd, state, statesize);
        MemAlloc(regsize, totsize = t->nregions*sizeof(size_t));
        s = regsize;
        for( reg = anchor; reg; reg = reg->next )
          *s++ = reg->size;
        StateWrite(fd, regsize, totsize);
        free(regsize);
        for( reg = anchor; reg; reg = reg->next )
          StateWrite(fd, reg, reg->size);
      } StateWriteClose(t, fd);
    }
  }

  for( tot = state->totals, comp = 0; tot < Tot; ++tot, ++comp ) {
    integral[comp] = tot->avg;
    error[comp] = tot->err;
    prob[comp] = ChiSquare(tot->chisq, state->df);
  }

#ifdef MLVERSION
  if( REGIONS ) {
    Vector(real, bounds, NDIM*2);

    MLPutFunction(stdlink, "List", 2);

    MLPutFunction(stdlink, "List", t->nregions);
    for( region = anchor; region; region = region->next ) {
      cResult *Res;
      real *d = bounds;
      for( B = (b = RegionBounds(region)) + t->ndim; b < B; ++b ) {
        *d++ = b->lower;
        *d++ = b->upper;
      }

      MLPutFunction(stdlink, "Cuba`Suave`region", 3);

      MLPutRealList(stdlink, bounds, 2*t->ndim);

      MLPutFunction(stdlink, "List", t->ncomp);
      for( Res = (res = region->result) + t->ncomp; res < Res; ++res ) {
        real r[] = {res->avg, res->err, res->chisq};
        MLPutRealList(stdlink, r, Elements(r));
      }

      MLPutInteger(stdlink, region->df);
    }
  }
#endif

abort:
  free(region);
  while( (region = anchor) ) {
    anchor = anchor->next;
    free(region);
  }
  WaitCores(t);
  ShmFree(t);

  StateRemove(t);

  return fail;
}

