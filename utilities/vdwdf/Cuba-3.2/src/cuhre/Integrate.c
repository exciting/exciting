/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Cuhre
		checkpointing by B. Chokoufe
		last modified 5 Aug 13 th
*/


#define POOLSIZE 1024

typedef struct pool {
  struct pool *next;
  Region region[];
} Pool;

typedef struct {
  signature_t signature;
  count nregions, ncur;
  number neval;
  Totals totals[];
} State;

static int Integrate(This *t, real *integral, real *error, real *prob)
{
  StateDecl;
  csize_t statesize = sizeof(State) + NCOMP*sizeof(Totals);
  Sized(State, state, statesize);
  csize_t regionsize = RegionSize;
  csize_t poolsize = sizeof(Pool) + POOLSIZE*regionsize;
  Vector(Result, result, NCOMP);
  Vector(char, out, 128*NCOMP + 256);

  Totals *tot, *Tot = state->totals + t->ncomp;
  Result *res, *resL, *resR;
  Bounds *b, *B;
  Pool *cur = NULL, *pool;
  Region *region;
  count comp, ipool, npool;
  int fail;

  if( VERBOSE > 1 ) {
    sprintf(out, "Cuhre input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  epsrel " REAL "\n  epsabs " REAL "\n"
      "  flags %d\n  mineval " NUMBER "\n  maxeval " NUMBER "\n"
      "  key " COUNT "\n"
      "  statefile \"%s\"",
      t->ndim, t->ncomp,
      t->epsrel, t->epsabs,
      t->flags, t->mineval, t->maxeval,
      t->key,
      t->statefile);
    Print(out);
  }

  if( BadComponent(t) ) return -2;
  if( BadDimension(t) ) return -1;

  t->epsabs = Max(t->epsabs, NOTZERO);

  RuleAlloc(t);
  t->mineval = IMax(t->mineval, t->rule.n + 1);
  FrameAlloc(t, ShmRm(t));
  ForkCores(t);

  if( (fail = setjmp(t->abort)) ) goto abort;

  StateSetup(t);

  if( StateReadTest(t) ) {
    StateReadOpen(t, fd) {
      Pool *prev = NULL;
      int size;
      if( read(fd, state, statesize) != statesize ||
          state->signature != StateSignature(t, 4) ) break;
      t->neval = state->neval;
      t->nregions = state->nregions;
      do {
        MemAlloc(cur, poolsize);
        cur->next = prev;
        prev = cur;
        size = read(fd, cur, poolsize);
      } while( size == poolsize );
      if( size != state->ncur*regionsize ) break;
    } StateReadClose(t, fd);
  }

  if( ini ) {
    MemAlloc(cur, poolsize);
    cur->next = NULL;
    state->ncur = t->nregions = 1;

    region = cur->region;
    region->div = 0;
    for( B = (b = region->bounds) + t->ndim; b < B; ++b ) {
      b->lower = 0;
      b->upper = 1;
    }

    t->neval = 0;
    Sample(t, region);

    for( res = RegionResult(region), tot = state->totals;
         tot < Tot; ++res, ++tot ) {
      tot->avg = tot->lastavg = tot->guess = res->avg;
      tot->err = tot->lasterr = res->err;
      tot->weightsum = 1/Max(Sq(res->err), NOTZERO);
      tot->avgsum = tot->weightsum*res->avg;
      tot->chisq = tot->chisqsum = tot->chisum = 0;
    }
  }

  /* main iteration loop */
  for( ; ; ) {
    count maxcomp, bisectdim;
    real maxratio, maxerr;
    Region *regionL, *regionR;
    Bounds *bL, *bR;

    if( VERBOSE ) {
      char *oe = out + sprintf(out, "\n"
        "Iteration " COUNT ":  " NUMBER " integrand evaluations so far",
        t->nregions, t->neval);
      for( tot = state->totals, comp = 0; tot < Tot; ++tot )
        oe += sprintf(oe, "\n[" COUNT "] "
          REAL " +- " REAL "  \tchisq " REAL " (" COUNT " df)",
          ++comp, tot->avg, tot->err, tot->chisq, t->nregions - 1);
      Print(out);
    }

    maxratio = -INFTY;
    maxcomp = 0;
    for( tot = state->totals, comp = 0; tot < Tot; ++tot, ++comp ) {
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
    regionL = cur->region;
    npool = state->ncur;
    for( pool = cur; pool; npool = POOLSIZE, pool = pool->next )
      for( ipool = 0; ipool < npool; ++ipool ) {
        Region *region = RegionPtr(pool, ipool);
        creal err = RegionResult(region)[maxcomp].err;
        if( err > maxerr ) {
          maxerr = err;
          regionL = region;
        }
      }

    if( state->ncur == POOLSIZE ) {
      Pool *prev = cur;
      MemAlloc(cur, poolsize);
      cur->next = prev;
      state->ncur = 0;
    }
    regionR = RegionPtr(cur, state->ncur++);

    regionR->div = ++regionL->div;
    FCopy(result, RegionResult(regionL));
    XCopy(regionR->bounds, regionL->bounds);

    bisectdim = result[maxcomp].bisectdim;
    bL = &regionL->bounds[bisectdim];
    bR = &regionR->bounds[bisectdim];
    bL->upper = bR->lower = .5*(bL->upper + bL->lower);

    Sample(t, regionL);
    Sample(t, regionR);

    for( res = result,
         resL = RegionResult(regionL),
         resR = RegionResult(regionR),
         tot = state->totals;
         tot < Tot; ++res, ++resL, ++resR, ++tot ) {
      real diff, err, w, avg, sigsq;

      tot->lastavg += diff = resL->avg + resR->avg - res->avg;

      diff = fabs(.25*diff);
      err = resL->err + resR->err;
      if( err > 0 ) {
        creal c = 1 + 2*diff/err;
        resL->err *= c;
        resR->err *= c;
      }
      resL->err += diff;
      resR->err += diff;
      tot->lasterr += resL->err + resR->err - res->err;

      tot->weightsum += w = 1/Max(Sq(tot->lasterr), NOTZERO);
      sigsq = 1/tot->weightsum;
      tot->avgsum += w*tot->lastavg;
      avg = sigsq*tot->avgsum;
      tot->chisum += w *= tot->lastavg - tot->guess;
      tot->chisqsum += w*tot->lastavg;
      tot->chisq = tot->chisqsum - avg*tot->chisum;

      if( LAST ) {
        tot->avg = tot->lastavg;
        tot->err = tot->lasterr;
      }
      else {
        tot->avg = avg;
        tot->err = sqrt(sigsq);
      }
    }
    ++t->nregions;

    if( StateWriteTest(t) ) {
      StateWriteOpen(t, fd) {
        Pool *prev = cur;
        state->signature = StateSignature(t, 4);
        state->nregions = t->nregions;
        state->neval = t->neval;
        StateWrite(fd, state, statesize);
        while( (prev = prev->next) ) StateWrite(fd, prev, poolsize);
        StateWrite(fd, cur, state->ncur*regionsize);
      } StateWriteClose(t, fd);
    }
  }

  for( tot = state->totals, comp = 0; tot < Tot; ++tot, ++comp ) {
    integral[comp] = tot->avg;
    error[comp] = tot->err;
    prob[comp] = ChiSquare(tot->chisq, t->nregions - 1);
  }

#ifdef MLVERSION
  if( REGIONS ) {
    MLPutFunction(stdlink, "List", 2);
    MLPutFunction(stdlink, "List", t->nregions);

    npool = state->ncur;
    for( pool = cur; pool; npool = POOLSIZE, pool = pool->next )
      for( ipool = 0; ipool < npool; ++ipool ) {
        Region const *region = RegionPtr(pool, ipool);
        Result *Res;

        MLPutFunction(stdlink, "Cuba`Cuhre`region", 2);
        MLPutRealList(stdlink, (real *)region->bounds, 2*t->ndim);

        MLPutFunction(stdlink, "List", t->ncomp);
        for( Res = (res = RegionResult(region)) + t->ncomp;
             res < Res; ++res ) {
          real r[] = {res->avg, res->err};
          MLPutRealList(stdlink, r, Elements(r));
        }
      }
  }
#endif

abort:
  while( (pool = cur) ) {
    cur = cur->next;
    free(pool);
  }

  WaitCores(t);
  FrameFree(t);
  RuleFree(t);

  StateRemove(t);

  return fail;
}

