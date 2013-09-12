/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Vegas
		last modified 8 Aug 13 th
*/


typedef struct {
  signature_t signature;
  count niter;
  number nsamples, neval;
  Cumulants cumul[];
} State;

static int Integrate(This *t, real *integral, real *error, real *prob)
{
  bin_t *bins;
  count dim, comp;
  int fail;

  StateDecl;
  csize_t statesize = sizeof(State) +
    NCOMP*sizeof(Cumulants) + NDIM*sizeof(Grid);
  Sized(State, state, statesize);
  Cumulants *c, *C = state->cumul + t->ncomp;
  Grid *state_grid = (Grid *)C;
  Array(Grid, margsum, NCOMP, NDIM);
  Vector(char, out, 128*NCOMP + 256);

  if( VERBOSE > 1 ) {
    sprintf(out, "Vegas input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  epsrel " REAL "\n  epsabs " REAL "\n"
      "  flags %d\n  seed %d\n"
      "  mineval " NUMBER "\n  maxeval " NUMBER "\n"
      "  nstart " NUMBER "\n  nincrease " NUMBER "\n"
      "  nbatch " NUMBER "\n  gridno %d\n"
      "  statefile \"%s\"",
      t->ndim, t->ncomp,
      t->epsrel, t->epsabs,
      t->flags, t->seed,
      t->mineval, t->maxeval,
      t->nstart, t->nincrease, t->nbatch,
      t->gridno, t->statefile);
    Print(out);
  }

  if( BadComponent(t) ) return -2;
  if( BadDimension(t) ) return -1;

  FrameAlloc(t, ShmRm(t));
  ForkCores(t);
  Alloc(bins, t->nbatch*t->ndim);

  if( (fail = setjmp(t->abort)) ) goto abort;

  IniRandom(t);

  StateSetup(t);

  if( StateReadTest(t) ) {
    StateReadOpen(t, fd) {
      if( read(fd, state, statesize) != statesize ||
          state->signature != StateSignature(t, 1) ) break;
    } StateReadClose(t, fd);
    t->neval = state->neval;
    t->rng.skiprandom(t, t->neval);
  }

  if( ini ) {
    state->niter = 0;
    state->nsamples = t->nstart;
    FClear(state->cumul);
    GetGrid(t, state_grid);
    t->neval = 0;
  }

  /* main iteration loop */
  for( ; ; ) {
    number nsamples = state->nsamples;
    creal jacobian = 1./nsamples;

    FClear(margsum);

    for( ; nsamples > 0; nsamples -= t->nbatch ) {
      cnumber n = IMin(t->nbatch, nsamples);
      real *w = t->frame;
      real *x = w + n;
      real *f = x + n*t->ndim;
      real *lastf = f + n*t->ncomp;
      bin_t *bin = bins;

      while( x < f ) {
        real weight = jacobian;

        t->rng.getrandom(t, x);

        for( dim = 0; dim < t->ndim; ++dim ) {
          creal pos = *x*NBINS;
          ccount ipos = (count)pos;
          creal prev = (ipos == 0) ? 0 : state_grid[dim][ipos - 1];
          creal diff = state_grid[dim][ipos] - prev; 
          *x++ = prev + (pos - ipos)*diff;
          *bin++ = ipos;
          weight *= diff*NBINS;
        }

        *w++ = weight;
      }

      DoSample(t, n, w, f, t->frame, state->niter + 1);

      bin = bins;
      w = t->frame;

      while( f < lastf ) {
        creal weight = *w++;
        Grid *m = &margsum[0][0];

        for( c = state->cumul; c < C; ++c ) {
          real wfun = weight*(*f++);
          if( wfun ) {
            c->sum += wfun;
            c->sqsum += wfun *= wfun;
            for( dim = 0; dim < t->ndim; ++dim )
              m[dim][bin[dim]] += wfun;
          }
          m += t->ndim;
        }

        bin += t->ndim;
      }
    }

    fail = 0;

    /* compute the integral and error values */

    for( c = state->cumul; c < C; ++c ) {
      real w = Weight(c->sum, c->sqsum, state->nsamples);
      real sigsq = 1/(c->weightsum += w);
      real avg = sigsq*(c->avgsum += w*c->sum);

      c->avg = LAST ? (sigsq = 1/w, c->sum) : avg;
      c->err = sqrt(sigsq);
      fail |= (c->err > MaxErr(c->avg));

      if( state->niter == 0 ) c->guess = c->sum;
      else {
        c->chisum += w *= c->sum - c->guess;
        c->chisqsum += w*c->sum;
      }
      c->chisq = c->chisqsum - avg*c->chisum;

      c->sum = c->sqsum = 0;
    }

    if( VERBOSE ) {
      char *oe = out + sprintf(out, "\n"
        "Iteration " COUNT ":  " NUMBER " integrand evaluations so far",
        state->niter + 1, t->neval);
      for( c = state->cumul, comp = 0; c < C; ++c )
        oe += sprintf(oe, "\n[" COUNT "] "
          REAL " +- " REAL "  \tchisq " REAL " (" COUNT " df)",
          ++comp, c->avg, c->err, c->chisq, state->niter);
      Print(out);
    }

    if( fail == 0 && t->neval >= t->mineval ) break;

    if( t->neval >= t->maxeval && !StateWriteTest(t) ) break;

    if( t->ncomp == 1 )
      for( dim = 0; dim < t->ndim; ++dim )
        RefineGrid(t, state_grid[dim], margsum[0][dim]);
    else {
      for( dim = 0; dim < t->ndim; ++dim ) {
        Grid wmargsum;
        Zap(wmargsum);
        for( comp = 0; comp < t->ncomp; ++comp ) {
          real w = state->cumul[comp].avg;
          if( w != 0 ) {
            creal *m = margsum[comp][dim];
            count bin;
            w = 1/Sq(w);
            for( bin = 0; bin < NBINS; ++bin )
              wmargsum[bin] += w*m[bin];
          }
        }
        RefineGrid(t, state_grid[dim], wmargsum);
      }
    }

    ++state->niter;
    state->nsamples += t->nincrease;

    if( StateWriteTest(t) ) {
      state->signature = StateSignature(t, 1);
      state->neval = t->neval;
      StateWriteOpen(t, fd) {
        StateWrite(fd, state, statesize);
      } StateWriteClose(t, fd);
      if( t->neval >= t->maxeval ) break;
    }
  }

  for( comp = 0; comp < t->ncomp; ++comp ) {
    cCumulants *c = &state->cumul[comp];
    integral[comp] = c->avg;
    error[comp] = c->err;
    prob[comp] = ChiSquare(c->chisq, state->niter);
  }

abort:
  PutGrid(t, state_grid);
  free(bins);
  WaitCores(t);
  FrameFree(t);

  StateRemove(t);

  return fail;
}

