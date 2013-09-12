/*
	Integrate.c
		partition the integration region until each region
		has approximately equal spread = 1/2 vol (max - min),
		then do a main integration over all regions
		this file is part of Divonne
		checkpointing by B. Chokoufe
		last modified 5 Aug 13 th
*/


typedef struct {
  signature_t signature;
  number neval, neval_opt, neval_cut;
  number nmin, nrand;
  count iregion, nregions;
  count phase, iter, pass, size;
  Totals totals[];
} State;

static int Integrate(This *t, real *integral, real *error, real *prob)
{
  StateDecl;
  csize_t statesize = sizeof(State) + NCOMP*sizeof(Totals);
  Sized(State, state, statesize);
  csize_t regionsize = RegionSize;
  Vector(char, out, 64*NDIM + 256*NCOMP);

  Totals *tot, *Tot = state->totals + t->ncomp;
  Bounds *b, *B;
  Result *res;
  count comp, err, iregion;
  number nwant;
  real nneed;
  ML_ONLY(number neff;)
  int fail;

  if( VERBOSE > 1 ) {
    sprintf(out, "Divonne input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  epsrel " REAL "\n  epsabs " REAL "\n"
      "  flags %d\n  seed %d\n"
      "  mineval " NUMBER "\n  maxeval " NUMBER "\n"
      "  key1 %d\n  key2 %d\n  key3 %d\n  maxpass " COUNT "\n"
      "  border " REAL "\n  maxchisq " REAL "\n  mindeviation " REAL "\n"
      "  ngiven " NUMBER "\n  nextra " NUMBER "\n"
      "  statefile \"%s\"",
      t->ndim, t->ncomp,
      t->epsrel, t->epsabs,
      t->flags, t->seed,
      t->mineval, t->maxeval,
      t->key1, t->key2, t->key3, t->maxpass,
      t->border.lower, t->maxchisq, t->mindeviation,
      t->ngiven, t->nextra,
      t->statefile);
    Print(out);
  }

  if( BadComponent(t) ) return -2;
  if( BadDimension(t, t->key1) ||
      BadDimension(t, t->key2) ||
      ((t->key3 & -2) && BadDimension(t, t->key3)) ) return -1;

  FORK_ONLY(t->nframe = 0;)

  AllocGiven(t);
  SamplesIni(&t->samples[0]);
  SamplesIni(&t->samples[1]);
  SamplesIni(&t->samples[2]);
  RuleAlloc(t);
  if( IsSobol(t->key1) | IsSobol(t->key2) | IsSobol(t->key3) )
    IniRandom(t);
  t->epsabs = Max(t->epsabs, NOTZERO);
  t->totals = state->totals;

  ForkCores(t);

  if( (fail = setjmp(t->abort)) ) goto abort;

  SamplesLookup(t, &t->samples[0], t->key1,
    (number)47, (number)INT_MAX, (number)0);
  SamplesAlloc(t, &t->samples[0]);

  StateSetup(t);

  if( StateReadTest(t) ) {
    StateReadOpen(t, fd) {
      if( read(fd, state, statesize) != statesize ||
          state->signature != StateSignature(t, 3) ) break;
      t->nregions = state->nregions;
      if( st.st_size != statesize + t->nregions*regionsize ) break;
      t->neval = state->neval;
      t->neval_opt = state->neval_opt;
      t->neval_cut = state->neval_cut;
      t->nrand = state->nrand;
      t->phase = state->phase;
      t->size = state->size;
      AllocRegions(t);
      StateRead(fd, t->region, t->nregions*regionsize);
    } StateReadClose(t, fd);

    if( IsSobol(t->key1) | IsSobol(t->key2) | IsSobol(t->key3) )
      t->rng.skiprandom(t, t->nrand);
  }

  if( ini ) {
#if MLVERSION
    t->neval = t->ngiven;
#else
    t->neval = 0;
#endif
    t->neval_opt = 0;
    t->neval_cut = 0;
    t->nrand = 0;

    t->size = CHUNKSIZE;
    AllocRegions(t);
    for( B = (b = t->region->bounds) + t->ndim; b < B; ++b ) {
      b->lower = 0;
      b->upper = 1;
    }
    t->nregions = 1;

    t->phase = 1;
    state->iter = 1;
    state->pass = 0;
    state->nmin = INT_MAX;
    state->iregion = 0;
    FClear(state->totals);
  }

  /* Step 1: partition the integration region */

  if( t->phase == 1 ) {
    if( VERBOSE ) Print("Partitioning phase:");

    if( ini ) Iterate(t, 0, INIDEPTH, 0, NULL);

    for( ; ; ++state->iter ) {
      Totals *maxtot;
      count valid;

      for( tot = state->totals; tot < Tot; ++tot ) {
        tot->avg = tot->spreadsq = 0;
        tot->spread = tot->secondspread = -INFTY;
      }

      for( iregion = 0; iregion < t->nregions; ++iregion )
        for( res = RegionResult(RegionPtr(iregion)), tot = state->totals;
             tot < Tot; ++res, ++tot ) {
          tot->avg += res->avg;
          tot->spreadsq += Sq(res->spread);
          if( res->spread > tot->spread ) {
            tot->secondspread = tot->spread;
            tot->spread = res->spread;
            tot->iregion = iregion;
          }
          else if( res->spread > tot->secondspread )
            tot->secondspread = res->spread;
        }

      valid = 0;
      for( maxtot = tot = state->totals, comp = 0; tot < Tot; ++tot, ++comp ) {
        integral[comp] = tot->avg;
        valid += tot->avg == tot->avg;
        if( tot->spreadsq > maxtot->spreadsq ) maxtot = tot;
        tot->spread = sqrt(tot->spreadsq);
        error[comp] = tot->spread/t->samples[0].neff;
      }

#define WriteState(t) \
if( StateWriteTest(t) ) { \
  StateWriteOpen(t, fd) { \
    state->signature = StateSignature(t, 3); \
    state->neval = t->neval; \
    state->neval_opt = t->neval_opt; \
    state->neval_cut = t->neval_cut; \
    state->nrand = t->nrand; \
    state->nregions = t->nregions; \
    state->phase = t->phase; \
    state->size = t->size; \
    StateWrite(fd, state, statesize); \
    StateWrite(fd, t->region, t->nregions*regionsize); \
  } StateWriteClose(t, fd); \
}

      WriteState(t);

      if( VERBOSE ) {
        char *oe = out + sprintf(out, "\n"
          "Iteration " COUNT " (pass " COUNT "):  " COUNT " regions\n"
          NUMBER7 " integrand evaluations so far,\n"
          NUMBER7 " in optimizing regions,\n"
          NUMBER7 " in finding cuts",
          state->iter, state->pass, t->nregions,
          t->neval, t->neval_opt, t->neval_cut);
        for( comp = 0; comp < t->ncomp; ++comp )
          oe += sprintf(oe, "\n[" COUNT "] "
            REAL " +- " REAL,
            comp + 1, integral[comp], error[comp]);
        Print(out);
      }

      if( valid == 0 ) goto abort;	/* all NaNs */

      if( t->neval > t->maxeval ) break;

      nneed = maxtot->spread/MaxErr(maxtot->avg);
      if( nneed < MAXPRIME ) {
        cnumber n = t->neval + t->nregions*(number)ceil(nneed);
        if( n < state->nmin ) {
          state->nmin = n;
          state->pass = 0;
        }
        else if( ++state->pass > t->maxpass && n >= t->mineval ) break;
      }

      Iterate(t, maxtot->iregion, DEPTH, -1, NULL);
    }
  }

  /* Step 2: do a "full" integration on each region */

/* nneed = t->samples[0].neff + 1; */
  nneed = 2*t->samples[0].neff;
  for( tot = state->totals; tot < Tot; ++tot ) {
    creal maxerr = MaxErr(tot->avg);
    tot->nneed = tot->spread/maxerr;
    nneed = Max(nneed, tot->nneed);
    tot->maxerrsq = Sq(maxerr);
    tot->mindevsq = tot->maxerrsq*Sq(t->mindeviation);
  }
  nwant = (number)Min(ceil(nneed), MARKMASK/40.);

  err = SamplesLookup(t, &t->samples[1], t->key2, nwant,
    (t->maxeval - t->neval)/t->nregions + 1, t->samples[0].n + 1);

  /* the number of points needed to reach the desired accuracy */
  fail = Unmark(err)*t->nregions;

  if( Marked(err) ) {
    if( VERBOSE ) Print("\nNot enough samples left for main integration.");
    for( comp = 0; comp < t->ncomp; ++comp )
      prob[comp] = -999;
    ML_ONLY(neff = t->samples[0].neff;)
  }
  else {
    bool can_adjust = (t->key3 == 1 &&
      t->samples[1].sampler != SampleRule &&
      (t->key2 < 0 || t->samples[1].neff < MAXPRIME));
    count df, nlimit;

    SamplesAlloc(t, &t->samples[1]);

    if( VERBOSE ) {
      sprintf(out, "\nMain integration on " COUNT
        " regions with " NUMBER " samples per region.",
        t->nregions, t->samples[1].neff);
      Print(out);
    }

    nlimit = t->maxeval - t->nregions*t->samples[1].n;
    df = 0;

#define CopyPhaseResults(f) \
  for( res = RegionResult(region), tot = state->totals; tot < Tot; ++res, ++tot ) { \
    PhaseResult *p = &tot->phase[f]; \
    p->avg = res->avg; \
    p->err = res->err; \
  }

#define Var2(f, r) Sq((r)->err ? (r)->err : res->spread/t->samples[f].neff)
#define Var(f) Var2(f, &tot->phase[f])

    while( state->iregion < t->nregions ) {
      Region *region;
      char *oe = out;
      int todo;

refine:
      region = RegionPtr(state->iregion);
      CopyPhaseResults(0);
      t->phase = 2;
      region->isamples = 1;
      t->samples[1].sampler(t, state->iregion);
      CopyPhaseResults(1);

      if( can_adjust )
        for( res = RegionResult(region), tot = state->totals;
             tot < Tot; ++res, ++tot )
          tot->spreadsq -= Sq(res->spread);

      nlimit += t->samples[1].n;
      todo = 0;

      for( res = RegionResult(region), tot = state->totals;
           tot < Tot; ++tot ) {
        if( t->neval < nlimit ) {
          creal avg2 = tot->phase[1].avg;
          creal diffsq = Sq(avg2 - tot->phase[0].avg);

          if( res->err*tot->nneed > res->spread ||
              diffsq > Max(t->maxchisq*(Var(0) + Var(1)), EPS*Sq(avg2)) ) {
            if( t->key3 && diffsq > tot->mindevsq ) {
              if( t->key3 == 1 ) {
                if( VERBOSE > 2 ) Print("\nSplit");
                t->phase = 1;
                Iterate(t, state->iregion, POSTDEPTH, 1, state->totals);

                if( can_adjust ) {
                  cnumber nnew = (tot->spreadsq/Sq(MARKMASK) > tot->maxerrsq) ?
                    MARKMASK :
                    (number)ceil(sqrt(tot->spreadsq/tot->maxerrsq));
                  if( nnew > nwant + nwant/64 ) {
                    ccount err = SamplesLookup(t, &t->samples[1], t->key2, nnew,
                      (t->maxeval - t->neval)/t->nregions + 1, t->samples[1].n);
                    fail += Unmark(err)*t->nregions;
                    nwant = nnew;
                    SamplesFree(&t->samples[1]);
                    SamplesAlloc(t, &t->samples[1]);

                    if( t->key2 > 0 && t->samples[1].neff >= MAXPRIME )
                      can_adjust = false;

                    if( VERBOSE > 2 ) {
                      sprintf(out, "Sampling remaining " COUNT
                        " regions with " NUMBER " points per region.",
                        t->nregions, t->samples[1].neff);
                      Print(out);
                    }
                  }
                }
                goto refine;
              }
              todo |= 3;
            }
            todo |= 1;
          }
        }
      }

      if( can_adjust )
        for( res = RegionResult(region), tot = state->totals;
             tot < Tot; ++res, ++tot )
          tot->maxerrsq -= Sq(res->spread/t->samples[1].neff);

      switch( todo ) {
      case 1:	/* get spread right */
        region->isamples = 1;
        ExploreSerial(t, state->iregion);
        break;

      case 3:	/* sample region again with more points */
        if( SamplesIniQ(&t->samples[2]) ) {
          SamplesLookup(t, &t->samples[2], t->key3,
            nwant, (number)INT_MAX, (number)0);
          SamplesAlloc(t, &t->samples[2]);
        }
        t->phase = 3;
        region->isamples = 2;
        t->samples[2].sampler(t, state->iregion);
        ExploreSerial(t, state->iregion);
        ++region->depth;	/* misused for df here */
        ++df;
      }

      if( VERBOSE > 2 ) {
        cchar *msg = "\nRegion (" REALF ") - (" REALF ")";
        for( B = (b = region->bounds) + t->ndim; b < B; ++b ) {
          oe += sprintf(oe, msg, b->lower, b->upper);
          msg = "\n       (" REALF ") - (" REALF ")";
        }
      }

      for( tot = state->totals, res = RegionResult(region), comp = 0;
           tot < Tot; ++tot, ++res ) {
        creal x1 = tot->phase[0].avg;
        creal v1 = Var(0);
        creal x2 = tot->phase[1].avg;
        creal v2 = Var(1);
        creal r2 = v1 ? v2/v1 :
          Sq(t->samples[1].neff/(real)t->samples[0].neff);

        real norm = 1 + r2;
        real avg = x2 + r2*x1;
        real sigsq = v2;
        real chisq = Sq(x2 - x1);
        real chiden = v1 + v2;

        if( todo == 3 ) {
          creal x3 = res->avg;
          creal v3 = Var2(2, res);
          creal r3 = v2 ? v3/v2 :
            Sq(t->samples[2].neff/(real)t->samples[1].neff);

          norm = 1 + r3*norm;
          avg = x3 + r3*avg;
          sigsq = v3;
          chisq = v1*Sq(x3 - x2) + v2*Sq(x3 - x1) + v3*chisq;
          chiden = v1*v2 + v3*chiden;
        }

        avg = LAST ? res->avg : (sigsq *= norm = 1/norm, avg*norm);
        if( chisq > EPS ) chisq /= Max(chiden, NOTZERO);

        if( VERBOSE > 2 ) {
#define Out2(f, r) (r)->avg, res->spread/t->samples[f].neff, (r)->err
#define Out(f) Out2(f, &tot->phase[f])
          oe += sprintf(oe, "\n[" COUNT "] "
            REAL " +- " REAL "(" REAL ")\n    "
            REAL " +- " REAL "(" REAL ")", ++comp, Out(0), Out(1));
          if( todo == 3 ) oe += sprintf(oe, "\n    "
            REAL " +- " REAL "(" REAL ")", Out2(2, res));
          oe += sprintf(oe, "  \tchisq " REAL, chisq);
        }

        tot->integral += avg;
        tot->sigsq += sigsq;
        tot->chisq += chisq;

        res->avg = avg;
        res->spread = sqrt(sigsq);
        res->chisq = chisq;
      }

      if( VERBOSE > 2 ) Print(out);
      ++state->iregion;

      WriteState(t);
    }

    df += t->nregions;

    for( tot = state->totals, comp = 0; tot < Tot; ++tot, ++comp ) {
      integral[comp] = tot->integral;
      error[comp] = sqrt(tot->sigsq);
      prob[comp] = ChiSquare(tot->chisq, df);
    }

    if( VERBOSE > 2 ) {
      char *oe = out + sprintf(out, "\nTotals:");
      for( tot = state->totals, comp = 0; tot < Tot; ++tot, ++comp )
        oe += sprintf(oe, "\n[" COUNT "] "
          REAL " +- " REAL "  \tchisq " REAL " (" COUNT " df)",
          comp + 1, integral[comp], error[comp], tot->chisq, df);
      Print(out);
    }

    ML_ONLY(neff = 1;)
  }

#ifdef MLVERSION
  if( REGIONS ) {
    Vector(real, bounds, t->ndim*2);

    MLPutFunction(stdlink, "List", 2);

    MLPutFunction(stdlink, "List", t->nregions);
    for( iregion = 0; iregion < t->nregions; ++iregion ) {
      Region *region = RegionPtr(iregion);
      cResult *Res;
      real *d = bounds;

      for( B = (b = region->bounds) + t->ndim; b < B; ++b ) {
        *d++ = b->lower;
        *d++ = b->upper;
      }

      MLPutFunction(stdlink, "Cuba`Divonne`region", 4);

      MLPutRealList(stdlink, bounds, 2*t->ndim);

      MLPutFunction(stdlink, "List", t->ncomp);
      for( Res = (res = RegionResult(region)) + t->ncomp; res < Res; ++res ) {
        real r[] = {res->avg, res->spread/neff, res->chisq};
        MLPutRealList(stdlink, r, Elements(r));
      }

      MLPutInteger(stdlink, region->depth + 1);  /* misused for df */
    }
  }
#endif

abort:
  WaitCores(t);
  FORK_ONLY(FrameFree(t, ShmRm(t));)

  RuleFree(t);
  SamplesFree(&t->samples[2]);
  SamplesFree(&t->samples[1]);
  SamplesFree(&t->samples[0]);
  free(t->region);
  free(t->xgiven);

  StateRemove(t);

  return fail;
}

