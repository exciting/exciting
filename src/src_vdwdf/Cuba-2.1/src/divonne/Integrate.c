/*
	Integrate.c
		partition the integration region until each region
		has approximately equal spread = 1/2 vol (max - min),
		then do a main integration over all regions
		this file is part of Divonne
		last modified 8 Jun 10 th
*/


#define INIDEPTH 3
#define DEPTH 5
#define POSTDEPTH 15

/*********************************************************************/

static int Integrate(This *t, real *integral, real *error, real *prob)
{
  TYPEDEFREGION;

  Totals totals[NCOMP];
  real nneed, weight;
  count dim, comp, iter, pass = 0, err, iregion;
  number nwant, nmin = INT_MAX;
  int fail = -99;

  if( VERBOSE > 1 ) {
    char s[512];
    sprintf(s, "Divonne input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  epsrel " REAL "\n  epsabs " REAL "\n"
      "  flags %d\n  seed %d\n"
      "  mineval " NUMBER "\n  maxeval " NUMBER "\n"
      "  key1 %d\n  key2 %d\n  key3 %d\n  maxpass " COUNT "\n"
      "  border " REAL "\n  maxchisq " REAL "\n  mindeviation " REAL "\n"
      "  ngiven " NUMBER "\n  nextra " NUMBER "\n",
      t->ndim, t->ncomp,
      t->epsrel, t->epsabs,
      t->flags, t->seed,
      t->mineval, t->maxeval,
      t->key1, t->key2, t->key3, t->maxpass,
      t->border.lower, t->maxchisq, t->mindeviation,
      t->ngiven, t->nextra);
    Print(s);
  }

  if( BadComponent(t) ) return -2;
  if( BadDimension(t, t->key1) ||
      BadDimension(t, t->key2) ||   
      ((t->key3 & -2) && BadDimension(t, t->key3)) ) return -1;

  t->neval_opt = t->neval_cut = 0;

  t->size = CHUNKSIZE;
  MemAlloc(t->voidregion, t->size*sizeof(Region));
  for( dim = 0; dim < t->ndim; ++dim ) {
    Bounds *b = &RegionPtr(0)->bounds[dim];
    b->lower = 0;
    b->upper = 1;
  }

  RuleIni(&t->rule7);
  RuleIni(&t->rule9);
  RuleIni(&t->rule11);
  RuleIni(&t->rule13);
  SamplesIni(&t->samples[0]);
  SamplesIni(&t->samples[1]);
  SamplesIni(&t->samples[2]);

  if( setjmp(t->abort) ) goto abort;

  t->epsabs = Max(t->epsabs, NOTZERO);

  /* Step 1: partition the integration region */

  if( VERBOSE ) Print("Partitioning phase:");

  if( IsSobol(t->key1) || IsSobol(t->key2) || IsSobol(t->key3) )
    IniRandom(t);

  SamplesLookup(t, &t->samples[0], t->key1,
    (number)47, (number)INT_MAX, (number)0);
  SamplesAlloc(t, &t->samples[0]);

  t->totals = totals;
  Zap(totals);
  t->phase = 1;

  Explore(t, 0, &t->samples[0], INIDEPTH, 1);

  for( iter = 1; ; ++iter ) {
    Totals *maxtot;
    count valid;

    for( comp = 0; comp < t->ncomp; ++comp ) {
      Totals *tot = &totals[comp];
      tot->avg = tot->spreadsq = 0;
      tot->spread = tot->secondspread = -INFTY;
    }

    for( iregion = 0; iregion < t->nregions; ++iregion ) {
      Region *region = RegionPtr(iregion);
      for( comp = 0; comp < t->ncomp; ++comp ) {
        cResult *r = &region->result[comp];
        Totals *tot = &totals[comp];
        tot->avg += r->avg;
        tot->spreadsq += Sq(r->spread);
        if( r->spread > tot->spread ) {
          tot->secondspread = tot->spread;
          tot->spread = r->spread;
          tot->iregion = iregion;
        }
        else if( r->spread > tot->secondspread )
          tot->secondspread = r->spread;
      }
    }

    maxtot = totals;
    valid = 0;
    for( comp = 0; comp < t->ncomp; ++comp ) {
      Totals *tot = &totals[comp];
      integral[comp] = tot->avg;
      valid += tot->avg == tot->avg;
      if( tot->spreadsq > maxtot->spreadsq ) maxtot = tot;
      tot->spread = sqrt(tot->spreadsq);
      error[comp] = tot->spread*t->samples[0].weight;
    }

    if( VERBOSE ) {
      char s[128 + 64*NCOMP], *p = s;

      p += sprintf(p, "\n"
        "Iteration " COUNT " (pass " COUNT "):  " COUNT " regions\n"
        NUMBER7 " integrand evaluations so far,\n"
        NUMBER7 " in optimizing regions,\n"
        NUMBER7 " in finding cuts",
        iter, pass, t->nregions, t->neval, t->neval_opt, t->neval_cut);

      for( comp = 0; comp < t->ncomp; ++comp )
        p += sprintf(p, "\n[" COUNT "] "
          REAL " +- " REAL,
          comp + 1, integral[comp], error[comp]);

      Print(s);
    }

    if( valid == 0 ) goto abort;	/* all NaNs */

    if( t->neval > t->maxeval ) break;

    nneed = maxtot->spread/MaxErr(maxtot->avg);
    if( nneed < MAXPRIME ) {
      cnumber n = t->neval + t->nregions*(number)ceil(nneed);
      if( n < nmin ) {
        nmin = n;
        pass = 0;
      }
      else if( ++pass > t->maxpass && n >= t->mineval ) break;
    }

    Split(t, maxtot->iregion, DEPTH);
  }

  /* Step 2: do a "full" integration on each region */

/* nneed = t->samples[0].neff + 1; */
  nneed = 2*t->samples[0].neff;
  for( comp = 0; comp < t->ncomp; ++comp ) {
    Totals *tot = &totals[comp];
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
    weight = t->samples[0].weight;
  }
  else {
    bool can_adjust = (t->key3 == 1 && t->samples[1].sampler != SampleRule &&
      (t->key2 < 0 || t->samples[1].neff < MAXPRIME));
    count df, nlimit;

    SamplesAlloc(t, &t->samples[1]);

    if( VERBOSE ) {
      char s[128];
      sprintf(s, "\nMain integration on " COUNT
        " regions with " NUMBER " samples per region.",
        t->nregions, t->samples[1].neff);
      Print(s);
    }

    ResClear(integral);
    ResClear(error);
    ResClear(prob);

    nlimit = t->maxeval - t->nregions*t->samples[1].n;
    df = 0;

    for( iregion = 0; iregion < t->nregions; ++iregion ) {
      Region *region = RegionPtr(iregion);
      char s[64*NDIM + 256*NCOMP], *p = s;
      int todo;

refine:
      t->phase = 2;
      t->samples[1].sampler(t, &t->samples[1], region->bounds, region->vol);

      if( can_adjust )
        for( comp = 0; comp < t->ncomp; ++comp )
          totals[comp].spreadsq -= Sq(region->result[comp].spread);

      nlimit += t->samples[1].n;
      todo = 0;

      for( comp = 0; comp < t->ncomp; ++comp ) {
        cResult *r = &region->result[comp];
        Totals *tot = &totals[comp];

        t->samples[0].avg[comp] = r->avg;
        t->samples[0].err[comp] = r->err;

        if( t->neval < nlimit ) {
          creal avg2 = t->samples[1].avg[comp];
          creal err2 = t->samples[1].err[comp];
          creal diffsq = Sq(avg2 - r->avg);

#define Var(s) Sq((s.err[comp] == 0) ? r->spread*s.weight : s.err[comp])

          if( err2*tot->nneed > r->spread ||
              diffsq > Max(t->maxchisq*(Var(t->samples[0]) + Var(t->samples[1])),
                           EPS*Sq(avg2)) ) {
            if( t->key3 && diffsq > tot->mindevsq ) {
              if( t->key3 == 1 ) {
                ccount xregion = t->nregions;

                if( VERBOSE > 2 ) Print("\nSplit");

                t->phase = 1;
                Explore(t, iregion, &t->samples[1], POSTDEPTH, 2);

                if( can_adjust ) {
                  number nnew;
                  count ireg, xreg;

                  for( ireg = iregion, xreg = xregion;
                       ireg < t->nregions; ireg = xreg++ ) {
                    cResult *result = RegionPtr(ireg)->result;
                    count c;
                    for( c = 0; c < t->ncomp; ++c )
                      totals[c].spreadsq += Sq(result[c].spread);
                  }

                  nnew = (tot->spreadsq/Sq(MARKMASK) > tot->maxerrsq) ?
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
                      char s[128];
                      sprintf(s, "Sampling remaining " COUNT
                        " regions with " NUMBER " points per region.",
                        t->nregions, t->samples[1].neff);
                      Print(s);
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

      if( can_adjust ) {
        for( comp = 0; comp < t->ncomp; ++comp )
          totals[comp].maxerrsq -=
            Sq(region->result[comp].spread*t->samples[1].weight);
      }

      switch( todo ) {
      case 1:	/* get spread right */
        Explore(t, iregion, &t->samples[1], 0, 2);
        break;

      case 3:	/* sample region again with more points */
        if( SamplesIniQ(&t->samples[2]) ) {
          SamplesLookup(t, &t->samples[2], t->key3,
            nwant, (number)INT_MAX, (number)0);
          SamplesAlloc(t, &t->samples[2]);
        }
        t->phase = 3;
        t->samples[2].sampler(t, &t->samples[2], region->bounds, region->vol);
        Explore(t, iregion, &t->samples[2], 0, 2);
        ++region->depth;	/* misused for df here */
        ++df;
      }

      ++region->depth;	/* misused for df here */

      if( VERBOSE > 2 ) {
        for( dim = 0; dim < t->ndim; ++dim ) {
          cBounds *b = &region->bounds[dim];
          p += sprintf(p,
            (dim == 0) ? "\nRegion (" REALF ") - (" REALF ")" :
                         "\n       (" REALF ") - (" REALF ")",
            b->lower, b->upper);
        }
      }

      for( comp = 0; comp < t->ncomp; ++comp ) {
        Result *r = &region->result[comp];

        creal x1 = t->samples[0].avg[comp];
        creal s1 = Var(t->samples[0]);
        creal x2 = t->samples[1].avg[comp];
        creal s2 = Var(t->samples[1]);
        creal r2 = (s1 == 0) ? Sq(t->samples[1].neff*t->samples[0].weight) : s2/s1;

        real norm = 1 + r2;
        real avg = x2 + r2*x1;
        real sigsq = s2;
        real chisq = Sq(x2 - x1);
        real chiden = s1 + s2;

        if( todo == 3 ) {
          creal x3 = t->samples[2].avg[comp];
          creal s3 = Var(t->samples[2]);
          creal r3 = (s2 == 0) ? Sq(t->samples[2].neff*t->samples[1].weight) : s3/s2;

          norm = 1 + r3*norm;
          avg = x3 + r3*avg;
          sigsq = s3;
          chisq = s1*Sq(x3 - x2) + s2*Sq(x3 - x1) + s3*chisq;
          chiden = s1*s2 + s3*chiden;
        }

        avg = LAST ? r->avg : (sigsq *= norm = 1/norm, avg*norm);
        if( chisq > EPS ) chisq /= Max(chiden, NOTZERO);

#define Out(s) s.avg[comp], r->spread*s.weight, s.err[comp]

        if( VERBOSE > 2 ) {
          p += sprintf(p, "\n[" COUNT "] "
            REAL " +- " REAL "(" REAL ")\n    "
            REAL " +- " REAL "(" REAL ")",
            comp + 1, Out(t->samples[0]), Out(t->samples[1]));
          if( todo == 3 ) p += sprintf(p, "\n    "
            REAL " +- " REAL "(" REAL ")",
            Out(t->samples[2]));
          p += sprintf(p, "  \tchisq " REAL, chisq);
        }

        integral[comp] += avg;
        error[comp] += sigsq;
        prob[comp] += chisq;

        r->avg = avg;
        r->spread = sqrt(sigsq);
        r->chisq = chisq;
      }

      if( VERBOSE > 2 ) Print(s);
    }

    for( comp = 0; comp < t->ncomp; ++comp )
      error[comp] = sqrt(error[comp]);

    df += t->nregions;

    if( VERBOSE > 2 ) {
      char s[16 + 128*NCOMP], *p = s;

      p += sprintf(p, "\nTotals:");

      for( comp = 0; comp < t->ncomp; ++comp )
        p += sprintf(p, "\n[" COUNT "] "
          REAL " +- " REAL "  \tchisq " REAL " (" COUNT " df)",
          comp + 1, integral[comp], error[comp], prob[comp], df);

      Print(s);
    }

    for( comp = 0; comp < t->ncomp; ++comp )
      prob[comp] = ChiSquare(prob[comp], df);

    weight = 1;
  }

#ifdef MLVERSION
  if( REGIONS ) {
    MLPutFunction(stdlink, "List", 2);
    MLPutFunction(stdlink, "List", t->nregions);
    for( iregion = 0; iregion < t->nregions; ++iregion ) {
      Region *region = RegionPtr(iregion);
      cBounds *b = region->bounds;
      real lower[NDIM], upper[NDIM];

      for( dim = 0; dim < t->ndim; ++dim ) {
        lower[dim] = b[dim].lower;
        upper[dim] = b[dim].upper;
      }

      MLPutFunction(stdlink, "Cuba`Divonne`region", 4);

      MLPutRealList(stdlink, lower, t->ndim);
      MLPutRealList(stdlink, upper, t->ndim);

      MLPutFunction(stdlink, "List", t->ncomp);
      for( comp = 0; comp < t->ncomp; ++comp ) {
        cResult *r = &region->result[comp];
        real res[] = {r->avg, r->spread*weight, r->chisq};
        MLPutRealList(stdlink, res, Elements(res));
      }

      MLPutInteger(stdlink, region->depth);  /* misused for df */
    }
  }
#endif

abort:
  SamplesFree(&t->samples[2]);
  SamplesFree(&t->samples[1]);
  SamplesFree(&t->samples[0]);
  RuleFree(&t->rule13);
  RuleFree(&t->rule11);
  RuleFree(&t->rule9);
  RuleFree(&t->rule7);

  free(t->voidregion);

  return fail;
}

