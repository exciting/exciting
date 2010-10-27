/*
	Fluct.c
		compute the fluctuation in the left and right half
		this file is part of Suave
		last modified 9 Feb 05 th
*/


#if defined(HAVE_LONG_DOUBLE) && defined(HAVE_POWL)

typedef long double xdouble;
#define XDBL_MAX_EXP LDBL_MAX_EXP
#define XDBL_MAX LDBL_MAX
#define powx powl

#else

typedef double xdouble;
#define XDBL_MAX_EXP DBL_MAX_EXP
#define XDBL_MAX DBL_MAX
#define powx pow

#endif

typedef struct {
  xdouble fluct;
  number n;
} Var;

/*********************************************************************/

static void Fluct(cThis *t, Var *var,
  cBounds *b, creal *w, number n, ccount comp, creal avg, creal err)
{
  creal *x = w + n;
  creal *f = x + n*t->ndim + comp;
  creal flat = 2/3./t->flatness;
  creal max = ldexp(1., (int)((XDBL_MAX_EXP - 2)/t->flatness));
  creal norm = 1/(err*Max(fabs(avg), err));
  count nvar = 2*t->ndim;

  Clear(var, nvar);

  while( n-- ) {
    count dim;
    const xdouble ft =
      powx(Min(1 + fabs(*w++)*Sq(*f - avg)*norm, max), t->flatness);

    f += t->ncomp;

    for( dim = 0; dim < t->ndim; ++dim ) {
      Var *v = &var[2*dim + (*x++ >= b[dim].mid)];
      const xdouble f = v->fluct + ft;
      v->fluct = (f > XDBL_MAX/2) ? XDBL_MAX/2 : f;
      ++v->n;
    }
  }

  while( nvar-- ) {
    var->fluct = powx(var->fluct, flat);
    ++var;
  }
}

