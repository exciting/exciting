#ifndef __MACHDEP_H__
#define __MACHDEP_H__

//#include <values.h>
#include <limits.h>
#include <float.h>
#define MAXINT INT_MAX

#ifndef MAXFLOAT
#define MAXFLOAT FLT_MAX
#endif

#define Real double
#ifdef OLD_COMPLEX
  #define Complex complex
#else
  #define Complex complex<double>
#endif
#define PATHSEP '/'
#include <unistd.h>
#include <stdlib.h>

#include "fixagg.h"

using namespace std;

#endif
