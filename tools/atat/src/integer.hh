#ifndef _INTEGER__H_
#define _INTEGER__H_
 
#include <math.h>
#include "misc.h"

extern Real zero_tolerance;

inline int near_zero(Real x) {
  return (fabs(x)<zero_tolerance);
}

inline int is_int(Real x) {
  return (fabs(x-rint(x)) < zero_tolerance);
}

int least_common_multiple(int a, int b);
int integer_ratio(int *p, int *q, Real x, Real epsilon);
int factorial(int n);

#endif
