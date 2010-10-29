#if SINGLE_PRECISION
#  define FLOAT float
#  define POW   powf
#  define LOG   logf
#  define ASINH asinhf
#  define ABS   fabsf
#  define XC(x) xc_s_ ## x
#  define XC_U(X) XC_S_ ## X
#  define FLOAT_EPSILON FLT_EPSILON
#  define FLOAT_MIN FLT_MIN
#  define FLOAT_MAX FLT_MAX
#else
#  define FLOAT double
#  define POW   pow
#  define LOG   log
#  define ASINH asinh
#  define ABS   fabs
#  define XC(x) xc_ ## x
#  define XC_U(X) XC_ ## X
#  define FLOAT_EPSILON DBL_EPSILON
#  define FLOAT_MIN DBL_MIN
#  define FLOAT_MAX DBL_MAX
#endif

#define XC_FC_FUNC2(a,b) FC_FUNC_(a,b) 
#define XC_FC_FUNC(a,b) XC_FC_FUNC2(XC(a), XC_U(b))
