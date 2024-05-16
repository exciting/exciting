! Macro for to adjust the limits taken into account the intervals for which the QP corrections are computed, i.e. [ibgw, nbgw]
! I use macro to imitate the C++ behavior of function returning a range object. Also, it prevents code repetition
! that might lead to errors
! So this can directly used in array indexes (requires the C preprocessor)
! The first ifs are to guard the macro against multiple inclusions
use modgw, only: ibgw, nbgw
#ifndef QP_ADJUST_RANGE_DEFINED
#define QP_ADJUST_RANGE_DEFINED
#define QP_ADJUST_RANGE(lowband,upband) max(ibgw,lowband):min(nbgw,upband)
#endif