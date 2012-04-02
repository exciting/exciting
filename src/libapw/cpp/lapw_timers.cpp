#include "lapw.h"

extern "C" void FORTRAN(lapw_timers)(void)
{
    timer::print();
}