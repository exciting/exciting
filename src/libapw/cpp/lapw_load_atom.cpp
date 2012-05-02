#include "lapw.h"

/// load information about atom into libapw
extern "C" void FORTRAN(lapw_load_atom)(int *is_, int *ic_)
{
    int is = *is_ - 1;
    int ic = *ic_ - 1;
    
    lapw_global.atoms.push_back(new Atom(lapw_global.species[is], ic));
}
