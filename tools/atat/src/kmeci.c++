#include "kmeci.h"
#include "linalg.h"

void MultiKSpaceECI::get_k_space_eci(Array<Array<Array<Array<Array<Complex> > > > > *p_ft_eci, const Array<Real> &x) {
  LinkedListIterator<KSpaceECI> i(*this);
  for ( ; i; i++) {
    i->get_k_space_eci(p_ft_eci,x);
  }
}
