#include "keci.h"
#include "linalg.h"

void MultiKSpaceECI::get_k_space_eci(Array2d<Complex> *pft_eci, const rVector3d &k, Real x) {
  if (theonlyone) {
    theonlyone->get_k_space_eci(pft_eci,k,x);
    return;
  }
  zero_array(pft_eci);
  iVector2d size=pft_eci->get_size();
  Array2d<Complex> mat;
  LinkedListIterator<KSpaceECI> i(*this);
  for ( ; i; i++) {
    i->get_k_space_eci(&mat,k,x);
    sum(pft_eci,*pft_eci,mat);
  }
}

void MultiKSpaceECI::check_only_one(void) {
  if (LinkedList<KSpaceECI>::get_size()==1) {
    LinkedListIterator<KSpaceECI> i(*this);
    theonlyone=i;
  }
}

