
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: fsmfield
! !INTERFACE:
subroutine fsmfield
! !USES:
use modmain
! !DESCRIPTION:
!   Updates the effective magnetic field, ${\bf B}_{\rm FSM}$, required for
!   fixing the spin moment to a given value, $\boldsymbol{\mu}_{\rm FSM}$. This
!   is done by adding a vector to the field which is proportional to the
!   difference between the moment calculated in the $i$th self-consistent loop
!   and the required moment:
!   $$ {\bf B}_{\rm FSM}^{i+1}={\bf B}_{\rm FSM}^i+\lambda\left(
!    \boldsymbol{\mu}^i-\boldsymbol{\mu}_{\rm FSM}\right), $$
!   where $\lambda$ is a scaling factor.
!
! !REVISION HISTORY:
!   Created March 2005 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ir,idm
real(8) v(3),t1
if ((.not.spinpol).or.(.not.fixspin)) return
if (ndmag.eq.3) then
  v(:)=momfix(:)
else
  v(1)=momfix(3)
end if
do idm=1,ndmag
  bfsmc(idm)=bfsmc(idm)+taufsm*(momtot(idm)-v(idm))
end do
t1=1.d0/y00
do idm=1,ndmag
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ir=1,nrmt(is)
        bxcmt(1,ir,ias,idm)=bxcmt(1,ir,ias,idm)+t1*bfsmc(idm)
      end do
    end do
  end do
  do ir=1,ngrtot
    bxcir(ir,idm)=bxcir(ir,idm)+bfsmc(idm)
  end do
end do
return
end subroutine
!EOC

