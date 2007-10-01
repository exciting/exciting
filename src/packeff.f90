
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: packeff
! !INTERFACE:
subroutine packeff(tpack,n,nu)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tpack : .true. for packing, .false. for unpacking (in,logical)
!   n     : total number of real values stored (out,integer)
!   nu    : packed potential (inout,real(*))
! !DESCRIPTION:
!   Packs/unpacks the muffin-tin and interstitial parts of the effective
!   potential and magnetic field into/from the single array {\tt nu}. This array
!   can then be passed directly to the mixing routine.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(out) :: n
real(8), intent(inout) :: nu(*)
! local variables
integer is,ia,ias,ir,lm,idm
n=0
! muffin-tin potential and field
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      do lm=1,lmmaxvr
        n=n+1
        if (tpack) then
          nu(n)=veffmt(lm,ir,ias)
        else
          veffmt(lm,ir,ias)=nu(n)
        end if
        do idm=1,ndmag
          n=n+1
          if (tpack) then
            nu(n)=bxcmt(lm,ir,ias,idm)
          else
            bxcmt(lm,ir,ias,idm)=nu(n)
          end if
        end do
      end do
    end do
  end do
end do
! interstitial potential and field
do ir=1,ngrtot
  n=n+1
  if (tpack) then
    nu(n)=veffir(ir)
  else
    veffir(ir)=nu(n)
  end if
  do idm=1,ndmag
    n=n+1
    if (tpack) then
      nu(n)=bxcir(ir,idm)
    else
      bxcir(ir,idm)=nu(n)
    end if
  end do
end do
return
end subroutine
!EOC

