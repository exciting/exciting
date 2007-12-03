
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
!   can then be passed directly to the mixing routine. See routine {\tt rfpack}.
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
integer idm,ias,lm1,lm2
integer ispn,jspn
n=0
call rfpack(tpack,n,1,veffmt,veffir,nu)
do idm=1,ndmag
  call rfpack(tpack,n,1,bxcmt(1,1,1,idm),bxcir(1,idm),nu)
end do
! pack the LDA+U potential if required
if (ldapu.ne.0) then
  do ias=1,natmtot
    do ispn=1,nspinor
      do jspn=1,nspinor
        do lm1=1,lmmaxlu
          do lm2=1,lmmaxlu
            n=n+1
            if (tpack) then
              nu(n)=dble(vmatlu(lm1,lm2,ispn,jspn,ias))
              n=n+1
              nu(n)=aimag(vmatlu(lm1,lm2,ispn,jspn,ias))
            else
              vmatlu(lm1,lm2,ispn,jspn,ias)=cmplx(nu(n),nu(n+1),8)
              n=n+1
            end if
          end do
        end do
      end do
    end do
  end do
end if
return
end subroutine
!EOC

