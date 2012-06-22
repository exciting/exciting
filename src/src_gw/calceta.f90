!BOP
!
! !ROUTINE: calceta
!
! !INTERFACE:
      function calceta() result(neta)
!
! !DESCRIPTION:
!
! This function calculates the optimal value of $\eta $ for the lattice
! summations needed to obtain the structure constants.
!
! !USES:

      use modmain, only: bvec, pi
      use modgw,   only: avec

! !LOCAL VARIABLES:

      implicit none

      integer(4) :: i

      real(8) :: mlg
      real(8) :: mlr
      real(8) :: neta
      real(8), dimension (3)   :: lgbs
      real(8), dimension (3)   :: lrbs
!
! !EXTERNAL ROUTINES: 
!
! !INTRINSIC ROUTINES: 
!
      intrinsic dsqrt
      intrinsic minval
      intrinsic isign
!
! !REVISION HISTORY:
!
! Created 7th. January 2004 by RGA
! Last Modified: 30th. March 2004 by RGA
!
!EOP
!BOC
        do i=1,3
          lrbs(i)=dsqrt(avec(1,i)*avec(1,i)+avec(2,i)*avec(2,i)+&
     &            avec(3,i)*avec(3,i))
          lgbs(i)=dsqrt(bvec(1,i)*bvec(1,i)+bvec(2,i)*bvec(2,i)+&
     &           bvec(3,i)*bvec(3,i))
        enddo
      mlr=minval(lrbs)
      mlg=minval(lgbs)
      neta=dsqrt(2.0d0*mlr/mlg)
      end function calceta
!EOC
