!BOP
! !ROUTINE: diagsgi
! !INTERFACE:
      subroutine diagsgi(iq)

! !DESCRIPTION:
!
! This subroutine generates the overlap matrix of the product functions in the
! interstitial region and diagonalizes it. 
! The output is the matrix $S_{\vec{G}i}$.
!
! !USES:

      use modmain, only: ivgig, ivg
      use modgw
      use modmpi,  only: rank
      
! !INPUT PARAMETERS:
       
      implicit none

      integer(4), intent(in) :: iq      

! !LOCAL VARIABLES:

      integer(4) :: ig
      integer(4) :: igq   ! Counter: Runs over igq's
      integer(4) :: jgq   ! Counter: Runs over igq's
      integer(4) :: info  ! Output status of the diagonalization routine.
      integer(4) :: lwork ! Size of the workspace array work.
      integer(4), dimension(3) :: iig    ! integer coordinates of G-G'

      real(8), allocatable :: rwork(:)   ! Workspace array for the
      real(8), allocatable :: epsipw(:)  ! the eigenvalues of sgi
                                         ! diagonalization routine.
      complex(8) :: cfact                ! Normalization factor
      complex(8), allocatable :: work(:) ! Workspace array for the
                                         ! diagonalization routine.

      real(8)   :: tstart, tend

! !EXTERNAL ROUTINES: 

      external zheev      ! Lapack diagonalization subroutine

! !INTRINSIC ROUTINES: 

      intrinsic dabs
      intrinsic dsqrt
!
! !REVISION HISTORY:
!
! Created Dec. 2003. by RGA
! Last Modification May. 2006 by RGA
! Revisited 10.05.2011 by DIN
!
!EOP
!BOC
      
      call cpu_time(tstart)
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'

!     Allocate the working space needed by the diagonalizing subroutine
      lwork = 2*ngq(iq)-1
      allocate(epsipw(ngq(iq)))
      allocate(work(lwork))
      allocate(rwork(3*ngq(iq)-2))
      if(allocated(sgi)) deallocate(sgi)
      allocate(sgi(ngq(iq),ngq(iq)))
!
!     Calculate the overlap matrix between product plane waves:
!
      sgi=0.0d0
      do igq=1,ngq(iq)
        ig=ivgig(0,0,0)
        !sgi(igq,igq)=conjg(cfunig(ig))
        sgi(igq,igq)=ipwint(ig)
!       Non diagonal elements:
        do jgq=igq+1,ngq(iq)
          iig(:)=ivg(:,igqig(igq,iq))-ivg(:,igqig(jgq,iq))
          ig=ivgig(iig(1),iig(2),iig(3))
          !sgi(igq,jgq)=conjg(cfunig(ig))
          sgi(igq,jgq)=ipwint(ig)
          sgi(jgq,igq)=conjg(sgi(igq,jgq))
        enddo ! jgq
      enddo ! igq
      
      if(DEBUG)then
        write(55,*) "### sgi-0 ###"
        do jgq=1,ngq(iq),ngq(iq)/10
          do igq=1,ngq(iq),ngq(iq)/10
            write(55,'(2i5,4f12.6)') igq,jgq,sgi(igq,jgq)
          enddo
        enddo
      endif 
     
!
!     Diagonalize sgi:
!
      call zheev('V','U',ngq(iq),sgi,ngq(iq),epsipw,work,lwork,rwork,info)
      if(info.ne.0) stop "diagsgi: Fail in calling zheev"

      if(DEBUG)then
        write(55,*) "### sgi-1 ###"
        do jgq=1,ngq(iq),ngq(iq)/10
          do igq=1,ngq(iq),ngq(iq)/10
            write(55,'(2i5,4f12.6)') igq,jgq,sgi(igq,jgq)
          enddo
        enddo
        write(55,*)
        write(55,*) "### epsipw ###"
        do igq=1,ngq(iq)
          write(55,'(i5,2f12.6)') igq,epsipw(igq)
        enddo 
      endif
!
!     Normalize sgi:
!
      do igq=1,ngq(iq)
        cfact=cmplx(1.0d0/sqrt(dabs(epsipw(igq))),0.0d0,8)
        sgi(:,igq)=cfact*sgi(:,igq)
      enddo
      
      deallocate(rwork)
      deallocate(work)
      deallocate(epsipw)

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      if (rank==0) call write_cputime(fgw,tend-tstart, 'DIAGSGI')
      
      return
      end subroutine diagsgi
!EOC

