!BOP
!
! !ROUTINE: calcwmix0
!
! !INTERFACE:
      subroutine calcwmix0

! !DESCRIPTION:
!
!This subroutine calculate the matrix elements $\mathcal{W}^i_0$
!      

! !USES:

      use modmain
      use modgw
      use modmpi, only: rank

! !LOCAL VARIABLES:

      implicit none
      
      integer(4) :: mixind  ! Counter, runs over mixed basis functions
      integer(4) :: ia
      integer(4) :: ias
      integer(4) :: is
      integer(4) :: imix    ! Counter, runs over radial mixed functions
      integer(4) :: l1      ! Angular momentum quantum number of the
!                             mixed function
      integer(4) :: ippw    ! Counter, runs over interstitial mixed basis
!                             functions
      real(8)    :: fact
      real(8)    :: tstart,tend
      
! !INTRINSIC ROUTINES: 

      intrinsic sqrt

! !REVISION HISTORY:
!
! Created 11.02.05 by RGA
! Revisited June 2011 by DIN
!
!EOP
!BOC
      call cpu_time(tstart)
      
      fact=dsqrt(4.0d0*pi*vi)  ! To do not forget: Y_{00}=1/sqrt(4*pi)
      
      if(allocated(wi0))deallocate(wi0)
      allocate(wi0(matsiz))
      wi0(1:matsiz)=zzero

      mixind=0
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          do imix=1,nmix(ias)
            l1=bigl(ias,imix)
            if(l1.eq.0)then
              mixind=mixind+1
              wi0(mixind)=fact*rtl(ias,imix)
            else
              mixind=mixind+2*l1+1
            endif
          enddo ! im
        enddo ! ieq
      enddo ! iat
      
      do ippw=1,ngq(1)
        mixind=locmatsiz+ippw
        wi0(mixind)=mpwipw(ippw,1)
      enddo  

!      do mixind=1,matsiz
!        write(82,*)mixind,wi0(mixind)
!      enddo  

      call cpu_time(tend)
      if (rank == 0) call write_cputime(fgw,tend-tstart,'CALCWMIX0')
      
      end subroutine calcwmix0
!EOC                
