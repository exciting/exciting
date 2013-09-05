!BOP
!
! !ROUTINE: init_freq
!
! !INTERFACE:
subroutine init_freq
      
! !DESCRIPTION:
!
! This subroutine generates the frequecy mesh for integration
!
! !USES:

      use modgw, only: fgw, freqs, freqmax, nomeg, wflag, womeg, debug

! !LOCAL VARIABLES:

      implicit none
      
      integer(4) :: i
      integer(4) :: n
      
      real(8), allocatable :: u(:)
      real(8), allocatable :: wu(:)
      

! !INTRINSIC ROUTINES: 

      intrinsic exp

! !EXTERNAL ROUTINES: 
      
      external gaulag
      external gauleg

! !REVISION HISTORY:
!
! Created 20.06.05 by RGA.
! Revisited 27.04.2011 by DIN
!
!EOP
!BOC
      call boxmsg(fgw,'-','Frequency grid')
      if(wflag.eq.3)then
        n=nomeg/2
        nomeg=n*2
      endif  
      if (allocated(freqs)) deallocate(freqs)
      allocate(freqs(1:nomeg))
      if (allocated(womeg)) deallocate(womeg)
      allocate(womeg(1:nomeg))
      
      select case(wflag) 

      case(1)    ! Equaly spaced mesh (for tests, not for integration)
        write(fgw,101)wflag
        write(fgw,102)'Equally spaced mesh'
        do i=1,nomeg
          freqs(i)=dble(i)*freqmax/dble(nomeg)
          womeg(i)=0.0d0
        enddo  
      
      case(2)    ! Grid for Gauss-Laguerre quadrature
        write(fgw,101)wflag
        write(fgw,102)'Grid for Gauss-Laguerre quadrature'
        allocate(wu(nomeg))
        call gaulag(freqs(1:nomeg),wu(1:nomeg),nomeg,-1.0d0)
        do i = 1, nomeg
          womeg(i) = wu(i)*exp(freqs(i))*freqs(i)
        enddo  
        deallocate(wu)

      case(3)    ! Grid for Gauss-Legendre quadrature from 0 to freqmax +
!                  Gauss-Legendre quadrature from freqmax to infinity
        write(fgw,101)wflag
        write(fgw,102)'Grid for double Gauss-Legendre quadrature,'
        write(fgw,102)'from 0 to freqmax and from freqmax to infinity'
        allocate(u(n))
        allocate(wu(n))
        call gauleg(-1.0d0,1.0d0,u,wu,n)
        do i = 1, n
          freqs(i)=freqmax*(u(i)+1.0d0)/2.0d0
          freqs(2*n-i+1)=2.0d0*freqmax/(u(i)+1.0d0)
          womeg(i)=freqmax*wu(i)/2.0d0
          womeg(2*n-i+1)=2.0d0*freqmax*wu(i)/((u(i)+1.0d0)*(u(i)+1.0d0))
        enddo ! i  
        deallocate(u)
        deallocate(wu)
          
      case(4)    ! Grid for Gauss-Legendre quadrature from 0 to omegamax
        write(fgw,101)wflag
        write(fgw,102)'Grid for Gauss-Legendre quadrature, from 0 to freqmax'
        call gauleg(0.0d0,freqmax,freqs(1:nomeg),womeg(1:nomeg),nomeg)

      end select
!
!     DEBUG
!
      if (debug) then
          write(fgw,103)nomeg,freqmax
          write(fgw,*)'      frequency       weight '
          do i=1,nomeg
            write(fgw,104)i,freqs(i),womeg(i)
          enddo
      end if
      
      call linmsg(fgw,'-','')

  101 format(' wflag=',i2,a)
  102 format(2x,a)
  103 format('number of frequencies:',i4,' freqmax=',f18.10)
  104 format(i4,1p,2g18.10)
      return

      end subroutine init_freq
!EOC
