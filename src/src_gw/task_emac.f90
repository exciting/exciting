!BOP
!
! !ROUTINE: testemac
!
! !INTERFACE:
      subroutine task_emac
      
! !DESCRIPTION:
!
! This subroutine calculates the macroscopic dielectric function
!
! !USES:

      use modmain
      use modgw
      use modmpi
      
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: iq
      integer(4) :: iom
      
      real(8)    :: tstart,tend
!
! !EXTERNAL ROUTINES: 
!
     
      
! !INTRINSIC ROUTINES: 

      intrinsic cpu_time      
      
! !REVISION HISTORY:
!
! Created 16.09.2005 by RGA
! Revisited: June 2011 by DIN
!
!EOP
!BOC            
      call cpu_time(tstart)
!
!     Calculate the integration weights using the linearized tetrahedron method
!
      call kintw
!
!     Calculate the momentum matrix elements
!
      if(.not.input%gw%rpmat)then
        call calcpmat        ! <--- original (RGA's) version
        !call calcpmatgw     ! <--- modified xs (SAG's) version
      else
        write(fgw,*)'PMAT and PMATC are read from file'
      end if

!     Only \Gamma point
      iq=1
      Gamma=gammapoint(iq)
!
!     Interstitial mixed basis functions
!
      matsiz=locmatsiz+ngq(iq)
      call diagsgi(iq)
      call calcmpwipw(iq)
!
!     Calculate the bare coulomb potential matrix and its square root
!
      call calcbarcmb(iq)
      call setbarcev(input%gw%BareCoul%barcevtol)
!
!     Calculate the q-dependent integration weights
!   
      if(convflg.eq.0)then
        call qdepwsum(iq)
      else
        call qdepwtet(iq)
      endif  
!
!     Calculate the dielectric function matrix
!      
      if(allocated(epsilon))deallocate(epsilon)
      allocate(epsilon(matsizmax,matsizmax,nomeg))
      epsilon=zzero
      call calcepsilon(iq,0)
!
!     Inverse of the dielectric function
!      
      if(allocated(inveps))deallocate(inveps)
      allocate(inveps(matsizmax,matsizmax,nomeg))

      call calcinveps(iq)

      deallocate(epsilon)
      deallocate(inveps)
      deallocate(head)
      deallocate(epsw1,epsw2)

!------------------------------------------------------------------------
!     Test the calculated dielectric function
!------------------------------------------------------------------------      

      open(80,file='EPSMACRO.OUT',form='formatted',status='unknown')
      write(80,*)'### Macroscopic dielectric constant with and without local field'
      write(80,*)'### freqs(iom)      eps+LF             eps_00           eps^-1'
      
      write(fgw,10)
      write(fgw,11)
      write(80,10)
      write(80,11)

      do iom = 1, nomeg
      
        write(fgw,12)freqs(iom),emac(1:2,iom),1.d0/emac(2,iom)
        write(80, 12)freqs(iom),emac(1:2,iom),1.d0/emac(2,iom)
      
      end do

 10   format(//,"# emac with and without local field effects")
 11   format("# \omega(au)",6x,"\eps_M      ",6x,          &
     &                      6x,"\eps_M(NLF) ",6x,          &
     &                      6x,"\eps^{-1}   ")
 12   format(10F12.5)

      call cpu_time(tend)
      call write_cputime(fgw,tend-tstart,'TESTEPS')
      
      return
      end subroutine
!EOC
