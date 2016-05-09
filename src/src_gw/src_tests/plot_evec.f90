!BOP
!
! !ROUTINE: plotevec
!
! !INTERFACE:
      subroutine plot_evec

! !DESCRIPTION:
!
! This subroutine calculates the real space representation of the
! (L)APW+lo basis functions in the line joining the two given atoms 
! for ploting
!
! !USES:

      use modinput
      use modmain
      use modgw
      use mod_rpath

! !INPUT PARAMETERS:

      implicit none

! !LOCAL VARIABLES:
      integer(4) :: ib
      integer(4) :: i
      integer(4) :: ia
      integer(4) :: ia1
      integer(4) :: ia2
      integer(4) :: irc
      integer(4) :: is
      integer(4) :: is1
      integer(4) :: is2
      integer(4) :: j
      integer(4) :: jrc
      integer(4) :: krc
      integer(4) :: np
      integer(4) :: ik
      
      character(len=64) :: filename

      real(8) :: irlen
      real(8) :: rd(3)
      real(8) :: rdlen
      real(8) :: rr(99)
      real(8), allocatable :: rs(:)
      
      complex(8), allocatable :: evec(:)

!
! !EXTERNAL ROUTINES: 
!
      external calcevecplot
      external ylm
! 
! !INTRINSIC ROUTINES: 
!
      intrinsic achar
      intrinsic dcos
      intrinsic dsin
      intrinsic mod
! 
! !REVISION HISTORY:
!
! Created: 9th. July 2004 by RGA
! Last modified: 15th. July 2004 by RGA
! Revisited 29.04.2011 by DIN
!
!EOP
!BOC
!
!     Set the name of the output file
!
      
      call boxmsg(6,'-','PLOTEVEC')
      
      ik=input%gw%iik
      write(*,*) 'Parameters:'
      write(*,*) 'k-point number (iik): ', ik
      write(*,*) 'lower bound for band number (ibmin): ', input%gw%ibmin
      write(*,*) 'upper bound for band number (ibmax): ', input%gw%ibmax
      write(*,*) 'atom 1 (at1): ', input%gw%at1
      write(*,*) 'atom 2 (at2): ', input%gw%at2
      write(*,*)

      if((input%gw%at1.lt.1).or.(input%gw%at1.gt.natmtot)) stop 'atom1 is wrong'
      if((input%gw%at2.lt.1).or.(input%gw%at2.gt.natmtot)) stop 'atom2 is wrong'

      call init_rpath(rpath,input%gw%at1,input%gw%at2)

!     Loop over eigenvectors
      do ib = input%gw%ibmin, input%gw%ibmax

         write(filename,5) ik, ib, input%gw%at1, input%gw%at2
         call str_strip(filename)
         open(unit=71,file=filename,status='unknown')

         allocate(evec(rpath%nptot))
         evec(:)=zzero
         call calcevecplot(ik,ib,evec)
!
!        write to file
!
         do i=1,np
           write(71,*) rs(i),real(evec(i)),aimag(evec(i))
         enddo
         deallocate(evec)
         
         close(71)
      end do ! ib

   5  format('evec-',i4,'-',i4,'-',i4,'-',i4,'.out')      
      return
      end subroutine
!EOC      
