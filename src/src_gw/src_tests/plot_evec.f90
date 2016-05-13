!BOP
!
! !ROUTINE: plotevec
!
! !INTERFACE:

subroutine plot_evec()

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
      integer(4) :: ib, ik, ir
      character(len=64) :: filename
      complex(8), allocatable :: evec(:)
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
      
      ik = input%gw%iik
      write(*,*) 'Parameters:'
      write(*,*) 'k-point number (iik): ', ik
      write(*,*) 'lower bound for band number (ibmin): ', input%gw%ibmin
      write(*,*) 'upper bound for band number (ibmax): ', input%gw%ibmax
      write(*,*) 'atom 1 (at1): ', input%gw%at1
      write(*,*) 'atom 2 (at2): ', input%gw%at2
      write(*,*)

      if ((input%gw%at1.lt.1).or.(input%gw%at1.gt.natmtot)) stop 'atom1 is wrong'
      if ((input%gw%at2.lt.1).or.(input%gw%at2.gt.natmtot)) stop 'atom2 is wrong'

      call init_rpath(rpath,input%gw%at1,input%gw%at2)

!     Loop over eigenvectors
      do ib = input%gw%ibmin, input%gw%ibmax

         write(filename,'("evec-",i4,"-",i4,"-",i4,"-",i4,".out")') &
         &     ik, ib, input%gw%at1, input%gw%at2
         call str_strip(filename)
         open(unit=71,file=filename,status='unknown')

         allocate(evec(rpath%nptot))
         evec(:) = zzero

         call calcevecplot(ik,ib,evec)

         do ir = 1, rpath%nptot
           write(71,*) rpath%r(ir,1), real(evec(ir)), aimag(evec(ir))
         enddo

         deallocate(evec)
         
         close(71)
      end do ! ib

      return
end subroutine
!EOC      
