
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: plot3d
! !INTERFACE:
subroutine plot3d(fnum,nf,lmax,ld,rfmt,rfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   fnum : plot file number (in,integer)
!   nf   : number of functions (in,integer)
!   lmax : maximum angular momentum (in,integer)
!   ld   : leading dimension (in,integer)
!   rfmt : real muffin-tin function (in,real(ld,nrmtmax,natmtot,nf))
!   rfir : real intersitial function (in,real(ngrtot,nf))
! !DESCRIPTION:
!   Produces a 3D plot of the real functions contained in arrays {\tt rfmt} and
!   {\tt rfir} spanning the number of unit cells in the global array {\tt np3d}.
!   See routine {\tt rfarray}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: fnum
integer, intent(in) :: nf
integer, intent(in) :: lmax
integer, intent(in) :: ld
real(8), intent(in) :: rfmt(ld,nrmtmax,natmtot,nf)
real(8), intent(in) :: rfir(ngrtot,nf)
! local variables
integer i,ip,ip1,ip2,ip3,i1,i2,i3
real(8) vl(3),vc(3)
! allocatable arrays
real(8), allocatable :: vpl(:,:)
real(8), allocatable :: fp(:,:)
if ((nf.lt.1).or.(nf.gt.4)) then
  write(*,*)
  write(*,'("Error(plot3d): invalid number of functions : ",I8)') nf
  write(*,*)
  stop
end if
! allocate local arrays
allocate(vpl(3,np3d(1)*np3d(2)*np3d(3)))
allocate(fp(np3d(1)*np3d(2)*np3d(3),nf))
! generate 3D grid
do ip3=0,np3d(3)-1
  do ip2=0,np3d(2)-1
    do ip1=0,np3d(1)-1
      ip=ip3*np3d(1)*np3d(2)+ip2*np3d(1)+ip1+1
      vpl(1,ip)=dble(ip1)/dble(np3d(1))
      vpl(2,ip)=dble(ip2)/dble(np3d(2))
      vpl(3,ip)=dble(ip3)/dble(np3d(3))
    end do
  end do
end do
do i=1,nf
  call rfarray(lmax,ld,rfmt(1,1,1,i),rfir(1,i),ip,vpl,fp(1,i))
end do
write(fnum,'(3I6," : grid size")') nup3d(1)*np3d(1),nup3d(2)*np3d(2), &
 nup3d(3)*np3d(3)
do i3=0,nup3d(3)-1
  do ip3=0,np3d(3)-1
    do i2=0,nup3d(2)-1
      do ip2=0,np3d(2)-1
        do i1=0,nup3d(1)-1
          do ip1=0,np3d(1)-1
            ip=ip3*np3d(1)*np3d(2)+ip2*np3d(1)+ip1+1
            vl(1)=dble(i1)+dble(ip1)/dble(np3d(1))
            vl(2)=dble(i2)+dble(ip2)/dble(np3d(2))
            vl(3)=dble(i3)+dble(ip3)/dble(np3d(3))
            vc(:)=vl(1)*avec(:,1)+vl(2)*avec(:,2)+vl(3)*avec(:,3)
            write(fnum,'(7G18.10)') vc(1),vc(2),vc(3),(fp(ip,i),i=1,nf)
          end do
        end do
      end do
    end do
  end do
end do
deallocate(vpl,fp)
return
end subroutine
!EOC
