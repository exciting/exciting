
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
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
!   {\tt rfir} in the parallelepiped defined by the corner vertices in the
!   global array {\tt vclp3d}. See routine {\tt rfarray}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!   Modified, October 2008 (F. Bultmark, F. Cricchio, L. Nordstrom)
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
integer np,ip,ip1,ip2,ip3,i
real(8) v1(3),v2(3),v3(3)
real(8) t1,t2,t3
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
v1(:)=vclp3d(:,2)-vclp3d(:,1)
v2(:)=vclp3d(:,3)-vclp3d(:,1)
v3(:)=vclp3d(:,4)-vclp3d(:,1)
ip=0
do ip3=0,np3d(3)-1
  t3=dble(ip3)/dble(np3d(3))
  do ip2=0,np3d(2)-1
    t2=dble(ip2)/dble(np3d(2))
    do ip1=0,np3d(1)-1
      t1=dble(ip1)/dble(np3d(1))
      ip=ip+1
      vpl(:,ip)=t1*v1(:)+t2*v2(:)+t3*v3(:)+vclp3d(:,1)
    end do
  end do
end do
np=ip
! evaluate the functions at the grid points
do i=1,nf
  call rfarray(lmax,ld,rfmt(:,:,:,i),rfir(:,i),np,vpl,fp(:,i))
end do
! write functions to file
write(fnum,'(3I6," : grid size")') np3d(:)
do ip=1,np
  call r3mv(avec,vpl(:,ip),v1)
  write(fnum,'(7G18.10)') v1(:),(fp(ip,i),i=1,nf)
end do
deallocate(vpl,fp)
return
end subroutine
!EOC
