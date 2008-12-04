
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: plot2d
! !INTERFACE:
subroutine plot2d(fnum,nf,lmax,ld,rfmt,rfir)
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
!   Produces a 2D plot of the real functions contained in arrays {\tt rfmt} and
!   {\tt rfir} on the parallelogram defined by the corner vertices in the global
!   array {\tt vclp2d}. See routine {\tt rfarray}.
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
integer i,ip,ip1,ip2
real(8) vl1(3),vl2(3),vc1(3),vc2(3)
real(8) d1,d2,d12,t1,t2,t3,t4
! allocatable arrays
real(8), allocatable :: vpl(:,:)
real(8), allocatable :: fp(:,:)
if ((nf.lt.1).or.(nf.gt.4)) then
  write(*,*)
  write(*,'("Error(plot2d): invalid number of functions : ",I8)') nf
  write(*,*)
  stop
end if
! allocate local arrays
allocate(vpl(3,np2d(1)*np2d(2)))
allocate(fp(np2d(1)*np2d(2),nf))
! generate 2D grid
vl1(:)=vclp2d(:,2)-vclp2d(:,1)
vl2(:)=vclp2d(:,3)-vclp2d(:,1)
vc1(:)=vl1(1)*avec(:,1)+vl1(2)*avec(:,2)+vl1(3)*avec(:,3)
vc2(:)=vl2(1)*avec(:,1)+vl2(2)*avec(:,2)+vl2(3)*avec(:,3)
d1=sqrt(vc1(1)**2+vc1(2)**2+vc1(3)**2)
d2=sqrt(vc2(1)**2+vc2(2)**2+vc2(3)**2)
if ((d1.lt.epslat).or.(d2.lt.epslat)) then
  write(*,*)
  write(*,'("Error(plot2d): zero length plotting vectors")')
  write(*,*)
  stop
end if
d12=(vc1(1)*vc2(1)+vc1(2)*vc2(2)+vc1(3)*vc2(3))/(d1*d2)
ip=0
do ip2=0,np2d(2)-1
  do ip1=0,np2d(1)-1
    ip=ip+1
    t1=dble(ip1)/dble(np2d(1))
    t2=dble(ip2)/dble(np2d(2))
    vpl(:,ip)=t1*vl1(:)+t2*vl2(:)+vclp2d(:,1)
  end do
end do
! evaluate the functions at the grid points
do i=1,nf
  call rfarray(lmax,ld,rfmt(:,:,:,i),rfir(:,i),ip,vpl,fp(:,i))
end do
! write the functions to file
write(fnum,'(2I6," : grid size")') np2d(:)
ip=0
do ip2=0,np2d(2)-1
  do ip1=0,np2d(1)-1
    ip=ip+1
    t1=dble(ip1)/dble(np2d(1))
    t2=dble(ip2)/dble(np2d(2))
    t3=t1*d1+t2*d2*d12
    t4=t2*d2*sqrt(abs(1.d0-d12**2))
    write(fnum,'(6G18.10)') t3,t4,(fp(ip,i),i=1,nf)
  end do
end do
deallocate(vpl,fp)
return
end subroutine
!EOC
