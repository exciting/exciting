
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: plot1d
! !INTERFACE:
subroutine plot1d(fnum1,fnum2,nf,lmax,ld,rfmt,rfir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   fnum1 : plot file number (in,integer)
!   fnum2 : vertex location file number (in,integer)
!   nf    : number of functions (in,integer)
!   lmax  : maximum angular momentum (in,integer)
!   ld    : leading dimension (in,integer)
!   rfmt  : real muffin-tin function (in,real(ld,nrmtmax,natmtot,nf))
!   rfir  : real intersitial function (in,real(ngrtot,nf))
! !DESCRIPTION:
!   Produces a 1D plot of the real functions contained in arrays {\tt rfmt} and
!   {\tt rfir} along the lines connecting the vertices in the global array
!   {\tt vvlp1d}. See routine {\tt rfarray}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: fnum1
integer, intent(in) :: fnum2
integer, intent(in) :: nf
integer, intent(in) :: lmax
integer, intent(in) :: ld
real(8), intent(in) :: rfmt(ld,nrmtmax,natmtot,nf)
real(8), intent(in) :: rfir(ngrtot,nf)
! local variables
integer i,ip,iv
real(8) fmin,fmax,t1
! allocatable arrays
real(8), allocatable :: fp(:,:)
if ((nf.lt.1).or.(nf.gt.4)) then
  write(*,*)
  write(*,'("Error(plot1d): invalid number of functions : ",I8)') nf
  write(*,*)
  stop
end if
allocate(fp(npp1d,nf))
! connect the plotting vertices
call connect(avec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
do i=1,nf
! evaluate function at each point
  call rfarray(lmax,ld,rfmt(:,:,:,i),rfir(:,i),npp1d,vplp1d,fp(:,i))
end do
fmin=fp(1,1)
fmax=fp(1,1)
do ip=1,npp1d
  do i=1,nf
    fmin=min(fmin,fp(ip,i))
    fmax=max(fmax,fp(ip,i))
  end do
! write the point distances and function to file
  write(fnum1,'(5G18.10)') dpp1d(ip),(fp(ip,i),i=1,nf)
end do
! write the vertex location lines
t1=0.5d0*(fmax-fmin)
do iv=1,nvp1d
  write(fnum2,'(2G18.10)') dvp1d(iv),fmax+t1
  write(fnum2,'(2G18.10)') dvp1d(iv),fmin-t1
  write(fnum2,'("     ")')
end do
deallocate(fp)
return
end subroutine
!EOC
