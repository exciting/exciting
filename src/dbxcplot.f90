



! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine dbxcplot
use modmain
use modinput
implicit none
! local variables
integer::idm, is, ia, ias, ir
! allocatable arrays
real(8), allocatable :: rvfmt(:, :, :, :)
real(8), allocatable :: rvfir(:, :)
real(8), allocatable :: rfmt(:, :, :)
real(8), allocatable :: rfir(:)
real(8), allocatable :: grfmt(:, :, :, :)
real(8), allocatable :: grfir(:, :)
if (.not.associated(input%groundstate%spin)) then
  write(*, *)
  write(*, '("Error(dbxcplot): spin-unpolarised field is zero")')
  write(*, *)
  stop
end if
! initialise universal variables
call init0
! read magnetisation from file
call readstate
allocate(rvfmt(lmmaxvr, nrmtmax, natmtot, 3))
allocate(rvfir(ngrtot, 3))
allocate(rfmt(lmmaxvr, nrmtmax, natmtot))
allocate(rfir(ngrtot))
allocate(grfmt(lmmaxvr, nrmtmax, natmtot, 3))
allocate(grfir(ngrtot, 3))
if (ncmag) then
! non-collinear
  rvfmt(:, :, :, :)=bxcmt(:, :, :, :)
  rvfir(:, :)=bxcir(:, :)
else
! collinear
  rvfmt(:, :, :, 1:2)=0.d0
  rvfir(:, 1:2)=0.d0
  rvfmt(:, :, :, 3)=bxcmt(:, :, :, 1)
  rvfir(:, 3)=bxcir(:, 1)
end if
rfmt(:, :, :)=0.d0
rfir(:)=0.d0
do idm=1, 3
  call gradrf(rvfmt(:, :, :, idm), rvfir(:, idm), grfmt, grfir)
  do is=1, nspecies
    do ia=1, natoms(is)
      ias=idxas(ia, is)
      do ir=1, nrmt(is)
	rfmt(:, ir, ias)=rfmt(:, ir, ias)+grfmt(:, ir, ias, idm)
      end do
    end do
  end do
  rfir(:)=rfir(:)+grfir(:, idm)
end do
   if(associated(input%properties%gradmvecfield%plot1d)) then

  call plot1d("DBXC", 1, input%groundstate%lmaxvr, lmmaxvr, rfmt, rfir,input%properties%gradmvecfield%plot1d)

  write(*, *)
  write(*, '("Info(dbxcplot):")')
  write(*, '(" 1D divergence of exchange-correlation field written to &
   &DBXC1D.OUT")')
  write(*, '(" vertex location lines written to DBXCLINES.OUT")')
endif
if(associated(input%properties%gradmvecfield%plot2d)) then

  call plot2d("DBXC", 1, input%groundstate%lmaxvr, lmmaxvr, rfmt, rfir,input%properties%gradmvecfield%plot2d)

  write(*, '("Info(dbxcplot):")')
  write(*, '(" 2D divergence of exchange-correlation field written to &
   &DBXC2D.OUT")')
endif
if(associated(input%properties%gradmvecfield%plot3d)) then
  call plot3d("DBXC", 1, input%groundstate%lmaxvr, lmmaxvr, rfmt, rfir,input%properties%gradmvecfield%plot3d)
  write(*, '("Info(dbxcplot):")')
  write(*, '(" 3D divergence of exchange-correlation field written to &
   &DBXC3D.OUT")')
endif
write(*, *)
deallocate(rvfmt, rvfir, rfmt, rfir, grfmt, grfir)
return
end subroutine
