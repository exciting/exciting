

! Copyright (C) 2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine mapxsparameters
	use modmain
use modinput
	use modxs
	implicit none
	input%groundstate%nosym=input%xs%nosym
	input%groundstate%ngkgrid(:)=input%xs%ngridk(:)
	input%groundstate%reducek=input%xs%reducek
	input%groundstate%vkloff(:)=input%xs%vkloff(:)
	ngridq(:)=input%xs%ngridq(:)
	if (associated(input%phonons))	input%phonons%reduceq=input%xs%reduceq
	input%groundstate%rgkmax=input%xs%rgkmax
	input%groundstate%swidth=input%xs%swidth
	input%groundstate%lmaxapw=input%xs%lmaxapw
	input%groundstate%lmaxmat=input%xs%lmaxmat
	input%groundstate%nempty=input%xs%nempty
end subroutine
