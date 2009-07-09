


! Copyright (C) 2002-2007 S. Sagmeister J. K. Dewhurst, S. Sharma and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: portstate
! !INTERFACE:


subroutine portstate(act)
! !USES:
use modinput
  use ioarray
! !DESCRIPTION:
!   Toggle file format of {\tt STATE.OUT}. If tb2a is true an ASCII
!   file with the name {\tt STATE.xml} is generated and the data
!   from {\tt STATE.OUT} is transferred. If tb2a is false the conversion
!   goes in the other direction. Based upon the routines {\tt readstate}
!   and {\tt writestate}.
!
! !REVISION HISTORY:
!   Created 2007 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  integer, intent(in) :: act
  ! local variables
  logical::tb2a, exist, spinpol_
  integer::natmtot, is
  integer::version_(3), nspecies_, lmmaxvr_, nrmtmax_
  integer::natoms_(10000), ngrid_(3)
  integer::ngrtot_, ngvec_, ndmag_, nspinor_, ldapu_, lmmaxlu_
  ! allocatable arrays
  integer, allocatable :: nrmt_(:)
  real(8), allocatable :: spr_(:, :)
  real(8), allocatable :: rhomt_(:, :, :)
  real(8), allocatable :: rhoir_(:)
  real(8), allocatable :: vclmt_(:, :, :)
  real(8), allocatable :: vclir_(:)
  real(8), allocatable :: vxcmt_(:, :, :)
  real(8), allocatable :: vxcir_(:)
  real(8), allocatable :: veffmt_(:, :, :)
  real(8), allocatable :: veffir_(:)
  real(8), allocatable :: magmt_(:, :, :, :)
  real(8), allocatable :: magir_(:, :)
  real(8), allocatable :: bxcmt_(:, :, :, :)
  real(8), allocatable :: bxcir_(:, :)
  complex(8), allocatable :: veffig_(:)
  complex(8), allocatable :: vmatlu_(:, :, :, :, :)
  select case(act)
  case(1, 2, -1, -2)
  case default
     write(*, *)
     write(*, '("Error(portstate): unknown action: ", i6)') act
     write(*, *)
     stop
  end select
  tb2a=(act.eq.1).or.(act.eq.-1)
  if (tb2a) then
     open(50, file='STATE.OUT', action='READ', form='UNFORMATTED', &
	  status='OLD')
     if (act.eq.-1) then
	write(*, *)
	write(*, '("Information on STATE.OUT file:")')
	write(*, *)
     end if
     inquire(file='STATE.xml', exist=exist)
     if (exist.and.(act.eq.1)) then
	write(*, *)
	write(*, '("Error(portstate): not overwriting existent STATE.xml file")')
	write(*, *)
	stop
     end if
     if (act.eq.1) open(51, file = 'STATE.xml', action = 'WRITE', form = 'FORMATTED', &
	  status = 'replace')
     read(50) version_
     if (act.eq.-1) then
	write(*, '("version:", 3i8)') version_
	write(*, *)
	return
     end if
     read(50) spinpol_
     read(50) nspecies_
     read(50) lmmaxvr_
     read(50) nrmtmax_
     write(51, '(a)') '<?xml version="1.0"?>'
     write(51, '(a)') '<state>'
     write(51, '(a)') '<data name = "version" type = "integer" dimension = "1" &
	  &shape = "3">'
     call ioarr(un=51, ioa='write', arr1di=version_)
     write(51, '(a)') '</data>'
     write(51, '(a)') '<data name = "spinpol" type = "logical" dimension = "1" &
	  &shape = "1">'
     write(51, *) spinpol_
     write(51, '(a)') '</data>'
     write(51, '(a)') '<data name = "nspecies" type = "integer" dimension = "1" &
	  &shape = "1">'
     write(51, *) nspecies_
     write(51, '(a)') '</data>'
     write(51, '(a)') '<data name = "lmmaxvr" type = "integer" dimension = "1" &
	  &shape = "1">'
     write(51, *) lmmaxvr_
     write(51, '(a)') '</data>'
     write(51, '(a)') '<data name = "nrmtmax" type = "integer" dimension = "1" &
	  &shape = "1">'
     write(51, *) nrmtmax_
     write(51, '(a)') '</data>'
  else
     open(50, file='STATE.xml', action='READ', form='FORMATTED', &
	  status='OLD')
     if (act.eq.-2) then
	write(*, *)
	write(*, '("Information on STATE.xml file:")')
	write(*, *)
     end if
     inquire(file='STATE.OUT', exist=exist)
     if (exist.and.(act.eq.2)) then
	write(*, *)
	write(*, '("Error(portstate): not overwriting existent STATE.OUT file")')
	write(*, *)
	stop
     end if
     if (act.eq.2) open(51, file = 'STATE.OUT', action = 'WRITE', form = 'UNFORMATTED', &
	  status = 'replace')
     read(50, *)
     read(50, *)
     read(50, *)
     call ioarr(un=50, ioa='read', arr1di=version_)
     if (act.eq.-2) then
	write(*, '("version:", 3i8)') version_
	write(*, *)
	return
     end if
     read(50, *)
     read(50, *)
     read(50, *) spinpol_
     read(50, *)
     read(50, *)
     read(50, *) nspecies_
     read(50, *)
     read(50, *)
     read(50, *) lmmaxvr_
     read(50, *)
     read(50, *)
     read(50, *) nrmtmax_
     read(50, *)
     write(51) version_
     write(51) spinpol_
     write(51) nspecies_
     write(51) lmmaxvr_
     write(51) nrmtmax_
  end if
  allocate(spr_(nrmtmax_, nspecies_))
  allocate(nrmt_(nspecies_))
  if (tb2a) then
     natmtot=0
     do is=1, nspecies_
	read(50) natoms_(is)
	read(50) nrmt_(is)
	read(50) spr_(1:nrmt_(is), is)
	write(51, '(a)') '<data name = "natoms" type = "integer" dimension = "1" &
	     &shape = "1" index = "species" indexval = "'//trim(i2str(is))//'">'
	write(51, *) natoms_(is)
	write(51, '(a)') '</data>'
	write(51, '(a)') '<data name = "nrmt" type = "integer" dimension = "1" &
	     &shape = "1" index = "species" indexval = "'//trim(i2str(is))//'">'
	write(51, *) nrmt_(is)
	write(51, '(a)') '</data>'
	write(51, '(a)') '<data name = "spr" type = "real(8)" dimension = "1" &
	     &shape = "'//trim(i2str(nrmt_(is)))//'" index = "species" indexval = "'&
	     //trim(i2str(is))//'">'
     call ioarr(un=51, ioa='write', arr1dr=spr_(1:nrmt_(is), is))
	write(51, '(a)') '</data>'
	natmtot=natmtot+natoms_(is)
     end do
     read(50) ngrid_
     read(50) ngvec_
     read(50) ndmag_
     ! versions > 0.9.131
     if ((version_(1).gt.0).or.(version_(2).gt.9).or.(version_(3).gt.131)) then
	read(50) nspinor_
	read(50) ldapu_
	read(50) lmmaxlu_
     else
	nspinor_=1
	if (spinpol_) nspinor_=2
	ldapu_=0
	lmmaxlu_=0
     end if
     write(51, '(a)') '<data name = "ngrid" type = "integer" dimension = "1" &
	  &shape = "3">'
     call ioarr(un=51, ioa='write', arr1di=ngrid_)
     write(51, '(a)') '</data>'
     write(51, '(a)') '<data name = "ngvec" type = "integer" dimension = "1" &
	  &shape = "1">'
     write(51, *) ngvec_
     write(51, '(a)') '</data>'
     write(51, '(a)') '<data name = "ndmag" type = "integer" dimension = "1" &
	  &shape = "1">'
     write(51, *) ndmag_
     write(51, '(a)') '</data>'
     write(51, '(a)') '<data name = "nspinor" type = "integer" dimension = "1" &
	  &shape = "1">'
     write(51, *) nspinor_
     write(51, '(a)') '</data>'
     write(51, '(a)') '<data name = "ldapu" type = "integer" dimension = "1" &
	  &shape = "1">'
     write(51, *) ldapu_
     write(51, '(a)') '</data>'
	  write(51, '(a)') '<data name = "lmmaxlu" type = "integer" dimension = "1" &
	  &shape = "1">'
     write(51, *) lmmaxlu_
     write(51, '(a)') '</data>'
  else
     natmtot=0
     do is=1, nspecies_
	read(50, *)
	read(50, *) natoms_(is)
	read(50, *)
	read(50, *)
	read(50, *) nrmt_(is)
	read(50, *)
	read(50, *)
	call ioarr(un=50, ioa='read', arr1dr=spr_(1:nrmt_(is), is))
	read(50, *)
	write(51) natoms_(is)
	write(51) nrmt_(is)
	write(51) spr_(1:nrmt_(is), is)
	natmtot=natmtot+natoms_(is)
     end do
     read(50, *)
     call ioarr(un=50, ioa='read', arr1di=ngrid_)
     read(50, *)
     read(50, *)
     read(50, *) ngvec_
     read(50, *)
     read(50, *)
     read(50, *) ndmag_
     read(50, *)
     read(50, *)
     read(50, *) nspinor_
     read(50, *)
     read(50, *)
     read(50, *) ldapu_
     read(50, *)
     read(50, *)
     read(50, *) lmmaxlu_
     read(50, *)
     write(51) ngrid_
     write(51) ngvec_
     write(51) ndmag_
     write(51) nspinor_
     write(51) ldapu_
     write(51) lmmaxlu_
  end if
  ngrtot_=ngrid_(1)*ngrid_(2)*ngrid_(3)
  allocate(rhomt_(lmmaxvr_, nrmtmax_, natmtot))
  allocate(rhoir_(ngrtot_))
  allocate(vclmt_(lmmaxvr_, nrmtmax_, natmtot))
  allocate(vclir_(ngrtot_))
  allocate(vxcmt_(lmmaxvr_, nrmtmax_, natmtot))
  allocate(vxcir_(ngrtot_))
  allocate(veffmt_(lmmaxvr_, nrmtmax_, natmtot))
  allocate(veffir_(ngrtot_))
  allocate(veffig_(ngvec_))
  if (spinpol_) then
     allocate(magmt_(lmmaxvr_, nrmtmax_, natmtot, ndmag_))
     allocate(magir_(ngrtot_, ndmag_))
     allocate(bxcmt_(lmmaxvr_, nrmtmax_, natmtot, ndmag_))
     allocate(bxcir_(ngrtot_, ndmag_))
  end if
  if (ldapu_.ne.0) then
     allocate(vmatlu_(lmmaxlu_, lmmaxlu_, nspinor_, nspinor_, natmtot))
  end if
  if (tb2a) then
     ! read muffin-tin density
     read(50) rhomt_, rhoir_
     ! read Coulomb potential (spin independent)
     read(50) vclmt_, vclir_
     ! read exchange-correlation potential
     read(50) vxcmt_, vxcir_
     ! read effective potential
     read(50) veffmt_, veffir_, veffig_
     ! write the density
     write(51, '(a)') '<data name = "rhomt" type = "real(8)" &
	  &dimension = "3" shape = "'//&
	  trim(i2str(lmmaxvr_))//', '//&
	  trim(i2str(nrmtmax_))//', '//&
	  trim(i2str(natmtot))//'">'
     call ioarr(un=51, ioa='write', arr3dr=rhomt_)
     write(51, '(a)') '</data>'
     write(51, '(a)') '<data name = "rhoir" type = "real(8)" &
	  &dimension = "1" shape = "'//&
	  trim(i2str(ngrtot_))//'">'
     call ioarr(un=51, ioa='write', arr1dr=rhoir_)
     write(51, '(a)') '</data>'
     ! write the Coulomb potential
     write(51, '(a)') '<data name = "vclmt" type = "real(8)" &
	  &dimension = "3" shape = "'//&
	  trim(i2str(lmmaxvr_))//', '//&
	  trim(i2str(nrmtmax_))//', '//&
	  trim(i2str(natmtot))//'">'
     call ioarr(un=51, ioa='write', arr3dr=vclmt_)
     write(51, '(a)') '</data>'
     write(51, '(a)') '<data name = "vclir" type = "real(8)" &
	  &dimension = "1" shape = "'//&
	  trim(i2str(ngrtot_))//'">'
     call ioarr(un=51, ioa='write', arr1dr=vclir_)
     write(51, '(a)') '</data>'
     ! write the exchange-correlation potential
     write(51, '(a)') '<data name = "vxcmt" type = "real(8)" &
	  &dimension = "3" shape = "'//&
	  trim(i2str(lmmaxvr_))//', '//&
	  trim(i2str(nrmtmax_))//', '//&
	  trim(i2str(natmtot))//'">'
     call ioarr(un=51, ioa='write', arr3dr=vxcmt_)
     write(51, '(a)') '</data>'
     write(51, '(a)') '<data name = "vxcir" type = "real(8)" &
	  &dimension = "1" shape = "'//&
	  trim(i2str(ngrtot_))//'">'
     call ioarr(un=51, ioa='write', arr1dr=vxcir_)
     write(51, '(a)') '</data>'
     ! write the effective potential
     write(51, '(a)') '<data name = "veffmt" type = "real(8)" &
	  &dimension = "3" shape = "'//&
	  trim(i2str(lmmaxvr_))//', '//&
	  trim(i2str(nrmtmax_))//', '//&
	  trim(i2str(natmtot))//'">'
     call ioarr(un=51, ioa='write', arr3dr=veffmt_)
     write(51, '(a)') '</data>'
     write(51, '(a)') '<data name = "veffir" type = "real(8)" &
	  &dimension = "1" shape = "'//&
	  trim(i2str(ngrtot_))//'">'
     call ioarr(un=51, ioa='write', arr1dr=veffir_)
     write(51, '(a)') '</data>'
     write(51, '(a)') '<data name = "veffig" type = "complex(8)" &
	  &dimension = "1" shape = "'//&
	  trim(i2str(ngvec_))//'">'
     call ioarr(un=51, ioa='write', arr1dc=veffig_)
     write(51, '(a)') '</data>'
     if (spinpol_) then
        ! read magnetisation and effective field
	read(50) magmt_, magir_
	read(50) bxcmt_, bxcir_
        ! write the magnetisation and effective magnetic fields
	write(51, '(a)') '<data name = "magmt" type = "real(8)" &
	     &dimension = "4" shape = "'//&
	  trim(i2str(lmmaxvr_))//', '//&
	  trim(i2str(nrmtmax_))//', '//&
	  trim(i2str(natmtot))//', '//&
	  trim(i2str(ndmag_))//'">'
	call ioarr(un=51, ioa='write', arr4dr=magmt_)
	write(51, '(a)') '</data>'
	write(51, '(a)') '<data name = "magir" type = "real(8)" &
	     &dimension = "2" shape = "'//&
	     trim(i2str(ngrtot_))//', '//&
	     trim(i2str(ndmag_))//'">'
	call ioarr(un=51, ioa='write', arr2dr=magir_)
	write(51, '(a)') '</data>'
	write(51, '(a)') '<data name = "bxcmt" type = "real(8)" &
	     &dimension = "4" shape = "'//&
	  trim(i2str(lmmaxvr_))//', '//&
	  trim(i2str(nrmtmax_))//', '//&
	  trim(i2str(natmtot))//', '//&
	  trim(i2str(ndmag_))//'">'
	call ioarr(un=51, ioa='write', arr4dr=bxcmt_)
	write(51, '(a)') '</data>'
	write(51, '(a)') '<data name = "bxcir" type = "real(8)" &
	     &dimension = "2" shape = "'//&
	     trim(i2str(ngrtot_))//', '//&
	     trim(i2str(ndmag_))//'">'
	call ioarr(un=51, ioa='write', arr2dr=bxcir_)
	write(51, '(a)') '</data>'
     end if
     if (ldapu_.ne.0) then
        ! read the LDA+U potential matrix elements
	read(50) vmatlu_
        ! write the LDA+U potential matrix elements
	write(51, '(a)') '<data name = "vmatlu" type = "complex(8)" &
	     &dimension = "5" shape = "'//&
	     trim(i2str(lmmaxlu_))//', '//&
	     trim(i2str(lmmaxlu_))//', '//&
	     trim(i2str(nspinor_))//', '//&
	     trim(i2str(nspinor_))//', '//&
	     trim(i2str(natmtot))//'">'
	call ioarr(un=51, ioa='write', arr5dc=vmatlu_)
	write(51, '(a)') '</data>'
     end if
     write(51, '(a)') '</state>'
  else
     ! read muffin-tin density
     read(50, *)
     call ioarr(un=50, ioa='read', arr3dr=rhomt_)
     read(50, *)
     read(50, *)
     call ioarr(un=50, ioa='read', arr1dr=rhoir_)
     read(50, *)
     ! read Coulomb potential (spin independent)
     read(50, *)
     call ioarr(un=50, ioa='read', arr3dr=vclmt_)
     read(50, *)
     read(50, *)
     call ioarr(un=50, ioa='read', arr1dr=vclir_)
     read(50, *)
     ! read exchange-correlation potential
     read(50, *)
     call ioarr(un=50, ioa='read', arr3dr=vxcmt_)
     read(50, *)
     read(50, *)
     call ioarr(un=50, ioa='read', arr1dr=vxcir_)
     read(50, *)
     ! read effective potential
     read(50, *)
     call ioarr(un=50, ioa='read', arr3dr=veffmt_)
     read(50, *)
     read(50, *)
     call ioarr(un=50, ioa='read', arr1dr=veffir_)
     read(50, *)
     read(50, *)
     call ioarr(un=50, ioa='read', arr1dc=veffig_)
     read(50, *)
     ! write the density
     write(51) rhomt_, rhoir_
     ! write the Coulomb potential
     write(51) vclmt_, vclir_
     ! write the exchange-correlation potential
     write(51) vxcmt_, vxcir_
     ! write the effective potential
     write(51) veffmt_, veffir_, veffig_
     if (spinpol_) then
        ! read magnetisation and effective field
	read(50, *)
	call ioarr(un=50, ioa='read', arr4dr=magmt_)
	read(50, *)
	read(50, *)
	call ioarr(un=50, ioa='read', arr2dr=magir_)
	read(50, *)
	read(50, *)
	call ioarr(un=50, ioa='read', arr4dr=bxcmt_)
	read(50, *)
	read(50, *)
	call ioarr(un=50, ioa='read', arr2dr=bxcir_)
	read(50, *)
        ! write the magnetisation and effective magnetic fields
	write(51) magmt_, magir_
	write(51) bxcmt_, bxcir_
     end if
     if (ldapu_.ne.0) then
        ! read the LDA+U potential matrix elements
	read(50, *)
	call ioarr(un=50, ioa='read', arr5dc=vmatlu_)
	read(50, *)
        ! write the LDA+U potential matrix elements
	write(51) vmatlu_
     end if
  end if
  close(50)
  close(51)
  deallocate(nrmt_, spr_, rhomt_, rhoir_, vclmt_, vclir_)
  deallocate(vxcmt_, vxcir_, veffmt_, veffir_, veffig_)
  if (spinpol_) deallocate(magmt_, magir_, bxcmt_, bxcir_)
  if (tb2a) then
     write(*, *)
     write(*, '("Info(portstate): generated portable ASCII file &
	  &STATE.xml from STATE.OUT file")')
     write(*, *)
  else
     write(*, *)
     write(*, '("Info(portstate): generated STATE.OUT file from portable &
	  &ASCII file STATE.xml")')
     write(*, *)
  end if
contains
  character(256) function i2str(i)
    ! arguments
    integer, intent(in) :: i
    ! local variables
    character(1024) :: str
    write(str, *) i
    i2str=trim(adjustl(str))
  end function i2str
end subroutine portstate
!EOC
