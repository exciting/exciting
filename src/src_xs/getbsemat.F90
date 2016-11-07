! Copyright(C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getbsemat
! !INTERFACE:
subroutine getbsemat(fname, ikkp, n1, n2, zmat)
! !USES:
  use modmpi
  use modinput
  use m_getunit
! !DESCRIPTION:
!   This routine is used for reading the screened coulomb interaction
!   and exchange interaction from file. It works for ou,ou combinations.
!
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich)
!EOP
!BOC

  implicit none

  ! Arguments
  character(*), intent(in) :: fname
  integer, intent(in) :: ikkp, n1, n2
  complex(8), intent(out) :: zmat(n1, n2, n1, n2)

  ! Local variables
  integer :: un, reclen
  integer :: ikkp_, iknr_, jknr_, iq_, iqr_
  integer :: n1_, n2_, n3_, n4_
  complex(8), allocatable :: zm(:, :, :, :)

  call getunit(un)

  ! Get sizes of stored matrix
  inquire(iolength=reclen) ikkp_, iknr_, jknr_, iq_, iqr_, n1_, n2_, n3_, n4_

  open(unit=un, file=trim(fname), form='unformatted',&
    & action='read', access='direct', recl=reclen)

  read(un, rec=1) ikkp_, iknr_, jknr_, iq_, iqr_, n1_, n2_, n3_, n4_

  close(un)

  ! Check if requested size can be retrieved
  if((n1 .gt. n1_) .or. (n2 .gt. n2_)) then
    write(*,*)
    write(*, '("Error(getbsemat): requested matrix size out of range")')
    write(*, '(" requested size : ", 2i8)') n1, n2
    write(*, '(" stored size	 : ", 2i8)') n1_, n2_
    write(*,*)
    call terminate
  end if

  ! Read matrix
  allocate(zm(n1_, n2_, n3_, n4_))

  inquire(iolength=reclen) ikkp_, iknr_, jknr_, iq_, iqr_, n1_, n2_, n3_, n4_, zm

  open(unit=un, file=trim(fname), form='unformatted',&
    & action='read', access='direct', recl=reclen)
  
  read(un, rec=ikkp) ikkp_, iknr_, jknr_, iq_, iqr_, n1_, n2_, n3_, n4_, zm

  close(un)

  ! Check kkp-index
  if(ikkp .ne. ikkp_) then
    write(*,*)
    write(*, '("Error(getbsemat): inconsistent(k, kp)-index")')
    write(*, '(" requested		       : ", i8)') ikkp
    write(*, '(" stored at requested position : ", i8)') ikkp_
    write(*,*)
    call terminate
  end if

  ! Cut matrix
  ! In case the saved matrix zm(n1_,n2_,n1_,n2_) is larger
  ! then the requested slice of z(n1, n2, n1, n2) then
  if(input%xs%bse%beyond) then
    ! n1 <-> nu, n2 <-> no
    !  for n1 take elements form lowest index to highes
    !  for n2 take elements form highest index to lowest
    zmat(:, :, :, :) = zm(:n1, n2_-n2+1:, :n1,  n2_-n2+1:)
  else
    ! n1 <-> no, n2 <-> nu
    !  for n1 take elements form highest index to lowest
    !  for n2 take elements form lowest index to highes
    zmat(:, :, :, :) = zm(n1_-n1+1:, :n2, n1_-n1+1:, :n2)
  end if
  ! The behaviour stems from the fact that one refers to occupied states
  ! and the other to unoccupied ones.
  deallocate(zm)
end subroutine getbsemat
!EOC
