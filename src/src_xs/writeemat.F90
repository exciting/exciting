! Copyright (C) 2006-2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine writeemat
  use modinput
  use modmpi
  use mod_misc, only: task
  use mod_kpoint, only: nkpt
  use mod_qpoint, only: nqpt
  use mod_APW_LO, only: lolmax
  use modxs, only: unitout, xsgnt
  use m_xsgauntgen
  use m_findgntn0
  use m_filedel
  use m_genfilname

  implicit none

  ! local variables
  character(*), parameter :: thisnam = 'writeemat'
  integer(4) :: iq

  ! initialise universal variables
  call init0
  call init1
  call init2

  ! k-point parallelization for tddft
  if((task .ge. 300) .and. (task .le. 399)) then
    call genparidxran('k', nkpt)
  end if

  ! q-point parallelization for screening
  if((task .ge. 400) .and. (task .le. 499)) then
    call genparidxran('q', nqpt)
  end if

  ! write q-point set
  if(rank .eq. 0) call writeqpts

  ! read fermi energy from file
  call readfermi

  ! save variables for the gamma q-point
  call xssave0

  ! generate gaunt coefficients
  call xsgauntgen(max(input%groundstate%lmaxapw, lolmax), &
     & input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))

  ! find indices for non-zero gaunt coefficients
  call findgntn0(max(input%xs%lmaxapwwf, lolmax), &
     & max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)

  write(unitout, '(a, 3i5)') 'Info(' // thisnam // '): Gaunt coeff&
     &icients generated within lmax values:', &
     & input%groundstate%lmaxapw, input%xs%lmaxemat, &
     & input%groundstate%lmaxapw
  write (unitout, '(a, i5)') 'Info(' // thisnam // '): Number of q-&
     &points: ', nqpt
  call flushifc(unitout)

  ! loop over q-points
  do iq = 1, nqpt
    ! call for q-point
    call ematq(iq)
    
    write(unitout, '(a, i5)') 'Info(' // thisnam // '): Matrix el&
        &ements of the exponentials finished for q - point:', iq
    call flushifc (unitout)
  end do

  ! synchronize
  call barrier

  write(unitout, '(a)') "Info(" // trim (thisnam) // "): Matrix el&
     &ements of exponential expression finished"

  call findgntn0_clear
  call genfilname(setfilext=.true.)

end subroutine writeemat
