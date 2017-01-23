! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: df
! !INTERFACE:
!
subroutine df
! !USES:
  use modinput, only: input
  use mod_APW_LO, only: lolmax
  use mod_qpoint, only: nqpt
  use modxs, only: tscreen, xsgnt, nwdf, qpari,&
                   & qparf, unitout
  use modmpi
  use m_writegqpts
  use m_xsgauntgen
  use m_findgntn0
  use m_genfilname
! !DESCRIPTION:
!   Control routine for setting up the Kohn-Sham response function or the
!   microscopic dielectric function/matrix for all specified ${\bf q}$-points.
!   Can be run with MPI parallelization for ${\bf q}$-points.
!
! !REVISION HISTORY:
!   Created March 2006 (Sagmeister)
!EOP
!BOC

  implicit none
  
  ! Local variables
  character(*), parameter :: thisnam = 'df'
  character(256) :: filex
  integer :: iq

  if(.not. tscreen) call genfilname(setfilext=.true.)

  call init0
  ! Initialise universal variables
  call init1
  ! Save gamma-point variables
  call xssave0
  ! Initialize q-point set
  ! For 'screen' task 430 this sets up the reduced q-point set
  ! in mod_qpoint.
  call init2

  if(tscreen) then

    ! Generate gaunt coefficients, and store them in modxs:xsgnt
    call xsgauntgen(max(input%groundstate%lmaxapw, lolmax),&
      & input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))

    ! Find indices for non-zero gaunt coefficients, and store
    ! relevant maps in the module m_findgntn0
    call findgntn0(max(input%xs%lmaxapwwf, lolmax),&
      & max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)

  end if

  ! Read Fermi energy (coming form 'screen' this reads from EFERMI_SCR.OUT)
  if(input%xs%dogroundstate .ne. "fromscratch") call readfermi

  ! w-point parallelization for dielectric function
  if(tscreen) then
    ! Only one frequency if BSE is used
    nwdf = 1
    ! Use q point parallelization instead (reduced set)
    call genparidxran('q', nqpt)
  else
    call genparidxran('w', nwdf)
  end if

  ! Set type of band combinations: ({v,x},{x,c})- and ({x,c},{v,x})-combiantions
  input%xs%emattype = 1

  ! Write out q-points
  if(rank == 0) then
    call writeqpts
  end if

  ! Loop over q-points 
  qloop: do iq = qpari, qparf

    ! Default is "formfile"
    if(input%xs%dogroundstate .eq. "fromscratch") then 
      call genfilname(iqmt=iq, fileext=filex, setfilext=.true.)
      call readfermi
    else
      call genfilname(iq=iq, fileext=filex)
    end if

    ! Call for q-point
    call dfq(iq)

    if(tscreen) call writegqpts(iq, filex)

    write(unitout, '(a, i8)') 'Info(' // thisnam // '): Kohn Sham&
      & response function finished for q - point:', iq

  end do qloop

  ! Synchronize
  call barrier

  if((procs .gt. 1) .and. ( .not. tscreen)) then
    call dfgather
  end if

  call barrier

  if(rank == 0) then
    write(unitout, '(a)') "Info(" // trim(thisnam) // "): Kohn-Sham&
      & response function finished"
  end if

  if( .not. tscreen) call genfilname(setfilext=.true.)

  if(tscreen) call findgntn0_clear

end subroutine df
!EOC
