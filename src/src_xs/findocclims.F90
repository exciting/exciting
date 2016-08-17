! Copyright(C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findocclims
! !INTERFACE:
subroutine findocclims(iq, iocc0, iocc, iunocc0, iunocc, io0, io, iu0, iu)
! !USES:
  use mod_kpoint, only: nkpt, vkl
  use mod_eigenvalue_occupancy, only: nstsv, occsv, evalsv, occmax,&
                                    & efermi
  use modinput, only: input
  use modxs, only: evalsv0, ikmapikq, occsv0, evlmin,&
                 & evlmax, evlmincut, evlmaxcut, evlhpo,&
                 & evllpu, ksgap, nstocc0,&
                 & vkl0, nstunocc0, unitout
  use m_genfilname
! !DESCRIPTION:
! For a given q-point index the routine inspects the occupancies and 
! second variational eigenvalues of the corresponding k+q and k set to determin
! the highest/lowest (at least partially) occupied/unoccupied 
! state index for all k+q and k points and the highest/lowest occupied/unoccupied
! state index over all k+q and k points.
! 
! !REVISION HISTORY:
!   Added description schema. And rudimentary description. (Aurich)
!EOP
!BOC

!NOTE: iunocc is always set to iunocc0 (lowest unoccupied state over all k points)

  implicit none

  ! Arguments
  integer, intent(in) :: iq
  integer, intent(out) :: iocc0, iocc, iunocc0, iunocc
  integer, intent(out) :: io0(nkpt), io(nkpt), iu0(nkpt), iu(nkpt)

  ! Local variables
  integer :: ik, ikq, i0, i
  logical :: t

  t = allocated(evalsv0)
  if( .not. t) allocate(evalsv0(nstsv, nkpt))

  k: do ik = 1, nkpt

    ! k+q-point set
    ikq = ik
    if(iq .ne. 0) ikq = ikmapikq(ik, iq)

    ! Get occupancies and eigenvalues
    call getoccsv(vkl(1, ikq), occsv(1, ikq))
    call getevalsv(vkl(1, ikq), evalsv(1, ikq))

    ! Count how, from low to high, many states at k+q are (partially) occupied,
    ! in the sense that they contribute to the electron density.
    ! epsocc defaults to 10^-8
    do i = 1, nstsv
      if(occsv(i, ikq) .lt. input%groundstate%epsocc) exit
    end do
    io(ik) = i - 1

    ! Count how many states are (partially) unoccupied. Count
    ! from highest state and stop when the states become fully
    ! occupied.
    do i = nstsv, 1, - 1
      if(occsv(i, ikq) .gt. (occmax-input%groundstate%epsocc)) exit
    end do
    iu(ik) = i + 1

    ! Do the same counting for the k point without q 
    if(iq .ne. 0) then
      ! k-point set (q=0)
      call getoccsv0(vkl0(1, ik), occsv0(1, ik))
      call getevalsv0(vkl0(1, ik), evalsv0(1, ik))
      do i0 = 1, nstsv
        if(occsv0(i0, ik) .lt. input%groundstate%epsocc) exit
      end do
      io0(ik) = i0 - 1
      do i0 = nstsv, 1, - 1
        if(occsv0(i0, ik) .gt. (occmax-input%groundstate%epsocc)) exit
      end do
      iu0(ik) = i0 + 1
    else
      io0(ik) = io(ik)
      iu0(ik) = iu(ik)
    end if

  end do k

  if(iq .ne. 0) then

    ! Lowest and highest valence energy 
    evlmin = min(minval(evalsv(1, :)), minval(evalsv0(1, :)))
    evlmax = max(maxval(evalsv(nstsv, :)), maxval(evalsv0(nstsv, :)))

    ! Lower and higher cutoff valence energy
    evlmincut = max(maxval(evalsv(1, :)), maxval(evalsv0(1, :)))
    evlmaxcut = min(minval(evalsv(nstsv, :)), minval(evalsv0(nstsv, :)))

  else

    ! Lowest and highest valence energy
    evlmin = minval(evalsv(1, :))
    evlmax = maxval(evalsv(nstsv, :))

    ! Lower and higher cutoff valence energy
    evlmincut = maxval(evalsv(1, :))
    evlmaxcut = minval(evalsv(nstsv, :))

  end if

  ! Overall highest (partially) occupied state
  iocc0 = maxval(io0)
  iocc = maxval(io)

  ! Overall lowest (partially) unoccupied state
  iunocc0 = minval(iu0)
  iunocc = minval(iu)

  ! The maximum/minimum value is used since a shifted (k+q)-mesh which is not
  ! commensurate can cause partially occupied states that are absent for the
  ! k-mesh
  iocc0 = max(iocc0, iocc)
  iocc = iocc0
  iunocc0 = min(iunocc0, iunocc)
  iunocc = iunocc0

  ! Determine if system has a gap in energy
  if(iq .ne. 0) then
    ! Highest (partially) occupied state energy
    evlhpo = max(maxval(evalsv(iocc0, :)), maxval(evalsv0(iocc0, :)))
    ! Lowest (partially) unoccupied state energy
    evllpu = min(minval(evalsv(iunocc0, :)), minval(evalsv0(iunocc0, :)))
  else
    ! Highest (partially) occupied state energy
    evlhpo = maxval(evalsv(iocc0, :))
    ! Lowest (partially) unoccupied state energy
    evllpu = minval(evalsv(iunocc0, :))
  end if

  ! Determine if system has a gap in energy
  ksgap = evlhpo .lt. efermi

  ! Assign nstocc0 and nstunocc0
  nstocc0 = iocc0
  nstunocc0 = nstsv - nstocc0
  if((iocc0 .ge. iunocc) .or. (iocc .ge. iunocc0)) then
    write(unitout, '(a)') 'Info(findocclims): Partially occupied states present'
  end if
  if(ksgap) then
    write(unitout, '(a)') 'Info(findocclims): System has kohn-sham gap'
  else
    write(unitout, '(a)') 'Info(findocclims): No kohn-sham gap found'
  end if

  ! Debug output
  if(input%xs%dbglev .gt. 0) then
    write(*, '(a)') 'debug(findocclims):'
    write(*, '(a)') ' iocc0, iocc, iunocc0, iunocc below:'
    write(*, '(4i8)') iocc0, iocc, iunocc0, iunocc
    write(*, '(a)') ' ik, io0, iu, diff, io, iu0, diff below:'
    do ik = 1, nkpt
       write(*, '(7i8)') ik, io0(ik), iu(ik), iu(ik) - io0(ik),&
         & io(ik), iu0(ik), iu0(ik) - io(ik)
    end do
    write(*,*)
  end if

  if( .not. t) deallocate(evalsv0)
end subroutine findocclims
!EOC
