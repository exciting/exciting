! Copyright(C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findocclims
! !INTERFACE:
subroutine findocclims(iq, iocc0, iocc, iunocc0, iunocc, io0, io, iu0, iu)
! !USES:
  use mod_constants, only: h2ev
  use mod_kpoint, only: nkpt, vkl
  use mod_eigenvalue_occupancy, only: nstsv, occsv, evalsv, occmax,&
                                    & efermi
  use modinput, only: input
  use modxs, only: evalsv0, ikmapikq, occsv0, evlmin,&
                 & evlmax, evlmincut, evlmaxcut, evlhpo,&
                 & evllpu, ksgap, nstocc0,&
                 & vkl0, nstunocc0, unitout
  use m_genfilname
! !INPUT/OUTPUT PARAMETERS:
! IN:
! integer :: iq ! q-point index 
! OUT:
! integer :: iocc0     ! Highest (partially) occupied state over all k-points
! integer :: iunocc0   ! Lowest (partially) unoccupied state over all k-points
! integer :: io0(nkpt) ! Highest (partially) occupied state for each k-point
! integer :: iu0(nkpt) ! Lowest (partially) unoccupied state for each k-point
! integer :: iocc, iunocc, io(nkpt), iu(nkpt) ! Same as above, but for k+q
! !NOTE: iocc=iocc0=max(iocc0,iocc) and iunocc0=iunocc=min(iunocc0,inunocc) was set! 
! INDIRECT OUT:
! integer :: modxs:nstocc0 = iocc0
! integer :: modxs:nstunocc0 = nstsv - iocc0
!
! !DESCRIPTION:
! For a given q-point index the routine inspects the occupancies and 
! second variational eigenvalues of the corresponding k+q and k set to determin
! the highest/lowest (at least partially) occupied/unoccupied 
! state index for all k+q and k points and the highest/lowest occupied/unoccupied
! state index over all k+q and k points.
! 
! !REVISION HISTORY:
!   Added description schema. And rudimentary description. (Aurich)
!   Added some more description. (Aurich)
!EOP
!BOC

!NOTE: iunocc is always set to iunocc0 (lowest unoccupied state over all k points)

  implicit none

  ! Arguments
  integer, intent(in) :: iq
  integer, intent(out) :: iocc0, iocc, iunocc0, iunocc
  integer, intent(out) :: io0(nkpt), io(nkpt), iu0(nkpt), iu(nkpt)

  ! Local variables
  real(8) :: thegap
  integer :: ik, ikq, i0, i
  logical :: t

  t = allocated(evalsv0)
  if( .not. t) allocate(evalsv0(nstsv, nkpt))

  k: do ik = 1, nkpt

    ! k+q-point set
    ikq = ik
    ! NOTE: ikmapikq has no iq=0 index
    if(iq .ne. 0) ikq = ikmapikq(ik, iq)

    ! Get occupancies and eigenvalues
    call getoccsv(vkl(1, ikq), occsv(1, ikq))
    call getevalsv(vkl(1, ikq), evalsv(1, ikq))

    ! Go from low to high through the sates at k+q, and
    ! check whether they are (partially) occupied, in the sense that
    ! they contribute to the electron density (epsocc defaults to 10^-8). 
    ! Save for each k+q point the index of the highest (partially) occupied state.
    do i = 1, nstsv
      if(occsv(i, ikq) .lt. input%groundstate%epsocc) exit
    end do
    io(ik) = i - 1

    ! Go from high to low through the sates at k+q, and
    ! check whether they are (partially) unoccupied and stop if the 
    ! state is fully occupied.
    ! Save for each k+q point the index of the lowest (partially) unoccupied state.
    do i = nstsv, 1, -1
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
    ! Lowest energy over all k+q and k
    evlmin = min(minval(evalsv(1, :)), minval(evalsv0(1, :)))
    ! Highest energy over all k+q and k
    evlmax = max(maxval(evalsv(nstsv, :)), maxval(evalsv0(nstsv, :)))

    ! Lower and higher cutoff valence energy
    ! Highest lowest energy over all k+q and k
    evlmincut = max(maxval(evalsv(1, :)), maxval(evalsv0(1, :)))
    ! Lowest highest energy over all k+q and k
    evlmaxcut = min(minval(evalsv(nstsv, :)), minval(evalsv0(nstsv, :)))

  else

    ! Lowest and highest valence energy
    ! Lowest energy over all k
    evlmin = minval(evalsv(1, :))
    ! Highest energy over all k
    evlmax = maxval(evalsv(nstsv, :))

    ! Lower and higher cutoff valence energy
    ! Highest lowest energy over all k
    evlmincut = maxval(evalsv(1, :))
    ! Lowest highest energy over all k
    evlmaxcut = minval(evalsv(nstsv, :))

  end if

  ! Find the highest (partially) occupied state over all k (k+q) 
  iocc0 = maxval(io0)
  iocc = maxval(io)

  ! Find the lowest (partially) unoccupied state over all k (k+q)
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

  ! Gap estimate 
  if(ksgap) then
    thegap = evllpu - evlhpo
  end if

  ! Assign nstocc0 and nstunocc0
  nstocc0 = iocc0
  nstunocc0 = nstsv - nstocc0
!!! iocc = iocc0 and iunocc = iunocc0 was set above !
  if((iocc0 .ge. iunocc) .or. (iocc .ge. iunocc0)) then
    write(unitout, '(a)') 'Info(findocclims): Partially occupied states present'
  end if
  if(ksgap) then
    write(unitout, '(a)') 'Info(findocclims): System has kohn-sham gap'
    write(unitout, '(a,E23.16)') 'Info(findocclims): Gap/H: ', thegap
    write(unitout, '(a,E23.16)') 'Info(findocclims): Gap/eV: ', thegap*h2ev
  else
    write(unitout, '(a)') 'Info(findocclims): No kohn-sham gap found'
  end if

  ! Debug output
  if(input%xs%dbglev .gt. 0) then
    write(*, '(a)') 'Debug(findocclims):'
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
