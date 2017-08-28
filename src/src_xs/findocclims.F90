! Copyright(C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findocclims
! !INTERFACE:
subroutine findocclims(iq, ikiq2ikp, iocc_common, iunocc_common, io0, io, iu0, iu)
! !USES:
  use mod_constants, only: h2ev
  use mod_kpoint, only: nkpt, vkl
  use mod_eigenvalue_occupancy, only: nstsv, occsv, evalsv, occmax,&
                                    & efermi
  use modinput, only: input
  use modxs, only: evalsv0, occsv0, evlmin,&
                 & evlmax, evlmincut, evlmaxcut, evlhpo,&
                 & evllpu, ksgap, ksgapval, qgap, nstocc0,&
                 & vkl0, nstunocc0, unitout
  use m_genfilname
! !INPUT/OUTPUT PARAMETERS:
! IN:
! integer :: iq ! q-point index 
! integer :: ikiq2ikp(nkpt) ! Index mapping from (ik,iq)-->ikp 
! OUT:
! integer :: iocc_common   ! Highest (partially) occupied state over all k-points and k+q-points
! integer :: iunocc_common ! Lowest (partially) unoccupied state over all k-points and k+q-ponts
! integer :: io0(nkpt) ! Highest (partially) occupied state for each k-point
! integer :: iu0(nkpt) ! Lowest (partially) unoccupied state for each k-point
! integer :: io(nkpt), iu(nkpt) ! Same as above, but for k+q
! INDIRECT OUT:
! integer :: modxs:nstocc0 = iocc0
! integer :: modxs:nstunocc0 = nstsv - iocc0
!
! !DESCRIPTION:
! For a given q-point index the routine inspects the occupancies and 
! second variational eigenvalues of the corresponding $k+q$ and $k$ set to determin
! the highest/lowest (at least partially) occupied/unoccupied 
! state index for all $k+q$ and $k$ points and the highest/lowest occupied/unoccupied
! state index over all $k+q$ and $k$ points.
! For the k set the files {\tt EVECSV\_QMT000.OUT} and {\tt OCCSV\_QMT000.OUT} is
! used, while for the $k+q$ set the files {\tt EVECSV\_QMTXXX.OUT} and 
! {\tt OCCSV\_QMTXXX.OUT} is used, where {\tt XXX} is the XXX'th element of the
! q-point list specified in the input file. For the input of $\text{iq}=0$ 
! the $k$ and $k+q$ sets are identical and the {\tt \_QMT000} files are used.
! 
! !REVISION HISTORY:
!   Added description schema. And rudimentary description. (Aurich)
!   Added some more description. (Aurich)
!   Added calculation of minimal (indirect) gap. (Aurich)
!   Added calculation of minimal gap for iq. (Aurich)
!   Change input output. (Aurich)
!EOP
!BOC

  implicit none

  ! Arguments
  integer, intent(in) :: iq
  integer, intent(in) :: ikiq2ikp(nkpt)
  integer, intent(out) :: iocc_common, iunocc_common
  integer, intent(out) :: io0(nkpt), io(nkpt), iu0(nkpt), iu(nkpt)

  ! Local variables
  real(8) :: t1, t0
  integer :: iocc0, iocc, iunocc0, iunocc
  integer :: ik, ikq, i0, i
  logical :: t

  !write(*,*)
  !write(*,*) "findocclims here"

  call timesec(t0)

  t = allocated(evalsv0)
  if( .not. t) allocate(evalsv0(nstsv, nkpt))
  ! Note: occsv0 is allocated during init1
  !   Only evalsv0 is not properly allocated/deallocated ?

  ! This routine is exclusively called withing the xs part at the moment.
  ! nkpt may reference the reduced or the non-reduced k set depending
  ! on input%xs%reducek. It defaults to false, so nkpt = nkptnr.
  ! Symmetry is not consistently used in XS.
  k: do ik = 1, nkpt

    ! k+q-point set
    ikq = ik
    if(iq .ne. 0) ikq = ikiq2ikp(ik)

    !write(*,*) "ik, ik'", ik, ikq
    !write(*,*) "findocclims: ik,ikq", ik, ikq
    !write(*,*) "findocclims: vkl0(:,ikq)", vkl0(:,ik)
    !write(*,*) "findocclims: vkl(:,ikq)", vkl(:,ikq)

    ! Get occupancies and eigenvalues form EVECSV and EVALSV 
    ! files that have file extensions that was set before calling
    ! findocclims.
    call getoccsv(vkl(1:3, ikq), occsv(1:nstsv, ikq))
    call getevalsv(vkl(1:3, ikq), evalsv(1:nstsv, ikq))

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
      call getoccsv0(vkl0(1:3, ik), occsv0(1:nstsv, ik))
      call getevalsv0(vkl0(1:3, ik), evalsv0(1:nstsv, ik))
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

  ! In case of no momentum transfer k and k+q quantities are the same.
  ! We do the copying here so we do not need to have different code for
  ! q=0 and q/=0 elsewhere.
  if(iq == 0) then 
    occsv0 = occsv
    evalsv0 = evalsv
  end if

  if(iq .ne. 0) then

    ! Lowest and highest valence energy 
    ! Lowest energy over all k+q and k
    evlmin = min(minval(evalsv(1, :)), minval(evalsv0(1, :)))
    !write(*,*) "findocclims: evlmin=", evlmin
    ! Highest energy over all k+q and k
    evlmax = max(maxval(evalsv(nstsv, :)), maxval(evalsv0(nstsv, :)))
    !write(*,*) "findocclims: evlmax=", evlmax

    ! Lower and higher cutoff valence energy
    ! Highest lowest energy over all k+q and k
    evlmincut = max(maxval(evalsv(1, :)), maxval(evalsv0(1, :)))
    !write(*,*) "findocclims: evlmincut=", evlmincut
    ! Lowest highest energy over all k+q and k
    evlmaxcut = min(minval(evalsv(nstsv, :)), minval(evalsv0(nstsv, :)))
    !write(*,*) "findocclims: evlmaxcut=", evlmaxcut

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
  !write(*,*) "findocclims: iocc0=", iocc0
  !write(*,*) "findocclims: iocc=", iocc

  ! Find the lowest (partially) unoccupied state over all k (k+q)
  iunocc0 = minval(iu0)
  iunocc = minval(iu)
  !write(*,*) "findocclims: iunocc0=", iunocc0
  !write(*,*) "findocclims: iunocc=", iunocc

  ! Calculate minimal q-gap (only reasonalble if system has a gap)
  if(iq /= 0) then 
    qgap = minval(evalsv(iunocc, ikiq2ikp(:)) - evalsv0(iocc0, :))
    !write(*,*) "findocclims: qgap=", qgap
  else
    qgap = minval(evalsv(iunocc, :) - evalsv0(iocc0, :))
  end if

  ! The maximum/minimum value is used since a shifted (k+q)-mesh which is not
  ! commensurate can cause partially occupied states that are absent for the
  ! k-mesh
  iocc_common = max(iocc0, iocc)
  !write(*,*) "findocclims: iocc_common =", iocc_common
  iunocc_common = min(iunocc0, iunocc)
  !write(*,*) "findocclims: iunocc_common =", iunocc_common

  ! Determine if system has a gap in energy
  if(iq .ne. 0) then
    ! Highest (partially) occupied state energy
    evlhpo = max(maxval(evalsv(iocc_common, :)), maxval(evalsv0(iocc_common, :)))
    !write(*,*) "findocclims: evlhpo=", evlhpo
    ! Lowest (partially) unoccupied state energy
    evllpu = min(minval(evalsv(iunocc_common, :)), minval(evalsv0(iunocc_common, :)))
    !write(*,*) "findocclims: evllpu=", evllpu
  else
    ! Highest (partially) occupied state energy
    evlhpo = maxval(evalsv(iocc_common, :))
    ! Lowest (partially) unoccupied state energy
    evllpu = minval(evalsv(iunocc_common, :))
  end if

  ! Determine if system has a gap in energy
  ksgap = evlhpo .lt. efermi
  !write(*,*) "findocclims: efermi=", efermi
  !write(*,*) "findocclims: ksgap=", ksgap

  ! Gap estimate 
  if(ksgap) then
    ksgapval = evllpu - evlhpo
  else
    ksgapval = 0.0d0
  end if

  ! Assign nstocc0 and nstunocc0
  nstocc0 = iocc_common
  nstunocc0 = nstsv - nstocc0

  if(iocc_common .ge. iunocc_common) then
    write(unitout, '(a)') 'Info(findocclims): Partially occupied states present'
  end if
  if(ksgap) then
    write(unitout, '(a)') 'Info(findocclims): System has Kohn-Sham gap'
    write(unitout, '(a,E23.16)') '  Gap/H: ', ksgapval
    write(unitout, '(a,E23.16)') '  Gap/eV: ', ksgapval*h2ev
    write(unitout, '(a)') 'Info(findocclims): Minimal gap for momentum tranfer:'
    write(unitout, '(a,I8)') '  iq: ', iq
    write(unitout, '(a,E23.16)') '  Gap(q)/H: ', qgap
    write(unitout, '(a,E23.16)') '  Gap(q)/eV: ', qgap*h2ev
  else
    write(unitout, '(a)') 'Info(findocclims): No Kohn-Sham gap found'
  end if

  ! Debug output
  if(input%xs%dbglev .gt. 0) then
    write(*, '(a)') 'Debug(findocclims):'
    write(*, '(a)') ' iocc0, iocc, iunocc0, iunocc, iocc_common, iunocc_common below:'
    write(*, '(4i8)') iocc0, iocc, iunocc0, iunocc, iocc_common, iunocc_common
    write(*, '(a)') ' ik, io0, iu, diff, io, iu0, diff below:'
    do ik = 1, nkpt
       write(*, '(7i8)') ik, io0(ik), iu(ik), iu(ik) - io0(ik),&
         & io(ik), iu0(ik), iu0(ik) - io(ik)
    end do
    write(*,*)
  end if

  if( .not. t) deallocate(evalsv0)

  call timesec(t1)
  write(unitout, '(a, f12.6)') 'Info(findocclims): Time needed/s = ', t1-t0

end subroutine findocclims
!EOC
