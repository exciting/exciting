! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: screen
! !INTERFACE:
subroutine screen
! !USES:
  use modmpi
!  use modxs, only: nwdf, unitout
  use m_genfilname
  use m_filedel

use mod_ematgrids
use m_ematqk
use mod_Gkvector, only: gkmax
use mod_eigensystem, only: nmatmax
use mod_eigenvalue_occupancy, only: nstfv
use mod_spin, only: nspnfv
use modxs
use m_b_ematqk
use mod_constants
use m_xsgauntgen
use m_findgntn0
use modinput
use mod_APW_LO
use m_writecmplxparts
use mod_symmetry, only: maxsymcrys
use mod_misc, only: task

use m_ematqk_orig
! !DESCRIPTION:
! This is a wrapper subroutine that backs up the number of 
! $ \omega $ points and calls the {\tt df} subroutine. Afterwards
! it restores them. It also sets the output file extensions for
! screening related quantities.
! 
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich)
!EOP
!BOC

  implicit none

  ! local variables
  integer :: nwdft


  !!+++++++TEST HACK ++++++++++++!!
  real(8) :: vqmtl(3)
  complex(8), allocatable :: eveck(:,:,:), eveckq(:,:,:), ematsepp(:,:,:), emat(:,:,:), ematsag(:,:,:)
  complex(8), allocatable :: evecsk(:,:,:,:)
  integer(4) :: i, ik, jk, ikkp, nkp, nkkp, ikq, iq, myngq, igq, ig, iqr, iqprev
  real(8) :: t0, t1, maxdiff, vq(3), vqr(3)
  integer(4) :: sc(maxsymcrys), ivgsc(3, maxsymcrys)
  integer(4), allocatable :: igqmap(:)
  integer(4) :: jsym, jsymi
  integer(4) :: nsc, ivgsym(3)
  type(bcbs) :: bc
  logical :: sag, sep

  write(*,*) "Screen here"
  call init0
  call init1
  call init2

  sep = .true.
  sag = .false.

  write(*,*) "Making grids"
  vqmtl = 0.0d0
  call ematgrids_init(vqmtl, gkmax)
  write(*,*) "Done"
  call ematgrids_write_grids()


  call genfilname(dotext='_SCR.OUT', setfilext=.true.)
  write(*,*) "Getting all evecs"
  allocate(evecsk(nmatmax, nstfv, nspnfv, kkqmtset%kset%nkpt))
  do ik = 1, kkqmtset%kset%nkpt
    write(*,*) "ik=", ik
    call getevecfv(kkqmtset%kset%vkl(1:3, ik), gkset%vgkl(1:3, 1:gkset%ngkmax, 1, ik),&
     & evecsk(:,:,:,ik))
  end do
  write(*,*) "Done"

  nkp = kkqmtset%kset%nkpt
  !nkp = 1

  nkkp = nkp*(nkp+1)/2

  if(sep) then 
    write(*,*) "Getting plane wave elements (Sepp modified)" 

    call timesec(t0)
    write(*,*) "emat init"
    call emat_initialize(lmaxapw_ = input%groundstate%lmaxmat, pre_rad_=.true., pre_fft_=.false.) 
    write(*,*) "emat init done"

    if(.false.) then
      write(*,*) "Precalculating ffts"
      do ik = 1, kkqmtset%kset%nkpt
        write(*,*) "FFT for ik=", ik
        call emat_precal_fft(ik, evecsk(:,:,1,ik))
      end do
    end if

    do i = 1, nkkp

      ikkp = qset%ikkp_qordered(i)
      write(*,*) "ikkp", ikkp
      call kkpmap(ikkp, nkp, ik, jk)
      write(*,*) "ik=", ik, " jk=", jk
      iq = qset%ikikp2iq_nr(ik, jk)
      myngq = gqset%ngknr(1,iq)
      ig = qset%ikikp2ig_nr(ik, jk)
      write(*,*) "iq=", iq, " ig=", ig
      ikq = qset%ikiq2ikp_nr(ik,iq)
      write(*,*) "ikq=", ikq

      if(allocated(emat)) deallocate(emat)
      allocate(emat(nstfv, nstfv, myngq))

     ! write(*,*) "ik=", ik
     ! write(*,*) "vkl=", kkqmtset%kset%vkl(1:3, ik)
     ! write(*,*) "iq=", iq
     ! write(*,*) "vql=", qset%qset%vkl(1:3, iq)
     ! write(*,*) "ikq=", ikq
     ! write(*,*) "vkql=", kkqmtset%kqmtset%vkl(1:3, ikq)
     ! write(*,*) "ngk=", gkset%ngk(1,ik)
     ! write(*,*) "ngkq=", gkqmtset%ngk(1,ikq)
     ! write(*,*) "ngq=", gqset%ngknr(1,iq)

      call emat_gen(ik, iq, evecsk(:,:,1,ik),evecsk(:,:,1,ikq),emat)

    end do
    call emat_finalize()
    call timesec(t1)
    write(*,*) "Finalize emat"
    write(*,*) "s=", t1-t0
    write(*,*) "done"
  end if

!  write(*,*) "Finalize ematgrids"
!  call ematgrids_finalize()
!  write(*,*) "Deallocating done"

  if(sag) then 
    write(*,*) "Getting plane wave elements (Sag)" 
    task = 440
    call init0
    call init1
    call init2
    call xssave0

    call timesec(t0)

    call xsgauntgen(max(input%groundstate%lmaxapw, lolmax),&
      & input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))
    call findgntn0(max(input%xs%lmaxapwwf, lolmax),&
      & max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)

    do iqr = 1, qset%qset%nkpt
      iq = qset%qset%ikp2ik(iqr)
      write(*,*) "Making radial integrals for iq=", iqr, "iqnr=", iq
      call putematrad(iqr,iq) 
    end do

    call ematqalloc
    iqprev=0

    do i = 1, nkkp

      ikkp = qset%ikkp_qordered(i)
      write(*,*) "ikkp", ikkp
      call kkpmap(ikkp, nkp, ik, jk)
      write(*,*) "ik=", ik, " jk=", jk
      iq = qset%ikikp2iq_nr(ik, jk)
      myngq = gqset%ngknr(1,iq)
      ig = qset%ikikp2ig_nr(ik, jk)
      write(*,*) "iq=", iq, " ig=", ig
      ikq = qset%ikiq2ikp_nr(ik,iq)
      write(*,*) "ikq=", ikq

      if(iqprev /= iq) then
        iqr= qset%qset%ik2ikp(iq)
        vq = qset%qset%vklnr(1:3,iq)
        vqr= qset%qset%vkl(1:3,iqr)

        ! Get radial integrals (previously calculated for reduced q set)
        write(*,*) "Fetching radial integral from disk for: iqnr=", iq, " iq=", iqr, "numgq=", myngq
        write(*,*) "vq=", vqr
        write(*,*) "vqnr=", vq
        call getematrad(iqr, iq)

        if(any(vqr /= vq)) then

          allocate(igqmap(myngq))
          call findsymeqiv(input%xs%bse%fbzq, vq, vqr, nsc, sc, ivgsc)
          write(*,*) "Number of symmetry operations for vq --> vqnr =", nsc
          ! Find the map that rotates the G-vectors
          call findgqmap(iq, iqr, nsc, sc, ivgsc, myngq, jsym, jsymi, ivgsym, igqmap)
          write(*,*) "map found."
          ! Rotate radial integrals calculated for the reduced q to get those for non-reduced q
          call rotematrad(myngq, igqmap)
          write(*,*) "integrals rotated"
          deallocate(igqmap)

        end if
      end if
      iqprev = iq

      if(allocated(ematsag)) deallocate(ematsag)
      allocate(ematsag(nstfv, nstfv, myngq))
      ematsag = zzero
      bc%n1 = nstfv
      bc%n2 = nstfv 
      bc%il1 = 1
      bc%il2 = 1 
      bc%iu1 = nstfv
      bc%iu2 = nstfv 
      call b_ematqk(iq, ik, ematsag, bc)
    end do
    call timesec(t1)
    write(*,*) "s=", t1-t0
    write(*,*) "done"

    !write(*,*) "Writing sagmeisters emat to text file"
    !call emat_write_plain('EMAT_SAG', ematsag)
    !write(*,*) "Done"
  end if

  if(sag .and. sep) then 
    write(*,*) "Maximal absolute difference:"
    write(*,*) maxval(abs(emat(:,:,:)-ematsag(:,:,:)))
  end if

  write(*,*) "Screen bye!"

  call terminate
  !!+++++++TEST HACK END ++++++++++++!!

  ! Back up number of energy points
  nwdft = nwdf

  ! Change file extension variable in mod_misc to '_SCR.OUT' for the
  ! following calculations
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  ! Call dielectric function with only one frequency point
  ! df is a wrapper for dfq(iq) which for 'screen' sets nwdf=1, i.e. static screening
  call df

  ! Alternative for checking only:
  nwdf = nwdft

  if(rank == 0) then
    write(unitout, '(a)') "Info(screen): Screening finished"
  end if
end subroutine screen
!EOC
