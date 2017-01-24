! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_screenlauncher
! !INTERFACE:
subroutine b_screenlauncher
! !USES:
  use modmpi
  use modinput, only: input
  use mod_misc, only: filext
  use mod_APW_LO, only: lolmax
  use mod_kpoint, only: nkpt, vkl, ikmap, ikmapnr
  use mod_qpoint, only: nqpt, vql
  use modxs, only: tscreen, xsgnt, nwdf, qpari,&
                   & qparf, unitout, nqmt, vqlmt, qvkloff,&
                   & gqdirname, eps0dirname, scrdirname, timingdirname,&
                   & ikmapikq, nkpt0, vkl0, usefilext0, filext0, iqmt0, iqmt1
  use mod_xsgrids
  use m_genfilname
  use m_filedel
  use m_writegqpts
  use m_xsgauntgen
  use m_findgntn0
  use mod_Gkvector, only: gkmax
! !DESCRIPTION:
! 
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich)
!EOP
!BOC
  implicit none

  integer(4) :: iqmt, iq, igq, qi, qf, ik
  logical :: firstisgamma
  logical :: fcoup, fti
  integer(4), parameter :: iqmtgamma = 1
  integer(4) :: ispin
  real(8), parameter :: epslat=1.d-6
  real(8) :: qgridoff(3)
  character(256) :: filex, syscommand
  character(*), parameter :: thisnam = 'b_screenlauncher'

  !write(*,*) "b_screenlauncher here at rank ", rank

  ! Initialise universal variables
  call init0
  ! Setting up k and G+k variables
  call init1
  ! Save k and G+k grid variables to 
  ! modxs (vkl0, ngk0, ...)
  call xssave0
  ! q-point and qmt-point setup
  !   Init 2 sets up (task 430):
  !   * A list of momentum transfer vectors form the q-point list (modxs::vqmtl)
  !   * The reduced unshifted q-grid (mod_qpoint::vql etc)
  !   * Offset of the k+q grid derived from k offset an q points (modxs::qvkloff)
  !     which is just equal to the k-offset, since it is the unshifted q grid
  !   * mapping between iknr,q and ik' grids (modxs::ikmapikq)
  !   * G+q quantities for the reduced unshifted q-grid (modxs)
  !   * The square root of the Coulomb potential for the reduced q points
  !   * Reads STATE.OUT
  !   * Generates radial functions (mod_APW_LO)
  call init2

  ! First Q-point in the list needs to be the Gamma point
  if(NORM2(vqlmt(1:3,1)) > epslat) then 
    firstisgamma = .false.
  else
    firstisgamma = .true.
  end if
  if(.not. firstisgamma) then 
    if(rank == 0) then 
      write(*,*) "First Q-point needs to be the gamma point."
    end if
    call terminate
  end if

  ! Write out q-points
  if(rank == 0) then
    call writeqpts
  end if

  !------------------------------------------------------------!
  ! Setting up folders for screening calculations              !
  !   Note: This is system dependend                           !
  !------------------------------------------------------------!
  ! Making folder for GQPOINTS info files
  gqdirname = 'GQPOINTS'
  if(rank == 0) then 
    syscommand = '[[ ! -e '//trim(adjustl(gqdirname))//' ]] && mkdir '//trim(adjustl(gqdirname))
    call system(trim(adjustl(syscommand)))
  end if
  ! Making folder for the binary screening EPS0 files
  eps0dirname = 'EPS0'
  if(rank == 0) then 
    syscommand = '[[ ! -e '//trim(adjustl(eps0dirname))//' ]] && mkdir '//trim(adjustl(eps0dirname))
    call system(trim(adjustl(syscommand)))
  end if
  ! Making folder for the ascii screening SCREEN files
  scrdirname = 'SCREEN'
  if(rank == 0) then
    syscommand = '[[ ! -e '//trim(adjustl(scrdirname))//' ]] && mkdir '//trim(adjustl(scrdirname))
    call system(trim(adjustl(syscommand)))
  end if
  ! Making folder for timing related info output
  timingdirname = 'TIMINGS'
  if(rank == 0) then
    syscommand = '[[ ! -e '//trim(adjustl(timingdirname))//' ]] && mkdir '//trim(adjustl(timingdirname))
    call system(trim(adjustl(syscommand)))
  end if
  call barrier
  !------------------------------------------------------------!


  !------------------------------------------------------------!
  ! Preparation for plane wave matrix elements                 !
  !------------------------------------------------------------!
  ! Generate gaunt coefficients, and store them in modxs:xsgnt
  call xsgauntgen(max(input%groundstate%lmaxapw, lolmax),&
    & input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))
  ! Find indices for non-zero gaunt coefficients, and store
  ! relevant maps in the module m_findgntn0
  call findgntn0(max(input%xs%lmaxapwwf, lolmax),&
    & max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)
  !------------------------------------------------------------!

  ! Only one frequency if BSE is used (w=0)
  nwdf = 1

  !------------------------------------------------------------!
  ! Make Screening for qmt=0 and calculate it for all q-points !
  ! q=k'-k (or the reduced q point set)                        !
  !------------------------------------------------------------!
  if(rank == 0) then 
    call printline(unitout, "+")
    write(unitout, '(a, i8)') 'Info(' // thisnam // '):&
      & Calculating screening for unshifted q-grid.'
    call printline(unitout, "+")
    write(unitout, *)
  end if

  ! Set *_SCR_QMT001.OUT as bra state file
  usefilext0 = .true.
  iqmt0 = iqmtgamma
  call genfilname(iqmt=iqmt0, scrtype='', setfilext=.true.)
  filext0 = filext
  !write(*,*) "filext0 =", trim(filext0)

  ! Set *_SCR_QMT001.OUT as ket state file
  iqmt1 = iqmtgamma
  call genfilname(iqmt=iqmt1, scrtype='', setfilext=.true.)
  !write(*,*) "filext=", trim(filext)

  ! Set type of band combinations: ({v,x},{x,c})- and ({x,c},{v,x})-combiantions
  input%xs%emattype = 1

  ! Use q point parallelization instead of w
  call genparidxran('q', nqpt)

  ! Loop over q-points 
  do iq = qpari, qparf

    ! Write q-point number to fileext, filext = "_SCR_QMT001_QXYZ.OUT"
    call genfilname(scrtype='', iqmt=iqmtgamma, iq=iq, fileext=filex)
    !write(*,*) "filex=", trim(filex)
    ! Write out G+q vectors to file "GQPOINTS_SCR_QMT001_QXYZ.OUT"
    call writegqpts(iq, filex, dirname=gqdirname)

    ! Generate screening for the given q-point
    call dfq(iq)

    write(unitout, '(a, i8)') 'Info(' // thisnam // '): Kohn Sham&
      & response function finished for q - point:', iq
    call printline(unitout, "-")

  end do

  if(rank == 0) then 
    call printline(unitout, "+")
    write(unitout, '(a, i8)') 'Info(' // thisnam // '):&
      & Screening for unshifted q-grid finished'
    call printline(unitout, "+")
    write(unitout, *)
  end if

  ! Synchronize
  call barrier
  !------------------------------------------------------------!

  !------------------------------------------------------------!
  ! In case that the coupling terms in the BSE matrix are      !
  ! taken into account, one potentially needs the screening    !
  ! also on shifted q-grids depending on vkloff, qmt and the   !
  ! coupling type.                                             !
  !------------------------------------------------------------!
  fcoup = input%xs%bse%coupling
  fti = input%xs%bse%ti
  ispin = 1

  if(fcoup) then 

    if(rank == 0) then
      call printline(unitout, "+")
      write(unitout, '(a)') 'Info(' // thisnam // '):&
        & Calculating screening on a shifted q-grid.'
      call printline(unitout, "+")
      write(unitout, *)
    end if

    ! Consider qmt vectors in qmt-list appart from gamma (already computed)
    do iqmt = 1, nqmt

      if(.not. fti .and. iqmt == 1) then 
        if(rank == 0) then 
          write(unitout, '(a, i3)') 'Info(' // thisnam // '):&
           & For vqmtl=0 no unshifted q grid is needed.'
          call printline(unitout, "-")
        end if
        cycle
      end if

      ! Generate q-qmt / -(q+qmt) grid offset
      call xsgrids_init(vqlmt(1:3, iqmt), gkmax)

      ! Reset k-grid variables to the unshifted (apart from xs%vkloff) k grid
      ! (because they get changed in dfq)
      call init1offs(k_kqmtp%kset%vkloff)
      ! Save k and G+k grid variables to 
      ! modxs (vkl0, ngk0, ...)
      call xssave0

      !------------------------------------------------------------!
      ! Generate q points and G+q quantities for a shifted q-grid. !
      ! The shift is determined by the k-grid offset and the       ! 
      ! considered momentum transver vector qmt.                   !
      !------------------------------------------------------------!
      if(rank == 0) then 
        write(unitout, '(a, i3)') 'Info(' // thisnam // '):&
         & Considering momentum transfer vector iqmt=', iqmt
        if(.not. fti) then 
          write(unitout, '(a)') 'Info(' // thisnam // '):&
            & Using q-qmt grid.'
        else
          write(unitout, '(a)') 'Info(' // thisnam // '):&
            & Using -(q+qmt) grid.'
        end if
        call printline(unitout, "-")
      end if

      if(.not. fti) then 
        ! q-qmt grid
        qgridoff = q_qmtm%qset%vkloff 
      else
        ! -(q+qmt) grid
        qgridoff = q_mqmtp%qset%vkloff 
      end if

      ! Free the xsgrids (only the offset was needed)
      call xsgrids_finalize()

      if(all(qgridoff == [0.0d0,0.0d0,0.0d0])) then 

        if(rank == 0) then 
          write(unitout, '(a, i3)') 'Info(' // thisnam // '):&
            & Shifted q-grid is identical to unshifted q-grid,&
            & skipping calculation.'
          call printline(unitout, "-")
        end if

        cycle

      end if

      ! Calculate the q vectors (mod_qpoint), G+q vectors with corresponding 
      ! structure factors and spherical harmonics and the square root
      ! of the Coulomb potential v^1/2(G,q) (modxs). Also creates ikmapikq
      ! that links (ik,iq) to ikp
      call init2offs(qgridoff, input%xs%reduceq)

      ! Set the file name extension to _SCR_QMTXYZ_mqmt.OUT / _SCR_QMTXYZ_m.OUT
      if(.not. fti) then 
        ! Set *_SCR_QMT001.OUT as bra state file
        usefilext0 = .true.
        iqmt0 = iqmtgamma
        call genfilname(iqmt=iqmt0, scrtype='', setfilext=.true.)
        filext0 = filext
        !write(*,*) "filext0 =", trim(filext0)

        ! Set *_SCR_QMTXYZ_mqmt.OUT as ket state file
        iqmt1 = iqmt 
        call genfilname(iqmt=iqmt1, scrtype='', auxtype='mqmt', setfilext=.true.)
        !write(*,*) "filext=", trim(filext)
      else
        ! Set *_SCR_QMT001.OUT as bra state file
        usefilext0 = .true.
        iqmt0 = iqmtgamma
        call genfilname(iqmt=iqmt0, scrtype='', setfilext=.true.)
        filext0 = filext
        !write(*,*) "filext0 =", trim(filext0)

        ! Set *_SCR_QMTXYZ_m.OUT as ket state file
        iqmt1 = iqmt 
        call genfilname(iqmt=iqmt1, scrtype='', auxtype='m', setfilext=.true.)
        !write(*,*) "filext=", trim(filext)
      end if

      ! Set type of band combinations: ({v,x},{x,c})- and ({x,c},{v,x})-combiantions
      input%xs%emattype = 1

      ! Use q point parallelization instead of w
      call genparidxran('q', nqpt)

      ! Loop over q-points 
      do iq = qpari, qparf

        ! Write q-point number to fileext, filext = _SCR_QMTXYZ_QXYZ_mqmt.OUT / _SCR_QMTXYZ_QXYZ_m.OUT
        if(.not. fti) then 
          call genfilname(scrtype='', iqmt=iqmt, auxtype='mqmt', iq=iq, fileext=filex)
        else
          call genfilname(scrtype='', iqmt=iqmt, auxtype='m', iq=iq, fileext=filex)
        end if
        !write(*,*) "filex=", trim(filex)
        ! Write out G+q vectors to file GQPOINTS_SCR_QMTXYZ_QXYZ_mqmt.OUT / GQPOINTS_SCR_QMTXYZ_QXYZ_m.OUT
        call writegqpts(iq, filex, dirname=gqdirname)

        ! Generate screening for the given q-point
        call dfq(iq)

        write(unitout, '(a, i4)') 'Info(' // thisnam // '): Kohn Sham&
          & response function finished for q - point:', iq
        call printline(unitout, "-")
        write(unitout, *)

      end do

      ! Synchronize
      call barrier

    end do

    if(rank == 0) then 
      call printline(unitout, "+")
      write(unitout, '(a)') 'Info(' // thisnam // '):&
        & Screening for shifted q-grid finished'
      call printline(unitout, "+")
      write(unitout, *)
    end if

  end if

  ! Delete gaunt maps
  call findgntn0_clear

  if(rank == 0) then
    write(unitout, '(a)') "Info(b_screenlauncher): Screening finished"
  end if
end subroutine b_screenlauncher
!EOC
