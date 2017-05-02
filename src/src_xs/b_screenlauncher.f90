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
  use mod_APW_LO, only: lolmax
  use mod_kpoint, only: nkpt
  use mod_qpoint, only: nqpt
  use modxs, only: xsgnt, nwdf, qpari,&
                   & qparf, unitout, nqmt, vqlmt, qvkloff,&
                   & gqdirname, eps0dirname, scrdirname, timingdirname,&
                   & ikmapikq, nkpt0, vkl0, usefilext0, filext0, filexteps, iqmt0, iqmt1
  use mod_xsgrids
  use m_genfilname
  use m_filedel
  use m_writegqpts
  use m_xsgauntgen
  use m_findgntn0
  use mod_Gkvector, only: gkmax
  use m_b_ematqk
! !DESCRIPTION:
! 
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich)
!EOP
!BOC
  implicit none

  integer(4) :: iqmt, iq, iqnr, iqmti, iqmtf
  logical :: firstisgamma
  logical :: fcoup, fti
  integer(4), parameter :: iqmtgamma = 1
  integer(4) :: ispin
  real(8), parameter :: epslat=1.d-8
  real(8) :: qgridoff(3), qgridoffgamma(3)
  character(256) :: filex, syscommand
  character(*), parameter :: thisname = 'b_screenlauncher'

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
      write(*,*) "Error(b_screenlauncher): First Q-point needs to be the gamma point."
    end if
    call terminate
  end if

  ! Write out q-points
  if(rank == 0) then
    call genfilname(iqmt=iqmtgamma, scrtype='', setfilext=.true.)
    call writeqpts
  end if

  !------------------------------------------------------------!
  ! Setting up folders for screening calculations              !
  !   Note: This is system dependend                           !
  !------------------------------------------------------------!
  ! Making folder for GQPOINTS info files
  gqdirname = 'GQPOINTS'
  if(rank == 0) then 
    syscommand = 'test ! -e '//trim(adjustl(gqdirname))//' && mkdir '//trim(adjustl(gqdirname))
    call system(trim(adjustl(syscommand)))
  end if
  ! Making folder for the binary screening EPS0 files
  eps0dirname = 'EPS0'
  if(rank == 0) then 
    syscommand = 'test ! -e '//trim(adjustl(eps0dirname))//' && mkdir '//trim(adjustl(eps0dirname))
    call system(trim(adjustl(syscommand)))
  end if
  ! Making folder for the ascii screening SCREEN files
  scrdirname = 'SCREEN'
  if(rank == 0) then
    syscommand = 'test ! -e '//trim(adjustl(scrdirname))//' && mkdir '//trim(adjustl(scrdirname))
    call system(trim(adjustl(syscommand)))
  end if
  ! Making folder for timing related info output
  timingdirname = 'TIMINGS'
  if(rank == 0) then
    syscommand = 'test ! -e '//trim(adjustl(timingdirname))//' && mkdir '//trim(adjustl(timingdirname))
    call system(trim(adjustl(syscommand)))
  end if
  call barrier(callername=trim(thisname))
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

  ! Read Fermi energy from file EFERMI
  ! Use EFERMI_SCR_QMT001.OUT (corresponding to the xs groundstate run for the unshifted k grid)
  call genfilname(iqmt=iqmtgamma, scrtype='', setfilext=.true.)
  call readfermi

  !------------------------------------------------------------!
  ! Make Screening for qmt=0 and calculate it for all q-points !
  ! q=k'-k (or the reduced q point set)                        !
  !------------------------------------------------------------!
  if(input%xs%screening%iqmtrange(1) == -1 .or. input%xs%screening%iqmtrange(1) == 1) then 

    if(rank == 0) then 
      call printline(unitout, "+")
      write(unitout, '(a, i8)') 'Info(' // thisname // '):&
        & Calculating screening for unshifted q-grid.'
      call printline(unitout, "+")
      write(unitout, *)
    end if

    ! Set *_SCR_QMT001.OUT as bra state file
    usefilext0 = .true.
    iqmt0 = iqmtgamma
    call genfilname(iqmt=iqmt0, scrtype='', fileext=filext0)
    !write(*,*) "filext0 =", trim(filext0)

    ! Set *_SCR_QMT001.OUT as ket state file
    iqmt1 = iqmtgamma
    call genfilname(iqmt=iqmt1, scrtype='', setfilext=.true.)
    !write(*,*) "filext=", trim(filext)

    ! Set *_QMT001.OUT as file extension for screening files
    call genfilname(iqmt=iqmt1, fileext=filexteps)
    !write(*,*) "filexteps=", trim(filexteps)

    ! Use <mk|e^{-i(G+q)r}|nk'> for q=k'-k in dfq
    emat_ccket = .false.
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

      write(unitout, '(a, i8)') 'Info(' // thisname // '): Kohn Sham&
        & response function finished for q - point:', iq
      call printline(unitout, "-")

    end do

    if(rank == 0) then 
      call printline(unitout, "+")
      write(unitout, '(a, i8)') 'Info(' // thisname // '):&
        & Screening for unshifted q-grid finished'
      call printline(unitout, "+")
      write(unitout, *)
    end if

    ! Synchronize
    call barrier(callername=trim(thisname))

  end if
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

    ! Save zero momentum transfer q grid offset
    call xsgrids_init(vqlmt(1:3, iqmtgamma), gkmax)
    if(fti) then 
      qgridoffgamma(1:3) = p_pqmtp%pset%vkloff
    else
      qgridoffgamma(1:3) = q_qmtm%qset%vkloff
    end if
    call xsgrids_finalize()

    if(rank == 0) then
      call printline(unitout, "+")
      write(unitout, '(a)') 'Info(' // thisname // '):&
        & Calculating screening on a shifted q-grid.'
      call printline(unitout, "+")
      write(unitout, *)
    end if

    iqmti = 1
    iqmtf = nqmt
    if(input%xs%screening%iqmtrange(1) /= -1) then 
      iqmti=input%xs%screening%iqmtrange(1)
      iqmtf=input%xs%screening%iqmtrange(2)
    end if

    ! Consider qmt vectors in qmt-list 
    do iqmt = iqmti, iqmtf

      if(.not. fti .and. iqmt == 1) then 
        if(rank == 0) then 
          write(unitout, '(a, i3)') 'Info(' // thisname // '):&
           & For vqmtl=0 no unshifted q grid is needed.'
          call printline(unitout, "-")
        end if
        cycle
      end if

      ! Generate q grid offset
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
        write(unitout, '(a, i3)') 'Info(' // thisname // '):&
         & Considering momentum transfer vector iqmt=', iqmt
        if(.not. fti) then 
          write(unitout, '(a)') 'Info(' // thisname // '):&
            & Using q = (kp-qmt)-k = q0-qmt grid.'
        else
          write(unitout, '(a)') 'Info(' // thisname // '):&
            & Using q = -(kp+qmt)-k grid.'
        end if
        call printline(unitout, "-")
      end if

      if(.not. fti) then 
        ! q=q0-qmt grid
        qgridoff = q_qmtm%qset%vkloff 
      else
        ! q=-(k'+qmt)-k grid
        qgridoff = p_pqmtp%pset%vkloff
      end if


      if(all(abs(qgridoff-qgridoffgamma) < epslat) .and. iqmt/=1) then 

        if(rank == 0) then 
          write(unitout, '("Info(",a,"):&
            & Shifted q-grid for iqmt=", i4, " is identical to q-grid,&
            & for iqmt=1. Skipping calculation.")') trim(thisname), iqmt
          call printline(unitout, "-")
        end if

        cycle

      end if

      ! Calculate the q vectors (mod_qpoint), G+q vectors with corresponding 
      ! structure factors and spherical harmonics and the square root
      ! of the Coulomb potential v^1/2(G,q) (modxs). 
      ! Also creates ikmapikq that links (ik,iq) to ikp ( k+q = k')
      ! and qvkloff which contains for each q point the offset of the k+q grid.
      call init2offs(qgridoff, input%xs%reduceq)
      if(fti) then
        ! Interested in q = -k'-k-qmt :
        ! <mk|e^{-i(G+q)r}|n k'> --> <mk|e^{-i(G+q)r}|n k+q> = <mk|e^{-i(G+q)r}|n -(k+qmt)>
        ! Exploiting time reversal symmetry:
        ! Using <mk|e^{-i(G+q)r}|(n k+qmt)^*> instead of <mk|e^{-i(G+q)r}|n -(k+qmt)>
        emat_ccket = .true.
        ! Then qvkloff needs to contain offset of k+qmt grid
        do iq=1,nqpt
          qvkloff(1:3,iq) = k_kqmtp%kqmtset%vkloff(1:3)
        end do
        ! and ikmapikq needs to link k,q -> k'+qmt instead of k,q -> -(k'+qmt)
        ! in order to use dfq with ematqk2 (i.e. <mk|e^{-i(G+q)r}|(n k+qmt)^*>)
        do iq=1,nqpt
          iqnr=p_pqmtp%pset%ikp2ik(iq)
          ikmapikq(1:nkpt,iq) = p_pqmtp%ikip2ikp_nr(1:nkpt,iqnr)
        end do
      else
        emat_ccket = .false.
      end if

      ! Free the xsgrids 
      call xsgrids_finalize()

      ! Write out q-points
      if(rank == 0) then
        if(fti) then 
          call genfilname(iqmt=iqmt, scrtype='', auxtype='m', setfilext=.true.)
        else
          call genfilname(iqmt=iqmt, scrtype='', auxtype='mqmt', setfilext=.true.)
        end if
        call writeqpts
      end if

      if(.not. fti) then 

        ! Set *_SCR_QMT001.OUT as bra state file
        usefilext0 = .true.
        iqmt0 = iqmtgamma
        call genfilname(iqmt=iqmt0, scrtype='', fileext=filext0)
        !write(*,*) "filext0 =", trim(filext0)

        ! Set *_SCR_QMTXYZ_mqmt.OUT as ket state file
        iqmt1 = iqmt 
        call genfilname(iqmt=iqmt1, scrtype='', auxtype='mqmt', setfilext=.true.)
        !write(*,*) "filext=", trim(filext)

        ! Set *_QMTXYZ_mqmt.OUT as filextension for the screening 
        call genfilname(iqmt=iqmt1, auxtype='mqmt', fileext=filexteps)
        !write(*,*) "filexteps=", trim(filexteps)

      else

        ! Set *_SCR_QMT001.OUT as bra state file
        usefilext0 = .true.
        iqmt0 = iqmtgamma
        call genfilname(iqmt=iqmt0, scrtype='', fileext=filext0)
        !write(*,*) "filext0 =", trim(filext0)

        ! Set *_SCR_QMTXYZ.OUT as ket state file
        iqmt1 = iqmt 
        call genfilname(iqmt=iqmt1, scrtype='', setfilext=.true.)
        !write(*,*) "filext=", trim(filext)

        ! Set *_QMTXYZ_m.OUT as filextension for the screening 
        call genfilname(iqmt=iqmt1, auxtype='m', fileext=filexteps)
        !write(*,*) "filexteps=", trim(filexteps)

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

        write(unitout, '(a, i4)') 'Info(' // thisname // '): Kohn Sham&
          & response function finished for q - point:', iq
        call printline(unitout, "-")
        write(unitout, *)

        if(mpiglobal%rank == 0) then
          write(6, '(a,"b_screenlauncher: Progess epsilon(q):", f10.3)', advance="no")&
            & achar( 13), 100.0d0*dble(iq-qpari+1)/dble(qparf-qpari+1)
          flush(6)
        end if

      end do

      if(mpiglobal%rank == 0) then
        write(6, *)
      end if

      ! Synchronize
      call barrier(callername=trim(thisname))

    end do

    if(rank == 0) then 
      call printline(unitout, "+")
      write(unitout, '(a)') 'Info(' // thisname // '):&
        & Screening for shifted q-grid finished'
      call printline(unitout, "+")
      write(unitout, *)
    end if

  end if

  ! Reset
  emat_ccket = .false.
  ! Delete gaunt maps
  call findgntn0_clear

  if(rank == 0) then
    write(unitout, '(a)') "Info(b_screenlauncher): Screening finished"
  end if
end subroutine b_screenlauncher
!EOC
