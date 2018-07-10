! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: screenlauncher
! !INTERFACE:
subroutine screenlauncher
! !USES:
  use modmpi
  use modinput, only: input
  use mod_APW_LO, only: lolmax
  use mod_kpoint, only: nkpt
  use mod_qpoint, only: nqpt
  use modxs, only: xsgnt, nwdf, qpari,&
                   & qparf, unitout, totalqlmt, qvkloff,&
                   & gqdirname, eps0dirname, scrdirname, timingdirname,&
                   & ikmapikq, nkpt0, vkl0, usefilext0, filext0, filexteps,&
                   & iqmt0, iqmt1, evalsv0, istocc0, istunocc0, isto0, isto,&
                   & istu0, istu, nst1, nst2, ngq
  use mod_xsgrids
  use m_genfilname
  use m_filedel
  use m_writegqpts
  use m_xsgauntgen
  use m_findgntn0
  use mod_Gkvector, only: gkmax
  use m_ematqk
! !DESCRIPTION:
!   This is a wrapper routine for the call of \texttt{dfq.f90} in
!   the screen task of the BSE calculation. It stats the calculation
!   of the Kohn-Sham dielectric matrix (in RPA) for the q-points
!   formed from the k-point differences $\vec{q} = \vec{k}'-\vec{k}$. 
!   If the BSE coupling blocks are included, the dielectric matrix is also calculated for
!   a shifted set of q-point which are formed form $q = -\vec{k}'-\vec{k}$.
! 
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich)
!EOP
!BOC
  implicit none

  integer(4) :: iq, iqnr
  logical :: firstisgamma
  logical :: fcoup
  integer(4), parameter :: iqmtgamma = 1
  integer(4) :: ispin
  real(8), parameter :: epslat=1.d-8
  real(8) :: pgridoff(3)
  character(256) :: filex, syscommand
  character(*), parameter :: thisname = 'screenlauncher'

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

  ! Clear modxs:ecalsv0 
  if(allocated(evalsv0)) deallocate(evalsv0)

  ! First Q-point in the list needs to be the Gamma point
  if(NORM2(totalqlmt(1:3,1)) > epslat) then 
    firstisgamma = .false.
  else
    firstisgamma = .true.
  end if
  if(.not. firstisgamma) then 
    if(rank == 0) then 
      write(*,*) "Error(screenlauncher): First Q-point needs to be the gamma point."
    end if
    call terminate
  end if

  ! Write out q-points
  if(rank == 0) then
    call genfilname(scrtype='', setfilext=.true.)
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

  ! Only one frequency is used in BSE (w=0)
  nwdf = 1

  ! Read Fermi energy from file EFERMI
  ! Use EFERMI_SCR.OUT (corresponding to the scr groundstate run
  ! for the reference k grid)
  call genfilname(scrtype='', setfilext=.true.)
  call readfermi

  !------------------------------------------------------------!
  ! Make Screening for for all q-points                        !
  ! q=k'-k (or the reduced q point set)                        !
  !------------------------------------------------------------!
  if(rank == 0) then 
    call printline(unitout, "+")
    write(unitout, '(a, i8)') 'Info(' // thisname // '):&
      & Calculating screening for unshifted q-grid.'
    write(unitout, '(a, i8)') 'Info(' // thisname // '):&
      & Number of q points:', nqpt
    call printline(unitout, "+")
    write(unitout, *)
  end if

  !! Set parameters for the plane-wave matrix element calculation,
  !! e.g. which files to use for bra and ket states.
  !
  ! Set *_SCR.OUT as bra state file
  usefilext0 = .true.
  iqmt0 = 1
  call genfilname(scrtype='', fileext=filext0)
  !
  ! Set *_SCR.OUT as ket state file
  iqmt1 = 1
  call genfilname(scrtype='', setfilext=.true.)
  !
  ! Use <mk|e^{-i(G+q)r}|nk'> for q=k'-k in dfq
  emat_ccket = .false.
  ! Set type of band combinations: ({v,x},{x,c})- and ({x,c},{v,x})-combiantions
  input%xs%emattype = 1

  ! Set *.OUT as file extension for screening files
  call genfilname(fileext=filexteps)
  ! write band combinations to INFOXS.OUT
  if (mpiglobal%rank == 0) then
    input%xs%emattype = 1
    do iq=1, nqpt
      call findocclims(iq, ikmapikq(:,iq), istocc0, istunocc0, isto0, isto, istu0, istu)
      call ematbdcmbs(input%xs%emattype)
      if (iq == 1) then
        write(unitout, '(a, 4i6)') 'Info(' // thisname // '):&
          & lowest (partially)  unoccupied state: ', istunocc0
        write(unitout, '(a, 4i6)') 'Info(' // thisname // '):&
          & highest (partially) occupied state  : ', istocc0
        write(unitout, '(a, i5)') 'Info(' // thisname // '):&
          & number of occupied states considered:',  nst1 
        write(unitout, '(a, i5)') 'Info(' // thisname // '):&
          & number of unoccupied states considered:',  nst2
      end if
      if (input%xs%BSE%outputlevelnumber == 1) then
        write(unitout, '(a, i6,a,i6)') 'Info(' // thisname // '):&
          & number of G + q vectors for q-point :', iq, ' : ', ngq(iq)
      end if
    end do
    call printline(unitout, "-")
  end if
  ! Use q point parallelization instead of frequency w points
  call genparidxran('q', nqpt)
  write(unitout,*)
  write(unitout, '(a, i4)') 'Info(' // thisname // '): Starting loop over q-points'
  write(unitout, *)
  ! Loop over q-points 
  do iq = qpari, qparf

    ! Write q-point number to fileext, filext = "_SCR_QXYZ.OUT"
    call genfilname(scrtype='', iq=iq, fileext=filex)

    ! Write out G+q vectors to file "GQPOINTS_SCR_QXYZ.OUT"
    call writegqpts(iq, filex, dirname=gqdirname)

    ! Generate screening for the given q-point
    call dfq(iq)

    if (input%xs%BSE%outputlevelnumber == 1) then
      write(unitout, '(a, i8)') 'Info(' // thisname // '): Kohn Sham&
        & response function finished for q - point:', iq
      call printline(unitout, "-")
    end if

#ifndef MPI
    if(mpiglobal%rank == 0) then
      write(6, '(a,"screenlauncher: Progess epsilon(q):", f10.3)', advance="no")&
        & achar( 13), 100.0d0*dble(iq-qpari+1)/dble(qparf-qpari+1)
      flush(6)
    end if
#endif

  end do

#ifndef MPI
  if(mpiglobal%rank == 0) then
    write(6, *)
  end if
#endif
  if(rank == 0) then 
    call printline(unitout, "+")
    write(unitout, '(a, i8)') 'Info(' // thisname // '):&
      & Screening for unshifted q-grid finished'
    call printline(unitout, "+")
    write(unitout, *)
  end if

  ! Synchronize
  call barrier(callername=trim(thisname))
  !------------------------------------------------------------!

  !------------------------------------------------------------!
  ! In case that the coupling terms in the BSE matrix are      !
  ! taken into account, one potentially needs the screening    !
  ! also on shifted q-grids (= p-grid) depending on vkloff.    !
  !------------------------------------------------------------!
  fcoup = input%xs%bse%coupling
  ispin = 1

  if(fcoup) then 

    ! Get offset of p=-k'-k grid
    call xsgrids_init(totalqlmt(1:3, 1), gkmax)
    pgridoff(1:3) = pqmt%pset%vkloff

    ! If no vkloff k-grid offset is set, then the p grid is 
    ! equal to the q grid.
    if(any(abs(pgridoff) > epslat)) then

      ! Reset k-grid variables to the unshifted (apart from xs%vkloff) k grid
      ! (because they get changed in dfq)
      call init1offs(k_kqmtp%kset%vkloff)
      ! Save k and G+k grid variables to 
      ! modxs (vkl0, ngk0, ...)
      call xssave0

      !------------------------------------------------------------!
      ! Generate p points and G+p quantities for a shifted q-grid. !
      ! The shift is determined by the k-grid offset.              !
      !------------------------------------------------------------!

      ! Calculate the q vectors (mod_qpoint), G+q vectors with corresponding 
      ! structure factors and spherical harmonics and the square root
      ! of the Coulomb potential v^1/2(G,q) (modxs). 
      ! Also creates ikmapikq that links (ik,iq) to ikp ( k+q = k')
      ! and qvkloff which contains for each q point the offset of the k+q grid.
      call init2offs(pgridoff, input%xs%reduceq)
      ! Interested in q = -k'-k :
      ! <mk|e^{-i(G+q)r}|n k'> --> <mk|e^{-i(G+q)r}|n k+q> = <mk|e^{-i(G+q)r}|n -k'>
      ! Exploiting time reversal symmetry:
      ! Using <mk|e^{-i(G+q)r}|(n k')^*> instead of <mk|e^{-i(G+q)r}|n -k'>
      emat_ccket = .true.
      ! Then qvkloff needs to contain offset of k grid
      do iq=1,nqpt
        qvkloff(1:3,iq) = k_kqmtp%kset%vkloff(1:3)
      end do
      ! and ikmapikq needs to link k,q -> k' instead of k,q -> -k'
      ! in order to use dfq with time reversal (i.e. <mk|e^{-i(G+q)r}|(n k')^*>)
      do iq=1,nqpt
        iqnr=pqmt%pset%ikp2ik(iq)
        ikmapikq(1:nkpt,iq) = pqmt%ikip2ikp_nr(1:nkpt,iqnr)
      end do

      ! Info output
      if(rank == 0) then
        call printline(unitout, "+")
        write(unitout, '(a)') 'Info(' // thisname // '):&
          & Calculating screening on a shifted q-grid (p-grid).'
        write(unitout, '(a)') 'Info(' // thisname // '):&
          & Using p = -kp-k grid.'
        write(unitout, '("q-grid offset derived form k offset:")')
        write(unitout, '("qvkloff = ", 3g18.10)') pgridoff(1:3)
        write(unitout, '("qvkloff/ngridk = ", 3g18.10)') pgridoff(1:3)/input%xs%ngridk
        write(unitout, '(a, i8)') 'Info(' // thisname // '):&
          & Number of q points:', nqpt
        call printline(unitout, "+")
        write(unitout, *)
      end if

      ! Free the xsgrids 
      call xsgrids_finalize()

      ! Write out q-points
      if(rank == 0) then
        call genfilname(scrtype='', auxtype='m', setfilext=.true.)
        call writeqpts
      end if

      !! Set which files to use in pwe construction
      !
      ! Set *_SCR.OUT as bra state file
      usefilext0 = .true.
      iqmt0 = 1
      call genfilname(scrtype='', fileext=filext0)
      !
      ! Set *_SCR.OUT as ket state file
      iqmt1 = 1
      call genfilname(scrtype='', setfilext=.true.)

      ! Set *_m.OUT as filextension for the screening 
      call genfilname(auxtype='m', fileext=filexteps)

      ! Set type of band combinations: ({v,x},{x,c})- and ({x,c},{v,x})-combiantions
      input%xs%emattype = 1

      ! Use q point parallelization instead of w
      call genparidxran('q', nqpt)

      write(unitout, '(a, i4)') 'Info(' // thisname // '): Starting loop over q-points'
      write(unitout, *)
      call printline(unitout, "-")
      ! Loop over q-points 
      do iq = qpari, qparf

        ! Write q-point number to fileext, filext = _SCR_QXYZ_m.OUT
        call genfilname(scrtype='', auxtype='m', iq=iq, fileext=filex)
        ! Write out G+q vectors to file GQPOINTS_SCR_QMTXYZ_QXYZ_m.OUT
        call writegqpts(iq, filex, dirname=gqdirname)

        ! Generate screening for the given q-point
        call dfq(iq)

        if (input%xs%BSE%outputlevelnumber == 1) then
          write(unitout, '(a, i4)') 'Info(' // thisname // '): Kohn Sham&
            & response function finished for q - point:', iq
          call printline(unitout, "-")
          write(unitout, *)
        end if
#ifndef MPI
        if((mpiglobal%rank == 0)) then
          write(6, '(a,"screenlauncher: Progess epsilon(q):", f10.3)', advance="no")&
            & achar( 13), 100.0d0*dble(iq-qpari+1)/dble(qparf-qpari+1)
          flush(6)
        end if
#endif MPI

      end do

#ifndef MPI
      if(mpiglobal%rank == 0) then
        write(6, *)
      end if
#endif MPI
      ! Synchronize
      call barrier(callername=trim(thisname))

      if(rank == 0) then 
        call printline(unitout, "+")
        write(unitout, '(a)') 'Info(' // thisname // '):&
          & Screening for shifted q-grid (p-grid) finished'
        call printline(unitout, "+")
        write(unitout, *)
      end if

    else

      if(rank == 0) then 
        call printline(unitout, "+")
        write(unitout, '(a)') 'Info(' // thisname // '):&
          & Shifted q-grid is idential to unshifted q-grid, skipping.'
        call printline(unitout, "+")
        write(unitout, *)
      end if

    end if

  end if

  ! Reset
  emat_ccket = .false.
  ! Delete gaunt maps
  call findgntn0_clear

  if(rank == 0) then
    write(unitout, '(a)') "Info(screenlauncher): Screening finished"
  end if

end subroutine screenlauncher
!EOC
