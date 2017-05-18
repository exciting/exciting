module m_setup_bse
      
  implicit none

  private

  public :: setup_bse_ti, setup_bse_full, setup_bse_block
  public :: setup_bse_ti_dist, setup_bse_block_dist

  contains

    !! Serial versions

    !BOP
    ! !ROUTINE: setup_bse_full
    ! !INTERFACE:
    subroutine setup_bse_full(ham, iqmt)
    ! !USES:
      use modmpi
      use modbse, only: hamsize
      use modxs, only: unitout
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   integer(4) :: iqmt  ! Index of momentum transfer Q
    !                       ! Note, currently only zero momentum transfer allowed
    ! In/Out:
    !   complex(8) :: ham(:,:) ! Full BSE matrix with coupling blocks
    ! 
    ! !DESCRIPTION:
    !   The routine sets up the full BSE matrix including
    !   resonant-resonant and resonant-antiresonant blocks.
    !
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !EOP
    !BOC

      ! I/O
      complex(8), intent(inout) :: ham(:, :)
      integer(4), intent(in) :: iqmt

      ! Timings
      real(8) :: ts0, ts1

      call timesec(ts0)
      write(unitout, '("Info(setup_bse_full): Setting up full hamiltonian")')

      ! RR
      call setup_bse_block(ham(1:hamsize,1:hamsize), iqmt, .false., .false.)

      ! RA
      call setup_bse_block(ham(1:hamsize,hamsize+1:hamsize*2), iqmt, .true., .false.)
      ! ham(1:hamsize,hamsize+1:hamsize*2) = zzero

      ! AR
      ! Note: AR part is the negative complex conjugate of RA even if qmt /= 0
      ham(hamsize+1:2*hamsize, 1:hamsize) = -conjg(ham(1:hamsize,hamsize+1:hamsize*2))
      ! ham(hamsize+1:2*hamsize, 1:hamsize) = zzero

      ! AA
      ! Note: AA part is the negative complex conjugate of RR ONLY if qmt /= 0
      ham(hamsize+1:2*hamsize, hamsize+1:2*hamsize) = -conjg(ham(1:hamsize,1:hamsize))

      call timesec(ts1)
      write(unitout, '(" Matrix build.")')
      write(unitout, '("Timing (in seconds)	   :", f12.3)') ts1 - ts0

    end subroutine setup_bse_full
    !EOC

    !BOP
    ! !ROUTINE: setup_bse_ti
    ! !INTERFACE:
    subroutine setup_bse_ti(iqmt, smat, cmat, cpmat)
    ! !USES:
      use modmpi
      use modinput, only: input
      use mod_constants, only: zzero, zone
      use modbse, only: hamsize
      use modxs, only: unitout
      use m_writecmplxparts
      use m_hesolver
      use m_sqrtzmat
      use m_invertzmat
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   integer(4) :: iqmt  ! Index of momentum transfer Q
    !
    ! Out:
    !   complex(8) :: smat(:,:)  ! Aux. EVP matrix S=(A-B)^{1/2}(A+B)(A-B)^{1/2}
    !   complex(8) :: cmat(:,:)  ! C = (A-B)^{1/2}
    !   complex(8) :: cpmat(:,:) ! C'= (A-B)^{-1/2}
    ! 
    ! !DESCRIPTION:
    !   The routine sets up the auxilliary EVP matrix.
    !
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !EOP
    !BOC

      ! I/O
      integer(4), intent(in) :: iqmt
      complex(8), allocatable, intent(inout), optional :: smat(:,:)
      complex(8), allocatable, intent(inout), optional :: cmat(:,:)
      complex(8), allocatable, intent(inout), optional :: cpmat(:,:)

      ! Local 
      character(*), parameter :: thisname = "setup_bse_ti"
      complex(8), allocatable :: rrmat(:, :)
      complex(8), allocatable :: ramat(:, :)
      complex(8), allocatable :: auxmat(:, :)


      integer(4) :: info, lwork
      integer(4), allocatable :: ipiv(:)
      complex(8), allocatable :: work(:)

      real(8) :: ts0, ts1, t1, t0
      integer(4) :: i, j

      real(8), allocatable :: evals(:)

      logical :: fwp
      logical :: fcheckpos, fsmat, fcmat, fcpmat

      fwp = input%xs%bse%writeparts
      fcheckpos = input%xs%bse%checkposdef

      fsmat = .false.
      fcmat = .false.
      fcpmat = .false.
      if(present(smat)) fsmat = .true.
      if(present(cmat)) fcmat = .true.
      if(present(cpmat)) fcpmat = .true.
      if(.not. (fsmat .or. fcmat .or. fcpmat)) then
        write(unitout, '("Error(",a,"):&
          & No matrix to build specified.")') trim(thisname)
        call terminate
      end if


      write(unitout, '("Info(setup_bse_ti): Setting up matrices for squared EVP")')
      call timesec(ts0)

      !===========================================================!
      ! Getting main and coupling blocks of BSE hamiltonian       !
      !===========================================================!
      ! Get RR part of BSE Hamiltonian

      write(unitout, '("Info(setup_bse_ti): Setting up RR Block of orignial BSE")')
      call timesec(t0)
      allocate(rrmat(hamsize,hamsize))
      call setup_bse_block(rrmat, iqmt, .false., .false.)
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0

      ! Get RA^{ti} part of BSE Hamiltonian
      write(unitout, '("Info(setup_bse_ti):&
        & Setting up RA^{ti} Block of orignial BSE")')
      call timesec(t0)
      allocate(ramat(hamsize,hamsize))
      call setup_bse_block(ramat, iqmt, .true., .true.)
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
      !===========================================================!

      !===========================================================!
      ! Make combination matrices                                 !
      !===========================================================!
      write(unitout, '("Info(",a,"):&
       & Setting up RR-RA matrix")') trim(thisname)
      ! RR - RA^{it}
      do j = 1, hamsize
        do i = 1, hamsize
          rrmat(i,j) = rrmat(i,j) - ramat(i,j)
        end do
      end do

      ! Construct S matrix - A+B only needed for S
      if(fsmat) then 
        write(unitout, '("Info(",a,"):&
         & Setting up RR+RA matrix")') trim(thisname)
        ! RR + RA^{it}
        do j = 1, hamsize
          do i = 1, hamsize
            ramat(i,j) = rrmat(i,j) + 2.0d0*ramat(i,j)
          end do
        end do
        if(fcheckpos) then 
          ! Check positive definitness of (A+B)
          write(unitout, '("Info(setup_bse_ti): Checking positve definitness of RR+RA")')
          call timesec(t0)
          allocate(auxmat(hamsize,hamsize))
          auxmat = ramat
          if(fwp) then
            call writecmplxparts("nd_apb_mat", remat=dble(auxmat), immat=aimag(auxmat))
          end if
          allocate(evals(hamsize))
          call hesolver(auxmat, evals)
          !write(*,*) "writing apm evals"
          if(fwp) then
            call writecmplxparts("nd_apb_evals", revec=evals, veclen=size(evals))
          end if
          if(any(evals < 0.0d0)) then 
            write(*,*) "Error(setup_bse_ti): A+B matrix is not positive definit"
            write(*,'(E10.3)') evals
            call terminate
          end if
          call timesec(t1)
          write(unitout, '("  RR+RA is positive definite")')
          write(unitout, '("  Time needed",f12.3,"s")') t1-t0
          deallocate(evals)
          deallocate(auxmat)
        else
          write(unitout, '("  RR+RA is assumed to be positive definite")')
        end if
      ! Do not construct S
      else
        ! B no longer needed if S is not to be constructed (only A-B)
        deallocate(ramat)
      end if

      ! Test writeout
      if(fwp) then 
        call writecmplxparts("nd_amb_mat", remat=dble(cpmat), immat=aimag(cpmat))
      end if

      !===========================================================!

      !===========================================================!
      ! Take the square root of (A-B)                             !
      ! Note: It is assumed to be positive definit.               !
      !===========================================================!
      write(unitout, '("Info(setup_bse_ti): Taking square root of RR-RA matrix")')
      call timesec(t0)
      ! rrmat -> (A-B)^1/2
      call sqrtzmat_hepd(rrmat)
      call timesec(t1)
      write(unitout, '("  RR-RA is positive definite")')
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
      if(fwp) then
        call writecmplxparts("nd_sqrtamb_mat", remat=dble(cpmat), immat=aimag(cpmat))
      end if
      !===========================================================!

      if(fsmat) then 
        !===========================================================!
        ! Construct S Matrix                                        !
        ! S = (A-B)^{1/2} (A+B) (A-B)^{1/2}                         !
        !===========================================================!

        write(unitout, '("Info(setup_bse_ti): Constructing S matrix")')
        call timesec(t0)

        !   C = (A-B)^{1/2} (A+B)
        allocate(auxmat(hamsize,hamsize))
        call zgemm('N','N', hamsize, hamsize, hamsize, zone, rrmat, hamsize,&
          & ramat, hamsize, zzero, auxmat, hamsize)
        deallocate(ramat)

        !   S = C (A-B)^{1/2}
        allocate(smat(hamsize, hamsize))
        call zgemm('N','N', hamsize, hamsize, hamsize, zone, auxmat, hamsize,&
          & rrmat, hamsize, zzero, smat, hamsize)
        deallocate(auxmat)

        !write(*,*) "printing s mat"
        if(fwp) then 
          call writecmplxparts('nd_s_mat', remat=dble(smat), immat=aimag(smat))
        end if

        call timesec(t1)
        write(unitout, '("  Time needed",f12.3,"s")') t1-t0
        !===========================================================!
      end if

      if(fcmat) then 
        !===========================================================!
        ! Make C Matrix                                             !
        ! Cmat = (A-B)^1/2                                          !
        ! Cmat is needed to build Oscillator strengths              !
        !===========================================================!
        write(unitout, '("Info(",a,"):&
          & Retruning C=(RR-RA)^1/2 matrix")') trim(thisname)
        allocate(cmat(hamsize,hamsize))
        do j = 1, hamsize
          do i = 1, hamsize
            cmat(i,j) = rrmat(i,j)
          end do
        end do
        !===========================================================!
      end if

      if(fcpmat) then 
        !===========================================================!
        ! Make C' Matrix                                            !
        ! Cpmat = (A-B)^{-1/2}                                      !
        ! Cpmat is needed additionally to build eigenvectors        !
        !===========================================================!
        write(unitout, '("Info(",a,"):&
          & Inverting C=(RR-RA)^1/2 matrix")') trim(thisname)
        call timesec(t0)

        call zinvert(rrmat)

        ! "Move" to output array (rrmat no longer needed)
        call move_alloc(rrmat, cpmat)

        call timesec(t1)
        write(unitout, '("  Time needed",f12.3,"s")') t1-t0
        !===========================================================!
      end if

      call timesec(ts1)
      write(unitout, '("Info(setup_bse_ti): Total time needed",f12.3,"s")') ts1-ts0

    end subroutine setup_bse_ti
    !EOC

    !BOP
    ! !ROUTINE: setup_bse_block
    ! !INTERFACE:
    subroutine setup_bse_block(ham, iqmt, fcoup, fti)
    ! !USES:
      use modmpi
      use modinput, only: input
      use mod_constants, only: zzero, zone
      use modxs, only: unitout
      use modbse, only: de, nou_bse_max, hamsize, kmap_bse_rg,&
                      & kousize, nkkp_bse, nk_bse, ofac,&
                      & scclifbasename, exclifbasename,&
                      & scclicfbasename, exclicfbasename,&
                      & scclictifbasename, infofbasename,&
                      & ensortidx
      use m_getunit
      use m_genfilname
      use m_putgetbsemat
      use m_writecmplxparts
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   integer(4) :: iqmt  ! Index of momentum transfer Q
    !   logical    :: fcoup ! If true, builds RA instead of RR block of BSE matrix 
    !   logical    :: fti   ! If true, uses time inverted anti-resonant basis
    ! In/Out:
    !   complex(8) :: ham(:,:) ! RR or RA block of BSE-Hamiltonian matrix
    ! 
    ! !DESCRIPTION:
    !   The routine sets up the resonant-resonant or resonant-antiresonant block of
    !   the BSE-Hamiltonian matrix. The routine reads {\tt EXCLI.OUT} and
    !   {\tt SCCLI.OUT} ({\tt EXCLIC.OUT} and {\tt SCCLIC.OUT}).
    !
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !EOP
    !BOC

      !! I/O
      complex(8), intent(inout) :: ham(:, :)
      integer(4), intent(in) :: iqmt
      logical, intent(in) :: fcoup, fti

      !! Local variables
      ! Indices
      integer(4) :: ikkp, ik, jk, iknr, jknr
      integer(4) :: inou, jnou
      integer(4) :: i1, i2, j1, j2
      ! Work arrays
      complex(8), dimension(nou_bse_max, nou_bse_max) :: excli_t, sccli_t

      ! Test write out vars 
      real(8), allocatable :: diag(:,:)
      complex(8), allocatable :: wint(:,:), vint(:,:), tmp(:,:)
      logical :: fwp

      ! Read in vars
      character(256) :: sfname, sinfofname
      character(256) :: efname, einfofname
      logical :: efcmpt, efid
      logical :: sfcmpt, sfid
      logical :: useexc, usescc

      ! Timings
      real(8) :: ts0, ts1
      

      call timesec(ts0)
      if(.not. fcoup) then 
        write(unitout, '("Info(setup_bse_block): Setting up RR part of hamiltonian")')
      else
        write(unitout, '("Info(setup_bse_block): Setting up RA part of hamiltonian")')
        if(fti) then 
          write(unitout, '("Info(setup_bse_block):&
            & Using time inverted anti-resonant basis")')
        end if
      end if

      usescc = .false.
      useexc = .false.
      select case(trim(input%xs%bse%bsetype))
        case('IP')
          usescc = .false.
          useexc = .false.
          if(mpiglobal%rank==0) then
            write(unitout, '("Info(setup_bse_block): Using IP")')
          end if
        case('singlet')
          usescc = .true.
          useexc = .true.
          if(mpiglobal%rank==0) then
            write(unitout, '("Info(setup_bse_block): Using singlet")')
          end if
        case('triplet')
          usescc = .true.
          useexc = .false.
          if(mpiglobal%rank==0) then
            write(unitout, '("Info(setup_bse_block): Using triplet")')
          end if
        case('RPA')
          usescc = .false.
          useexc = .true.
          if(mpiglobal%rank==0) then
            write(unitout, '("Info(setup_bse_block): Using RPA")')
          end if
      end select

      ! Write out W and V (testing feature)
      fwp = input%xs%bse%writeparts
      if(fwp) then
        allocate(diag(hamsize, 1))
        allocate(wint(hamsize, hamsize))
        allocate(vint(hamsize, hamsize))
        allocate(tmp(hamsize, hamsize))
      end if

      ! Zero ham
      ham = zzero

      ! Select form which files W and V are to be read
      if(.not. fcoup) then
        call genfilname(basename=scclifbasename, iqmt=iqmt, filnam=sfname)
        call genfilname(basename=exclifbasename, iqmt=iqmt, filnam=efname)
      else
        if(fti) then
          call genfilname(basename=scclictifbasename, iqmt=iqmt, filnam=sfname)
          call genfilname(basename=exclifbasename, iqmt=iqmt, filnam=efname)
        else
          call genfilname(basename=scclicfbasename, iqmt=iqmt, filnam=sfname)
          call genfilname(basename=exclicfbasename, iqmt=iqmt, filnam=efname)
        end if
      end if

      sinfofname = trim(infofbasename)//'_'//trim(sfname)
      einfofname = trim(infofbasename)//'_'//trim(efname)

      if(usescc) write(unitout, '("  Reading form info from ", a)')&
        & trim(sinfofname)
      if(useexc) write(unitout, '("  Reading form info from ", a)')&
        & trim(einfofname)

      ! Check saved Quantities for compatiblity
      sfcmpt=.false.
      sfid=.false.
      if(usescc) call b_getbseinfo(trim(sinfofname), iqmt,&
        & fcmpt=sfcmpt, fid=sfid)
      efcmpt=.false.
      efid=.false.
      if(useexc) call b_getbseinfo(trim(einfofname), iqmt,&
        & fcmpt=efcmpt, fid=efid)
      if(usescc .and. useexc) then 
        if(efcmpt .neqv. sfcmpt .or. efid .neqv. sfid) then 
          write(*, '("Error(setup_bse_block): Info files differ")')
          write(*, '("  efcmpt, efid, sfcmpt, sfcmpt")') efcmpt, efid, sfcmpt, sfcmpt
          call terminate
        end if
      end if

      if(usescc) write(unitout, '("  Reading form W from ", a)') trim(sfname)
      if(usescc) write(unitout, '("  compatible:",l," identical:",l)') sfcmpt, sfid
      if(useexc) write(unitout, '("  Reading form V from ", a)') trim(efname)
      if(usescc) write(unitout, '("  compatible:",l," identical:",l)') efcmpt, efid

      ! Set up kkp blocks of RR or RA Hamiltonian
      !! Note: If the Hamilton matrix 
      !! has the elements H^RR_{i,j}(qmt) and the indices enumerate the
      !! states according to
      !! i = {u1o1k1, u2o1k1, ..., uM(k1)o1k1,
      !!      uM(k1)o2k1, ..., uM(k1)oN(k1)k1, u1o1k2,
      !!      ..., uM(kO)oN(kO)kO} -> {1,...,\Prod_{i=1}^O M(ki)*N(ki)}
      !! then because of H^RR_{j,i} = H^RR*_{i,j} (or H^RA_{j,i} = H^RA_{i,j})
      !! only jk = ik,..,kO is needed.
      do ikkp = 1, nkkp_bse

        ! Get index ik jk form combination index ikkp
        call kkpmap(ikkp, nk_bse, ik, jk)

        ! Get global index iknr
        iknr = kmap_bse_rg(ik)
        inou = kousize(iknr)

        ! Get global index jknr
        jknr = kmap_bse_rg(jk)
        jnou = kousize(jknr)

        ! Read corresponding ikkp blocks of W and V from file
        if(usescc) then 
          call b_getbsemat(trim(sfname), iqmt, ikkp,&
            & sccli_t(1:inou,1:jnou), check=.false., fcmpt=sfcmpt, fid=sfid)
        end if

        if(useexc) then 
          ! Read RR/RA part of exchange interaction v_{iuioik,jujojk}(qmt)
          call b_getbsemat(trim(efname), iqmt, ikkp,&
            & excli_t(1:inou,1:jnou), check=.false., fcmpt=efcmpt, fid=efid)
        end if

        if(fwp) then 
          if(ikkp == 1 .and. fcoup .eqv. .false.) then 
            if(usescc) call writecmplxparts('ikkp1_read_Wrr',dble(sccli_t(1:inou,1:jnou)), aimag(sccli_t(1:inou,1:jnou)))
            if(useexc) call writecmplxparts('ikkp1_read_Vrr',dble(excli_t(1:inou,1:jnou)), aimag(excli_t(1:inou,1:jnou)))
          end if

          if(ikkp == 1 .and. fcoup .eqv. .true. .and. fti .eqv. .false.) then 
            if(usescc) call writecmplxparts('ikkp1_read_Wra',dble(sccli_t(1:inou,1:jnou)), aimag(sccli_t(1:inou,1:jnou)))
            if(useexc) call writecmplxparts('ikkp1_read_Vra',dble(excli_t(1:inou,1:jnou)), aimag(excli_t(1:inou,1:jnou)))
          end if

          if(ikkp == 1 .and. fcoup .eqv. .true. .and. fti .eqv. .true.) then 
            if(usescc) call writecmplxparts('ikkp1_read_Wra_ti',dble(sccli_t(1:inou,1:jnou)), aimag(sccli_t(1:inou,1:jnou)))
            if(useexc) call writecmplxparts('ikkp1_read_Vra_ti',dble(excli_t(1:inou,1:jnou)), aimag(excli_t(1:inou,1:jnou)))
          end if
        end if

        ! Position of ikkp block in global matrix
        i1 = sum(kousize(1:iknr-1)) + 1
        j1 = sum(kousize(1:jknr-1)) + 1
        i2 = i1+inou-1
        j2 = j1+jnou-1

        !! RR and RA part
        ! Add correlation term and optionally exchange term
        ! (2* v_{iu io ik, ju jo jk}(qmt)
        ! - W_{iu io ik, ju jo jk})(qmt)
        ! * sqrt(abs(f_{io ik} - f_{ju jk+qmt}))
        ! * sqrt(abs(f_{jo jk}-f_{ju jk+qmt})) or sqrt(abs(f_{jo jk}-f_{ju jk-qmt}))
        if(usescc .and. useexc) then 
          if(fwp) then 
            !write(*,*) "i1:i2,j1:j2", i1,i2,j1,j2
            call setint(ham(i1:i2,j1:j2),&
              & ofac(i1:i2), ofac(j1:j2),&
              & scc=sccli_t(1:inou,1:jnou), exc=excli_t(1:inou,1:jnou),&
              & w=wint(i1:i2,j1:j2), v=vint(i1:i2,j1:j2))
          else
            call setint(ham(i1:i2,j1:j2),&
              & ofac(i1:i2), ofac(j1:j2),&
              & scc=sccli_t(1:inou,1:jnou), exc=excli_t(1:inou,1:jnou))
          end if
        else if(usescc) then
          if(fwp) then 
            call setint(ham(i1:i2,j1:j2),&
              & ofac(i1:i2), ofac(j1:j2),&
              & scc=sccli_t(1:inou,1:jnou),&
              & w=wint(i1:i2,j1:j2))
          else
            call setint(ham(i1:i2,j1:j2),&
              & ofac(i1:i2), ofac(j1:j2),&
              & scc=sccli_t(1:inou,1:jnou))
          end if
        else if(useexc) then 
          if(fwp) then 
            call setint(ham(i1:i2,j1:j2),&
              & ofac(i1:i2), ofac(j1:j2),&
              & exc=excli_t(1:inou,1:jnou),&
              & v=vint(i1:i2,j1:j2))
          else
            call setint(ham(i1:i2,j1:j2),&
              & ofac(i1:i2), ofac(j1:j2),&
              & exc=excli_t(1:inou,1:jnou))
          end if
        end if

        !! RR only
        ! For blocks on the diagonal, add the KS transition
        ! energies to the diagonal of the block.
        if(.not. fcoup) then 
          if(iknr .eq. jknr) then
            if(fwp) then 
              call addkstransdiag(i1, ham(i1:i2,j1:j2), diag(i1:i2,1))
            else
              call addkstransdiag(i1, ham(i1:i2,j1:j2))
            end if
          end if
        end if

      ! ikkp loop end
      end do

      !write(*,*) "ham(1,1)", ham(1,1)
      !write(*,*) "dble(ham(1,1))", dble(ham(1,1))
      !write(*,*) "real(ham(1,1),8)", real(ham(1,1),8)

      ! Write lower triangular part in case it is explicitly needed
      do i1 = 1, hamsize
        if(.not.(fcoup .and. .not. fti)) then
          ham(i1, i1) = cmplx(dble(ham(i1,i1)), 0.0d0, 8)
        end if
        do i2 = i1+1, hamsize
          if(fcoup .and. .not. fti) then 
            ham(i2,i1) = ham(i1,i2)
          else
            ham(i2,i1) = conjg(ham(i1,i2))
          end if
        end do
      end do

      !write(*,*) "ham(1,1)", ham(1,1)

      call timesec(ts1)
      write(unitout, '(" Matrix build.")')
      write(unitout, '("Timing (in seconds)	   :", f12.3)') ts1 - ts0

      ! Test output
      if(fwp) then 
        call timesec(ts0)
        if(mpiglobal%rank == 0) then 
          if(.not. fcoup) then 
            write(unitout, '("Info(setup_bse_block): Writing RR: H, W, V and E to file")')
          else
            if(fti) then 
              write(unitout, '("Info(setup_bse_block): Writing RA^ti: H, W, V to file")')
            else
              write(unitout, '("Info(setup_bse_block): Writing RA: H, W, V to file")')
            end if
          end if
        end if
        ! Write unsymmetrised
        if(fcoup) then 
          if(fti) then 
            if(usescc) call writecmplxparts('WCti_unsym', dble(wint), immat=aimag(wint))
            if(useexc) call writecmplxparts('VCti_unsym', dble(vint), immat=aimag(vint))
          else
            if(usescc) call writecmplxparts('WC_unsym', dble(wint), immat=aimag(wint))
            if(useexc) call writecmplxparts('VC_unsym', dble(vint), immat=aimag(vint))
          end if
        else
          call writecmplxparts('KS', diag)
          if(usescc) call writecmplxparts('W_unsym', dble(wint), immat=aimag(wint))
          if(useexc) call writecmplxparts('V_unsym', dble(vint), immat=aimag(vint))
        end if
        ! Make W,V hermitian (RR or RA^ti) or symmetric (RA)
        do i1 = 1, hamsize
          do i2 = i1, hamsize
            if(fcoup) then 
              if(fti) then 
                if(usescc) wint(i2,i1) = conjg(wint(i1,i2))
                if(useexc) vint(i2,i1) = conjg(vint(i1,i2))
              else
                if(usescc) wint(i2,i1) = wint(i1,i2)
                if(useexc) vint(i2,i1) = vint(i1,i2)
              end if
            else
              if(usescc) wint(i2,i1) = conjg(wint(i1,i2))
              if(useexc) vint(i2,i1) = conjg(vint(i1,i2))
            end if
          end do
        end do
        if(fcoup) then 
          if(fti) then 
            call writecmplxparts('HamCti', dble(ham), immat=aimag(ham))
            if(usescc) call writecmplxparts('WCti', dble(wint), immat=aimag(wint))
            if(useexc) call writecmplxparts('VCti', dble(vint), immat=aimag(vint))
          else
            call writecmplxparts('HamC', dble(ham), immat=aimag(ham))
            if(usescc) call writecmplxparts('WC', dble(wint), immat=aimag(wint))
            if(useexc) call writecmplxparts('VC', dble(vint), immat=aimag(vint))
          end if
        else
          call writecmplxparts('Ham', dble(ham), immat=aimag(ham))
          if(usescc) call writecmplxparts('W', dble(wint), immat=aimag(wint))
          if(useexc) call writecmplxparts('V', dble(vint), immat=aimag(vint))
        end if

        ! Order for energy
        do i1 = 1, hamsize
          do i2 = 1, hamsize
            tmp(i1, i2) = ham(ensortidx(i1), ensortidx(i2))
          end do
        end do
        if(fcoup) then 
          if(fti) then 
            call writecmplxparts('HamCti_sorted', dble(tmp), immat=aimag(tmp))
          else
            call writecmplxparts('HamC_sorted', dble(tmp), immat=aimag(tmp))
          end if
        else
          call writecmplxparts('Ham_sorted', dble(tmp), immat=aimag(tmp))
        end if
        do i1 = 1, hamsize
          tmp(i1, 1) = diag(ensortidx(i1), 1)
        end do
        diag(:,1) = tmp(:,1)
        if(.not. fcoup) then 
          call writecmplxparts('KS_sorted', diag)
        end if
        do i1 = 1, hamsize
          do i2 = 1, hamsize
            tmp(i1, i2) = wint(ensortidx(i1), ensortidx(i2))
          end do
        end do
        if(fcoup) then 
          if(fti) then 
            if(usescc) call writecmplxparts('WCti_sorted', dble(tmp), immat=aimag(tmp))
          else
            if(usescc) call writecmplxparts('WC_sorted', dble(tmp), immat=aimag(tmp))
          end if
        else
          if(usescc) call writecmplxparts('W_sorted', dble(tmp), immat=aimag(tmp))
        end if
        do i1 = 1, hamsize
          do i2 = 1, hamsize
            tmp(i1, i2) = vint(ensortidx(i1), ensortidx(i2))
          end do
        end do
        if(fcoup) then 
          if(fti) then 
            if(useexc) call writecmplxparts('V_sorted', dble(tmp), immat=aimag(tmp))
          else
            if(useexc) call writecmplxparts('VC_sorted', dble(tmp), immat=aimag(tmp))
          end if
        else
          if(useexc) call writecmplxparts('V_sorted', dble(tmp), immat=aimag(tmp))
        end if
        deallocate(tmp)
        deallocate(wint)
        deallocate(vint)
        deallocate(diag)
        call timesec(ts1)
        write(unitout, '(" Parts writtten.")')
        write(unitout, '("Timing (in seconds)	   :", f12.3)') ts1 - ts0
      end if

      contains

        subroutine setint(hamblock, oc1, oc2, scc, exc, w, v)
          complex(8), intent(out) :: hamblock(:,:)
          real(8), intent(in) :: oc1(:), oc2(:)
          complex(8), intent(in), optional :: scc(:,:)
          complex(8), intent(in), optional :: exc(:,:)
          complex(8), intent(out), optional :: w(:,:), v(:,:)
          
          integer(4) :: i, j

          if(present(exc) .and. present(scc)) then 
            if (input%xs%bse%xas) then
              do j= 1, size(hamblock,2)
                do i= 1, size(hamblock,1)
                  hamblock(i,j) = oc1(i)*oc2(j) * (exc(i,j) - scc(i,j))
                end do
              end do
            else
              do j= 1, size(hamblock,2)
                do i= 1, size(hamblock,1)
                  hamblock(i,j) = oc1(i)*oc2(j) * (2.0d0 * exc(i,j) - scc(i,j))
                !write(*,*) "i,j", i, j
                !write(*,*) "oc1,oc2", oc1(i), oc1(j)
                !write(*,*) "exc,scc", exc(i,j), scc(i,j)
                !write(*,*) "hamblock", hamblock(i,j)
                end do
              end do
            end if
          else if(present(scc)) then 
            do j= 1, size(hamblock,2)
              do i= 1, size(hamblock,1)
                hamblock(i,j) = -oc1(i)*oc2(j) * scc(i,j)
              end do
            end do
          else if(present(exc)) then 
            if (input%xs%bse%xas) then
              do j= 1, size(hamblock,2)
                do i= 1, size(hamblock,1)
                  hamblock(i,j) = oc1(i)*oc2(j) * exc(i,j)
                end do
              end do
            else
              do j= 1, size(hamblock,2)
                do i= 1, size(hamblock,1)
                  hamblock(i,j) = oc1(i)*oc2(j) * 2.0d0 * exc(i,j)
                end do
              end do
            end if
          end if

          if(present(w)) then 
            w(:,:) = scc(:,:)
          end if
          if(present(v)) then 
            v(:,:) = exc(:,:)
          end if
        end subroutine setint

        subroutine addkstransdiag(ig, hamblock, d)
          integer(4), intent(in) :: ig
          complex(8), intent(inout) :: hamblock(:,:)
          real(8), intent(out), optional :: d(:)

          integer(4) :: i

          ! Calculate ks energy differences
          do i = 1, size(hamblock,1)
            hamblock(i,i) = hamblock(i,i)&
              & + cmplx(de(ig+i-1), 0.0d0, 8)
            if(present(d)) then 
              d(i) = cmplx(de(ig+i-1), 0.0d0, 8)
            end if
          end do
        end subroutine addkstransdiag

    end subroutine setup_bse_block
    !EOC

    !! Distributed versions

    !BOP
    ! !ROUTINE: setup_bse_ti_dist
    ! !INTERFACE:
    subroutine setup_bse_ti_dist(iqmt, binfo, smat, cmat, cpmat)
    ! !USES:
      use modmpi
      use modscl
      use modinput, only: input
      use mod_constants, only: zzero, zone
      use modxs, only: unitout
      use modbse, only: hamsize
      use m_dhesolver
      use m_sqrtzmat
      use m_invertzmat
      use m_writecmplxparts
      use m_dzmatmult
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   integer(4) :: iqmt  ! Index of momentum transfer Q
    !   type(blacsinfo) :: binfo  ! Info struct for BLACS grid on which the
    !   matrices are distributed
    !
    ! In/Out:
    !   type(dzmat), optional :: smat(:,:)  ! Aux. EVP matrix S=(A-B)^{1/2}(A+B)(A-B)^{1/2}
    !   type(dzmat), optional :: cmat(:,:)  ! C = (A-B)^{1/2}
    !   type(dzmat), optional :: cpmat(:,:) ! C'= (A-B)^{-1/2}
    ! 
    ! !DESCRIPTION:
    !   The routine sets up the auxilliary EVP matrix using block cyclic
    !   distributed matrices for ScaLapack.
    !
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !EOP
    !BOC

      ! I/O
      integer(4), intent(in) :: iqmt
      type(blacsinfo), intent(in) :: binfo
      type(dzmat), intent(inout), optional :: smat
      type(dzmat), intent(inout), optional :: cmat
      type(dzmat), intent(inout), optional :: cpmat

      ! Local 
      character(*), parameter :: thisname = "setup_bse_ti_dist"
      type(dzmat) :: rrmat
      type(dzmat) :: ramat
      type(dzmat) :: auxmat

      logical :: fcheckpos, fsmat, fcmat, fcpmat

      real(8) :: ts0, ts1, t1, t0
      integer(4) :: i, j

      real(8), allocatable :: evals(:)

      ! Test writeout
      logical :: fwp 
      complex(8), allocatable :: localmat(:,:)

      fwp = input%xs%bse%writeparts
      fcheckpos = input%xs%bse%checkposdef

      fsmat = .false.
      fcmat = .false.
      fcpmat = .false.
      if(present(smat)) fsmat = .true.
      if(present(cmat)) fcmat = .true.
      if(present(cpmat)) fcpmat = .true.
      if(.not. (fsmat .or. fcmat .or. fcpmat)) then
        write(unitout, '("Error(",a,"):&
          & No matrix to build specified.")') trim(thisname)
        call terminate
      end if

      write(unitout, '("Info(setup_bse_ti_dist):&
        & Setting up distributed matrices for squared EVP")')
      call timesec(ts0)

      !===========================================================!
      ! Getting main and coupling blocks of BSE hamiltonian       !
      !===========================================================!
      ! Get RR part of BSE Hamiltonian
      write(unitout, '("Info(",a,"):&
        & Setting up RR Block of orignial BSE")') trim(thisname)
      call timesec(t0)
      call new_dzmat(rrmat, hamsize, hamsize, binfo)
      call setup_bse_block_dist(rrmat, iqmt, .false., .false., binfo)
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
      ! Get RA^{ti} part of BSE Hamiltonian
      write(unitout, '("Info(",a,"):&
        & Setting up RA^{ti} Block of orignial BSE")') trim(thisname)
      call timesec(t0)
      call new_dzmat(ramat, hamsize, hamsize, binfo)
      call setup_bse_block_dist(ramat, iqmt, .true., .true., binfo)
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
      !===========================================================!

      !===========================================================!
      ! Make combination matrices                                 !
      !===========================================================!
      write(unitout, '("Info(",a,"):&
       & Setting up RR-RA matrix")') trim(thisname)
      ! rrmat = RR - RA^{it} = A-B
      do j = 1, rrmat%ncols_loc
        do i = 1, rrmat%nrows_loc
          rrmat%za(i,j) = rrmat%za(i,j) - ramat%za(i,j)
        end do
      end do

      ! Construct S matrix - A+B only needed for S
      if(fsmat) then 
        write(unitout, '("Info(",a,"):&
         & Setting up RR+RA matrix")') trim(thisname)
        ! ramat = RR + RA^{it} = A+B = A - B + 2B
        do j = 1, rrmat%ncols_loc
          do i = 1, rrmat%nrows_loc
            ramat%za(i,j) = rrmat%za(i,j) + 2.0d0*ramat%za(i,j)
          end do
        end do
        ! Check positive definitness of (A+B) (should usually be positive)
        ! Defaults to false.
        if(fcheckpos) then 
          write(unitout, '("Info(",a,"):&
            & Checking positve definitness of RR+RA")') trim(thisname)
          call timesec(t0)
          ! Get eigenvalues of A+B
          allocate(evals(hamsize))
          call new_dzmat(auxmat, hamsize, hamsize, binfo)
          auxmat%za = ramat%za
          call dhesolver(auxmat, evals, binfo)
          if(fwp) then
            if(mpiglobal%rank == 0) then 
              write(*,*) "Writing evals for A+B"
              call writecmplxparts('apb_evals', revec=evals, veclen=size(evals))
              call barrier(callername=trim(thisname))
            else
              call barrier(callername=trim(thisname))
            end if
          end if
          if(any(evals < 0.0d0)) then 
            write(*,'("Error(",a,"):&
              & RR+RA matrix is not positive definit")') trim(thisname)
            write(*,'(E10.3)') evals
            call terminate
          end if
          deallocate(evals)
          call del_dzmat(auxmat)
          call timesec(t1)
          write(unitout, '("  RR+RA is positive definite")')
          write(unitout, '("  Time needed",f12.3,"s")') t1-t0
        else
          write(unitout, '("  RR+RA is assumed to be positive definite")')
        end if
      ! Do not construct S
      else
        ! B no longer needed if S is not to be constructed (only A-B) 
        call del_dzmat(ramat)
      end if

      ! Test writeouts
      if(fwp) then 
        if(fsmat) then 
          ! A+B
          call dzmat_send2global_root(localmat, ramat, binfo)
          if(mpiglobal%rank == 0) then
            call writecmplxparts('apb_mat', remat=dble(localmat), immat=aimag(localmat))
          end if
        end if
        ! A-B
        call dzmat_send2global_root(localmat, rrmat, binfo)
        if(mpiglobal%rank == 0) then
          call writecmplxparts('amb_mat', remat=dble(localmat), immat=aimag(localmat))
        end if
      end if
      !===========================================================!

      !===========================================================!
      ! Take the square root of (A-B)                             !
      ! Note: It is assumed to be positive definit.               !
      !===========================================================!
      write(unitout, '("Info(",a,"):&
        & Taking square root of RR-RA matrix")') trim(thisname)
      call timesec(t0)
      ! rrmat -> (A-B)^1/2
      call sqrtdzmat_hepd(rrmat, binfo, eecs=input%xs%bse%eecs)
      call timesec(t1)
      write(unitout, '("  RR-RA is positive definite")')
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
      ! Test writeout
      if(fwp) then 
        call dzmat_send2global_root(localmat, rrmat, binfo)
        if(mpiglobal%rank == 0) then
          call writecmplxparts('sqrtamb_mat', remat=dble(localmat), immat=aimag(localmat))
        end if
      end if
      !===========================================================!

      if(fsmat) then 
        !===========================================================!
        ! Construct S Matrix                                        !
        ! S = (A-B)^{1/2} (A+B) (A-B)^{1/2}                         !
        !===========================================================!
        write(unitout, '("Info(",a,"):&
          & Constructing S matrix")') trim(thisname)
        call timesec(t0)
        !   C = (A-B)^{1/2} (A+B)
        call new_dzmat(auxmat, hamsize, hamsize, binfo)
        call dzmatmult(rrmat, ramat, auxmat)
        call del_dzmat(ramat)
        !   S = C (A-B)^{1/2}
        call new_dzmat(smat, hamsize, hamsize, binfo)
        call dzmatmult(auxmat, rrmat, smat)
        call del_dzmat(auxmat)
        call timesec(t1)
        write(unitout, '("  Time needed",f12.3,"s")') t1-t0
        ! Test writeout
        if(fwp) then 
          call dzmat_send2global_root(localmat, smat, binfo)
          if(mpiglobal%rank == 0) then
            call writecmplxparts('s_mat', remat=dble(localmat), immat=aimag(localmat))
          end if
        end if
        !===========================================================!
      end if

      if(fcmat) then
        !===========================================================!
        ! Make C Matrix                                             !
        ! Cmat = (A-B)^1/2                                          !
        ! Cmat is needed to build Oscillator strengths              !
        !===========================================================!
        write(unitout, '("Info(",a,"):&
          & Retruning C=(RR-RA)^1/2 matrix")') trim(thisname)
        call new_dzmat(cmat, hamsize, hamsize, binfo)
        do j = 1, rrmat%ncols_loc
          do i = 1, rrmat%nrows_loc
            cmat%za(i,j) = rrmat%za(i,j)
          end do
        end do
        !===========================================================!
      end if

      if(fcpmat) then 
        !===========================================================!
        ! Make C' Matrix                                            !
        ! Cpmat = (A-B)^{-1/2}                                      !
        ! Cpmat is needed additionally to build eigenvectors        !
        !===========================================================!
        write(unitout, '("Info(",a,"):&
          & Inverting C=(RR-RA)^1/2 matrix")') trim(thisname)
        call timesec(t0)
        call dzinvert(rrmat)
        call new_dzmat(cpmat, hamsize, hamsize, binfo)
        do j = 1, rrmat%ncols_loc
          do i = 1, rrmat%nrows_loc
            cpmat%za(i,j) = rrmat%za(i,j)
          end do
        end do
        call timesec(t1)
        write(unitout, '("  Time needed",f12.3,"s")') t1-t0
        !===========================================================!
      end if

      call del_dzmat(rrmat)
      call timesec(ts1)
      write(unitout, '("Info(",a,"):&
        & Total time needed",f12.3,"s")') trim(thisname), ts1-ts0

    end subroutine setup_bse_ti_dist
    !EOC

    !BOP
    ! !ROUTINE: setup_bse_block_dist
    ! !INTERFACE:
    subroutine setup_bse_block_dist(ham, iqmt, fcoup, fti, binfo)
    ! !USES:
      use modmpi
      use modscl
      use modinput, only: input
      use mod_constants, only: zzero, zone
      use modxs, only: unitout
      use modbse, only: nou_bse_max, kmap_bse_rg,&
                      & kousize, nkkp_bse, nk_bse, ofac,&
                      & scclifbasename, exclifbasename,&
                      & scclicfbasename, exclicfbasename,&
                      & scclictifbasename, infofbasename
      use m_getunit
      use m_genfilname
      use m_putgetbsemat
      use m_writecmplxparts
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   integer(4) :: iqmt ! Index of momentum transfer Q
    !   logical :: fcoup   ! If true, builds RA instead of RR block of BSE matrix 
    !   logical :: fti     ! If true, use time inverted antiresonant basis
    !   type(blacsinfo) :: binfo ! Info type of the BLACS grid
    ! In/Out:
    !   type(dzmat) :: ham ! 2D block cyclic distributed RR or RA
    !                      ! block of BSE-Hamiltonian matrix
    ! 
    ! !DESCRIPTION:
    !   The routine sets up the content of the local array of the 
    !   2d block cyclic distributed resonant-resonant or resonant-antiresonant block of
    !   the BSE-Hamiltonian matrix. Process 0 reads {\tt EXCLI.OUT} and {\tt SCCLI.OUT}
    !   ({\tt EXCLIC.OUT} and {\tt SCCLIC.OUT} or {\tt SCCLICTI.OUT}) for each {\tt ikkp}
    !   record and send the data block-wise to the responsible processes.
    !   If the matrix is not to be distributed the routine calls 
    !   {\tt setup\_bse} instead.
    !
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !EOP
    !BOC

      implicit none

      !! I/O
      type(dzmat), intent(inout) :: ham
      integer(4), intent(in) :: iqmt
      logical, intent(in) :: fcoup, fti
      type(blacsinfo), intent(in) :: binfo

      !! Local variables
      integer(4) :: ikkp, ik, jk, iknr, jknr
      integer(4) :: inou, jnou

      ! Arrays to read V and W from file
      complex(8), allocatable, dimension(:,:) :: excli_t, sccli_t

      ! Send buffers for distributing V and W 
      complex(8), allocatable, dimension(:,:) :: ebuff, sbuff
      complex(8), allocatable, dimension(:,:) :: ebuff2, sbuff2

      ! BLACS context
      integer(4) :: context
      ! Process row and column coordinates
      integer(4) :: pr, pc
      integer(4) :: pr2, pc2

      ! Row and column index of ham (Global view)
      integer(4) :: ig, jg
      ! Row and column block length of ham
      integer(4) :: iblck, jblck
      ! Position of ikjk sub matrix within ham (Global view)
      integer(4) :: ii, jj
      ! Block lengths of sub matrix
      integer(4) :: ib, jb
      ! Row and column index of the sub matrix (Global view)
      integer(4) :: i, j
      ! Row and column index (Local view)
      integer(4) :: il, jl
      integer(4) :: il2, jl2

      ! Read in vars
      character(256) :: sfname, sinfofname
      character(256) :: efname, einfofname
      logical :: efcmpt, efid
      logical :: sfcmpt, sfid

      logical :: useexc, usescc

      logical :: fwp

      fwp = input%xs%bse%writeparts

      if(ham%isdistributed) then 

#ifdef SCAL
        !!*****************************!!
        !! ONLY PROCESS 0 READS DATA   !!
        !! AND SENDS THE CORRESPONDING !!
        !! DATA CHUNKS TO THE OTHERS.  !!
        !!*****************************!!

        usescc = .false.
        useexc = .false.
        select case(trim(input%xs%bse%bsetype))
          case('IP')
            usescc = .false.
            useexc = .false.
            if(mpiglobal%rank==0) then
              write(unitout, '("Info(setup_bse_block_dist): Using IP")')
            end if
          case('singlet')
            usescc = .true.
            useexc = .true.
            if(mpiglobal%rank==0) then
              write(unitout, '("Info(setup_bse_block_dist): Using singlet")')
            end if
          case('triplet')
            usescc = .true.
            useexc = .false.
            if(mpiglobal%rank==0) then
              write(unitout, '("Info(setup_bse_block_dist): Using triplet")')
            end if
          case('RPA')
            usescc = .false.
            useexc = .true.
            if(mpiglobal%rank==0) then
              write(unitout, '("Info(setup_bse_block_dist): Using RPA")')
            end if
        end select

        ! Check whether the saved data is usable and 
        ! allocate read in work arrays
        if(binfo%isroot) then 

          if(.not. fcoup) then 
            write(unitout, '("Info(setup_bse_block_dist):&
              & Setting up RR part of hamiltonian")')
          else
            write(unitout, '("Info(setup_bse_block_dist):&
              & Setting up RA part of hamiltonian")')
            if(fti) then 
              write(unitout, '("Info(setup_bse_block_dist):&
                & Using time inverted anti-resonant basis")')
            end if
          end if

          ! Select form which files W and V are to be read
          if(.not. fcoup) then
            call genfilname(basename=scclifbasename, iqmt=iqmt, filnam=sfname)
            call genfilname(basename=exclifbasename, iqmt=iqmt, filnam=efname)
          else
            if(fti) then
              call genfilname(basename=scclictifbasename, iqmt=iqmt, filnam=sfname)
              call genfilname(basename=exclifbasename, iqmt=iqmt, filnam=efname)
            else
              call genfilname(basename=scclicfbasename, iqmt=iqmt, filnam=sfname)
              call genfilname(basename=exclicfbasename, iqmt=iqmt, filnam=efname)
            end if
          end if

          sinfofname = trim(infofbasename)//'_'//trim(sfname)
          einfofname = trim(infofbasename)//'_'//trim(efname)

          if(usescc) write(unitout, '("  Reading form info from ", a)')&
            & trim(sinfofname)
          if(useexc) write(unitout, '("  Reading form info from ", a)')&
            & trim(einfofname)

          ! Check saved quantities for compatibility
          sfcmpt=.false.
          sfid=.false.
          if(usescc) call b_getbseinfo(trim(sinfofname), iqmt,&
            & fcmpt=sfcmpt, fid=sfid)
          efcmpt=.false.
          efid=.false.
          if(useexc) call b_getbseinfo(trim(einfofname), iqmt,&
            & fcmpt=efcmpt, fid=efid)

          if(usescc .and. useexc) then 
            if(efcmpt .neqv. sfcmpt .or. efid .neqv. sfid) then 
              write(*, '("Error(setup_bse_block_dist): Info files differ")')
              write(*, '("  efcmpt, efid, sfcmpt, sfcmpt")') efcmpt, efid, sfcmpt, sfcmpt
              call terminate
            end if
          end if

          if(usescc) write(unitout, '("  Reading form W from ", a)') trim(sfname)
          if(usescc) write(unitout, '("  compatible:",l," identical:",l)') sfcmpt, sfid
          if(useexc) write(unitout, '("  Reading form V from ", a)') trim(efname)
          if(usescc) write(unitout, '("  compatible:",l," identical:",l)') efcmpt, efid

          allocate(excli_t(nou_bse_max, nou_bse_max))
          allocate(sccli_t(nou_bse_max, nou_bse_max))

        end if

        ! Block sizes of global matrix
        iblck = ham%mblck
        jblck = ham%nblck 

        ! Context
        context = ham%context
        if(context /= binfo%context) then 
          write(*, '("Error(setup_bse_block_dist): Context mismatch:", 2i4)')&
            & context, binfo%context
          call terminate
        end if

        ! Send/receive buffer
        !   Upper triangular part
        allocate(ebuff(iblck, jblck))
        allocate(sbuff(iblck, jblck))
        !   Lower triangular part 
        !   (ham is either hermitian or symmetric)
        allocate(ebuff2(jblck, iblck))
        allocate(sbuff2(jblck, iblck))

        ! Loop over ikkp blocks of the global 
        ! BSE Hamilton matrix.
        ! Note: Since W and V are hermitian or symmetric 
        !       only ikjk blocks are saved with jk >= ik
        !       so ikkp labels ikjk blocks of the upper triangular part of the matrix
        do ikkp = 1, nkkp_bse

          ! Get index ik jk form combination index ikkp
          call kkpmap(ikkp, nk_bse, ik, jk)

          ! Get total non reduced k point index iknr (jknr) form 
          ! BSE k point selection index ik (jk)
          iknr = kmap_bse_rg(ik)
          jknr = kmap_bse_rg(jk)

          ! Position of ikkp block in global matrix
          ii = sum(kousize(1:iknr-1)) + 1
          jj = sum(kousize(1:jknr-1)) + 1

          ! Size of sub-matrix
          inou = kousize(iknr)
          jnou = kousize(jknr)

          !!***********!!
          !! READ DATA !!
          !!***********!!
          ! Root reads a ikjk block of W and V from file
          if(binfo%isroot) then

            ! Read in screened coulomb interaction for ikkp
            if(usescc) then 
              ! Read RR/RA part of screened coulomb interaction W_{iuioik,jujojk}(qmt)
              call b_getbsemat(trim(sfname), iqmt, ikkp,&
                & sccli_t(1:inou,1:jnou), check=.false., fcmpt=sfcmpt, fid=sfid)
            end if

            ! Read in exchange interaction for ikkp
            if(useexc) then 
              ! Read RR/RA part of exchange interaction v_{iuioik,jujojk}(qmt)
              call b_getbsemat(trim(efname), iqmt, ikkp,&
                & excli_t(1:inou,1:jnou), check=.false., fcmpt=efcmpt, fid=efid)
            end if

            if(fwp) then
              if(ikkp == 1 .and. fcoup .eqv. .false.) then 
                if(usescc) then 
                  call writecmplxparts('ikkp1_read_Wrr',dble(sccli_t(1:inou,1:jnou)), aimag(sccli_t(1:inou,1:jnou)))
                end if
                if(useexc) then 
                  call writecmplxparts('ikkp1_read_Vrr',dble(excli_t(1:inou,1:jnou)), aimag(excli_t(1:inou,1:jnou)))
                end if
              end if

              if(ikkp == 1 .and. fcoup .eqv. .true. .and. fti .eqv. .false.) then 
                if(usescc) then 
                  call writecmplxparts('ikkp1_read_Wra',dble(sccli_t(1:inou,1:jnou)), aimag(sccli_t(1:inou,1:jnou)))
                end if
                if(useexc) then 
                  call writecmplxparts('ikkp1_read_Vra',dble(excli_t(1:inou,1:jnou)), aimag(excli_t(1:inou,1:jnou)))
                end if
              end if

              if(ikkp == 1 .and. fcoup .eqv. .true. .and. fti .eqv. .true.) then 
                if(usescc) then 
                  call writecmplxparts('ikkp1_read_Wra_ti',dble(sccli_t(1:inou,1:jnou)), aimag(sccli_t(1:inou,1:jnou)))
                end if
                if(useexc) then 
                  call writecmplxparts('ikkp1_read_Vra_ti',dble(excli_t(1:inou,1:jnou)), aimag(excli_t(1:inou,1:jnou)))
                end if
              end if
            end if

            ! Make ik=jk blocks explicitly symmetric/hermitian
            if(iknr == jknr) then 
              if(fcoup .and. .not. fti) then 
                do j = 1, jnou
                  do i = 1, j
                    if(usescc) then 
                      sccli_t(j,i) = sccli_t(i,j)
                    end if
                    if(useexc) then
                      excli_t(j,i) = excli_t(i,j)
                    end if
                  end do
                end do
              else
                do j = 1, jnou
                  do i = 1, j
                    if(usescc) then
                      sccli_t(j,i) = conjg(sccli_t(i,j))
                    end if
                    if(useexc) then 
                      excli_t(j,i) = conjg(excli_t(i,j))
                    end if
                  end do
                end do
              end if
            end if

          end if

          !!******************!!
          !! SEND DATA CHUNKS !!
          !!******************!!
          ! Root distributes blocks of the ikjk matrix to all processes
          ! according to the block cyclic distribution scheme of ScaLAPACK.
          ! Note:
          ! Also the jkik block is distributed explicitly to avoid later 
          ! complications, in case that the lower triangular part
          ! of the matrix is also needed.

          ! Column index of global ikkp sub-matrix
          j = 1
          do while(j <= jnou)

            ! Column index of global matrix
            jg = jj + j - 1

            ! Calculate column block size
            if(j == 1) then 
              ! First column block size of global sub-matrix
              ! Adjust for possible truncation of first column block size 
              jb = jblck - mod(jj-1, jblck)
              jb = min(jb, jnou-j+1)
            else
              ! Adjust for possible truncation of last column block size 
              jb = min(jblck, jnou-j+1)
            end if


            !! We have a sub block of the 
            !! global sub-matrix at col j of colsize jb
            !! that needs to be sent do one process only.
            ! Get process grid column coordinate of responsible process.
            pc = indxg2p( jg, jblck, binfo%mypcol, 0, binfo%npcols)
            ! Get column position of ib*jb block in local matrix
            jl = indxg2l( jg, jblck, pc, 0, binfo%npcols)

            ! Get corresponding quantities in the lower 
            ! triangular part of the matrix.
            if(iknr /= jknr) then 
              ! Get process grid row coordinate
              pr2 = indxg2p( jg, iblck, binfo%myprow, 0, binfo%nprows)
              ! Get row position of jb*ib block in local matrix
              il2 = indxg2l( jg, iblck, pr2, 0, binfo%nprows)

            end if

            ! Row index of global sub-matrix
            i = 1
            do while(i <= inou)

              ! Row index of global matrix
              ig = ii + i - 1

              if( i == 1 ) then 
                ! First row block size of global sub-matrix
                ! Adjust for possible truncation of first row block size 
                ib = iblck - mod(ii-1, iblck)
                ib = min(ib, inou-i+1)
              else
                ! Adjust for possible truncation of last row block size 
                ib = min(iblck, inou-i+1)
              end if

              !! We have a sub block of the 
              !! global sub-matrix at coordinates i,j of size ib*jb
              !! that needs to be sent do one process only.
              ! Get process grid row coordinate of responsible process.
              pr = indxg2p( ig, iblck, binfo%myprow, 0, binfo%nprows)
              ! Get row position of ib*jb block in local matrix
              il = indxg2l( ig, iblck, pr, 0, binfo%nprows)

              ! Get corresponding quantities in the lower 
              ! triangular part of the matrix.
              if(iknr /= jknr) then 
                ! Get process grid column coordinate
                pc2 = indxg2p( ig, jblck, binfo%mypcol, 0, binfo%npcols)
                ! Get column position of jb*ib block in local matrix
                jl2 = indxg2l( ig, jblck, pc2, 0, binfo%npcols)

              end if


              ! Root does data sending
              if(binfo%isroot) then 

                ! Prepare send packages
                if(usescc) sbuff(1:ib,1:jb) = sccli_t(i:i+ib-1, j:j+jb-1)
                if(useexc) ebuff(1:ib,1:jb) = excli_t(i:i+ib-1, j:j+jb-1)

                if(iknr /= jknr) then 
                  ! Also build lower part of hermitian/symmetric matrix
                  if(fcoup .and. .not. fti) then 
                    if(usescc) sbuff2(1:jb,1:ib) = transpose(sccli_t(i:i+ib-1, j:j+jb-1))
                    if(useexc) ebuff2(1:jb,1:ib) = transpose(excli_t(i:i+ib-1, j:j+jb-1))
                  else
                    if(usescc) sbuff2(1:jb,1:ib) = conjg(transpose(sccli_t(i:i+ib-1, j:j+jb-1)))
                    if(useexc) ebuff2(1:jb,1:ib) = conjg(transpose(excli_t(i:i+ib-1, j:j+jb-1)))
                  end if
                end if

                ! No send needed, root is responsible for that block

                ! Upper triangular part
                if( pr == 0 .and. pc == 0) then

                  ! Singlet
                  if(useexc .and. usescc) then 
                    call buildham(fcoup, ham%za(il:il+ib-1, jl:jl+jb-1),&
                      & ig, jg, ib, jb,&
                      & occ1=ofac(ig:ig+ib-1), occ2=ofac(jg:jg+jb-1),&
                      & scc=sbuff(1:ib,1:jb), exc=ebuff(1:ib,1:jb))
                  ! Triplet
                  else if(usescc) then 
                    call buildham(fcoup, ham%za(il:il+ib-1, jl:jl+jb-1),&
                      & ig, jg, ib, jb,&
                      & occ1=ofac(ig:ig+ib-1), occ2=ofac(jg:jg+jb-1),&
                      & scc=sbuff(1:ib,1:jb))
                  ! RPA
                  else if(useexc) then
                    call buildham(fcoup, ham%za(il:il+ib-1, jl:jl+jb-1),&
                      & ig, jg, ib, jb,&
                      & occ1=ofac(ig:ig+ib-1), occ2=ofac(jg:jg+jb-1),&
                      & exc=ebuff(1:ib,1:jb))
                  ! IP
                  else
                    call buildham(fcoup, ham%za(il:il+ib-1, jl:jl+jb-1),&
                      & ig, jg, ib, jb,&
                      & occ1=ofac(ig:ig+ib-1), occ2=ofac(jg:jg+jb-1))
                  end if

                ! Send data 
                else

                  if(usescc) call zgesd2d(context, ib, jb, sbuff, iblck, pr, pc)
                  if(useexc) call zgesd2d(context, ib, jb, ebuff, iblck, pr, pc)

                end if

                ! Lower triangular part
                if(iknr /= jknr) then 

                  if( pr2 == 0 .and. pc2 == 0) then 

                    if(useexc .and. usescc) then 
                      call buildham(fcoup, ham%za(il2:il2+jb-1, jl2:jl2+ib-1),&
                        & jg, ig, jb, ib,&
                        & occ1=ofac(jg:jg+jb-1), occ2=ofac(ig:ig+ib-1),&
                        & scc=sbuff2(1:jb,1:ib), exc=ebuff2(1:jb,1:ib))
                    else if(usescc) then
                      call buildham(fcoup, ham%za(il2:il2+jb-1, jl2:jl2+ib-1),&
                        & jg, ig, jb, ib,&
                        & occ1=ofac(jg:jg+jb-1), occ2=ofac(ig:ig+ib-1),&
                        & scc=sbuff2(1:jb,1:ib))
                    else if(useexc) then
                      call buildham(fcoup, ham%za(il2:il2+jb-1, jl2:jl2+ib-1),&
                        & jg, ig, jb, ib,&
                        & occ1=ofac(jg:jg+jb-1), occ2=ofac(ig:ig+ib-1),&
                        & exc=ebuff2(1:jb,1:ib))
                    else
                      call buildham(fcoup, ham%za(il2:il2+jb-1, jl2:jl2+ib-1),&
                        & jg, ig, jb, ib,&
                        & occ1=ofac(jg:jg+jb-1), occ2=ofac(ig:ig+ib-1))
                    end if

                  ! Send data 
                  else

                    if(usescc) call zgesd2d(context, jb, ib, sbuff2, jblck, pr2, pc2)
                    if(useexc) call zgesd2d(context, jb, ib, ebuff2, jblck, pr2, pc2)

                  end if

                end if

              ! All others only receive
              else
               
                ! Upper triangular part
                if(pr == binfo%myprow  .and. pc == binfo%mypcol) then

                  ! Receive block
                  if(usescc) call zgerv2d(context, ib, jb, sbuff, iblck, 0, 0)
                  if(useexc) call zgerv2d(context, ib, jb, ebuff, iblck, 0, 0)

                  ! Assemble sub-block in local Hamilton matrix
                  if(useexc .and. usescc) then 
                    call buildham(fcoup, ham%za(il:il+ib-1, jl:jl+jb-1),&
                      & ig, jg, ib, jb,&
                      & occ1=ofac(ig:ig+ib-1), occ2=ofac(jg:jg+jb-1),&
                      & scc=sbuff(1:ib,1:jb), exc=ebuff(1:ib,1:jb))
                  else if(usescc) then
                    call buildham(fcoup, ham%za(il:il+ib-1, jl:jl+jb-1),&
                      & ig, jg, ib, jb,&
                      & occ1=ofac(ig:ig+ib-1), occ2=ofac(jg:jg+jb-1),&
                      & scc=sbuff(1:ib,1:jb))
                  else if(useexc) then
                    call buildham(fcoup, ham%za(il:il+ib-1, jl:jl+jb-1),&
                      & ig, jg, ib, jb,&
                      & occ1=ofac(ig:ig+ib-1), occ2=ofac(jg:jg+jb-1),&
                      & exc=ebuff(1:ib,1:jb))
                  else
                    call buildham(fcoup, ham%za(il:il+ib-1, jl:jl+jb-1),&
                      & ig, jg, ib, jb,&
                      & occ1=ofac(ig:ig+ib-1), occ2=ofac(jg:jg+jb-1))
                  end if

                end if

                ! Lower triangular part
                if(iknr /= jknr) then

                  if(pr2 == binfo%myprow  .and. pc2 == binfo%mypcol) then

                    ! Receive block
                    if(usescc) call zgerv2d(context, jb, ib, sbuff2, jblck, 0, 0)
                    if(useexc) call zgerv2d(context, jb, ib, ebuff2, jblck, 0, 0)

                    ! Lower triangular part
                    if(useexc .and. usescc) then 
                      call buildham(fcoup, ham%za(il2:il2+jb-1, jl2:jl2+ib-1),&
                        & jg, ig, jb, ib,&
                        & occ1=ofac(jg:jg+jb-1), occ2=ofac(ig:ig+ib-1),&
                        & scc=sbuff2(1:jb,1:ib), exc=ebuff2(1:jb,1:ib))
                    else if(usescc) then
                      call buildham(fcoup, ham%za(il2:il2+jb-1, jl2:jl2+ib-1),&
                        & jg, ig, jb, ib,&
                        & occ1=ofac(jg:jg+jb-1), occ2=ofac(ig:ig+ib-1),&
                        & scc=sbuff2(1:jb,1:ib))
                    else if(useexc) then
                      call buildham(fcoup, ham%za(il2:il2+jb-1, jl2:jl2+ib-1),&
                        & jg, ig, jb, ib,&
                        & occ1=ofac(jg:jg+jb-1), occ2=ofac(ig:ig+ib-1),&
                        & exc=ebuff2(1:jb,1:ib))
                    else
                      call buildham(fcoup, ham%za(il2:il2+jb-1, jl2:jl2+ib-1),&
                        & jg, ig, jb, ib,&
                        & occ1=ofac(jg:jg+jb-1), occ2=ofac(ig:ig+ib-1))
                    end if

                  end if

                end if

              end if

              ! Next row block
              i = i + ib

              call blacs_barrier(binfo%context, 'A')

            ! i while loop
            end do

            ! Next column block
            j = j + jb

            call blacs_barrier(binfo%context, 'A')

          ! j while loop
          end do

          call blacs_barrier(binfo%context, 'A')

        ! ikkp loop
        end do
#else
        write(*,*) "Error(setup_bse_block_dist): Scalapack needed."
        call terminate
#endif

      else

        call setup_bse_block(ham%za, iqmt, fcoup, fti)

      end if

    end subroutine setup_bse_block_dist
    !EOC

    !BOP
    ! !ROUTINE: buildham
    ! !INTERFACE:
    subroutine buildham(fc, hamblck, ig, jg, ib, jb, occ1, occ2, scc, exc)
    ! !USES:
      use modbse, only: de
      use modinput, only: input
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! logical :: fc ! Build RA instead of RR
    ! integer(4) :: ig, jg ! Position of sub block in global matrix
    ! integer(4) :: ib, jb ! Sub block size
    ! real(8)    :: occ1(ib), occ2(jb)    ! Occupation factors
    ! complex(8), optional :: scc(ib, jb) ! Screened Coulomb interaction 
    ! complex(8), optional :: exc(ib, jb) ! Exchange interaction 
    ! In/Out:
    ! complex(8) :: hamblck(ib,jb)       ! Sub block of BSE-Hamiltonian
    !
    ! !DESCRIPTION:
    !   The routine returns a sub block of the distributed BSE-Hamiltonian matrix:\\
    !   $H(i_g:ig+ib-1, j_g:j_g+jb-1)$ where each entry is computed according to \\
    !   $H(i, j) = E(i, j) + F(i) \left( - W(i, j) + 2*V(i, j) \right) F(j)$\\
    !   Only if the sub block contains diagonal elements of the matrix the kohn sham
    !   transition energies $E$ will be added (RR case only).
    !   From the transition energies the gap energy is subtracted and
    !   the scissor is added. The exchange term $V$ is added optionally. \\
    !
    !   The matrix indices correspond to combined indices $\alpha$: \\
    !   Where $\alpha = \{ \vec{k}_\alpha, o_\alpha, u_\alpha \}$, so that: \\
    !   $F_{\alpha} = \sqrt{ \left| f_{\vec{k}_{\alpha_1} o_{\alpha_1}} 
    !                        - f_{\vec{k}_{\alpha_1} u_{\alpha_1}} \right|}$ \\
    !   $V_{\alpha_1,\alpha_2} = 
    !    V_{\vec{k}_{\alpha_1} o_{\alpha_1} u_{\alpha_1},
    !        \vec{k}_{\alpha_2} o_{\alpha_2} u_{\alpha_2}}$ \\
    !   $W_{\alpha_1,\alpha_2} = 
    !    W_{\vec{k}_{\alpha_1} o_{\alpha_1} u_{\alpha_1},
    !        \vec{k}_{\alpha_2} o_{\alpha_2} u_{\alpha_2}}$ 
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
      integer(4), intent(in) :: ig, jg, ib, jb
      logical, intent(in) :: fc
      real(8), intent(in) :: occ1(ib), occ2(jb)
      complex(8), intent(in), optional :: scc(ib,jb)
      complex(8), intent(in), optional :: exc(ib,jb)
      complex(8), intent(inout) :: hamblck(ib,jb)
      
      integer(4) :: r, c
      complex(8), parameter :: ztwo=(2.0d0,0.0d0) 
      complex(8), parameter :: zone=(1.0d0,0.0d0) 
      
      do c = 1, jb
        do r = 1, ib
          if(present(exc) .and. present(scc)) then
            ! Singlet case with exchange interaction
            if(input%xs%bse%xas) then
              hamblck(r, c) = occ1(r) * (zone * exc(r, c) - scc(r, c)) * occ2(c)
            else
              hamblck(r, c) = occ1(r) * (ztwo * exc(r, c) - scc(r, c)) * occ2(c)
            end if
          else if(present(scc)) then
            ! Triplet case without exchange interaction
            hamblck(r, c) = -occ1(r) * scc(r, c) * occ2(c)
          else if(present(exc)) then
            ! RPA
            if(input%xs%bse%xas) then 
              hamblck(r, c) = occ1(r) * zone * exc(r, c) * occ2(c)
            else
              hamblck(r, c) = occ1(r) * ztwo * exc(r, c) * occ2(c)
            end if
          else
            ! IP
            hamblck(r, c) = (0.0d0,0.0d0)
          end if
          ! Add KS transition energies
          if(.not. fc) then 
            if(ig+r-1 == jg+c-1) then 
              hamblck(r, c) = hamblck(r, c) + cmplx(de(ig+r-1), 0.0d0, 8)
            end if
          end if
        end do
      end do

    end subroutine buildham
    !EOC

end module m_setup_bse
