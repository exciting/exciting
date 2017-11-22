module m_setup_bse
      
  implicit none

  private

  public :: setup_bse_tr, setup_bse_block
  public :: setup_bse_tr_dist, setup_bse_block_dist

  contains

    !! Serial versions

    !BOP
    ! !ROUTINE: setup_bse_tr
    ! !INTERFACE:
    subroutine setup_bse_tr(iqmt, smat, cmat, cpmat)
    ! !USES:
      use modmpi
      use modinput, only: input
      use mod_constants, only: zzero, zone
      use modbse, only: hamsize
      use modxs, only: unitout
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
      character(*), parameter :: thisname = "setup_bse_tr"
      complex(8), allocatable :: rrmat(:, :)
      complex(8), allocatable :: ramat(:, :)
      complex(8), allocatable :: auxmat(:, :)


      real(8) :: ts0, ts1, t1, t0
      integer(4) :: i, j

      real(8), allocatable :: evals(:)

      logical :: fcheckpos, fsmat, fcmat, fcpmat

      ! Check if A+B is positive definite
      fcheckpos = input%xs%bse%checkposdef

      ! What to build
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

      write(unitout, '("Info(setup_bse_tr): Setting up matrices for squared EVP")')
      call timesec(ts0)

      !===========================================================!
      ! Getting main and coupling blocks of BSE hamiltonian       !
      !===========================================================!
      ! Get RR part of BSE Hamiltonian

      write(unitout, '("Info(setup_bse_tr): Setting up RR Block of orignial BSE")')
      call timesec(t0)
      allocate(rrmat(hamsize,hamsize))
      call setup_bse_block(rrmat, iqmt, .false.)
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0

      ! Get RA^{ti} part of BSE Hamiltonian
      write(unitout, '("Info(setup_bse_tr):&
        & Setting up RA^{ti} Block of orignial BSE")')
      call timesec(t0)
      allocate(ramat(hamsize,hamsize))
      call setup_bse_block(ramat, iqmt, .true.)
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
          write(unitout, '("Info(setup_bse_tr): Checking positve definitness of RR+RA")')
          call timesec(t0)
          allocate(auxmat(hamsize,hamsize))
          auxmat = ramat
          allocate(evals(hamsize))
          call hesolver(auxmat, evals)
          if(any(evals < 0.0d0)) then 
            write(*,*) "Error(setup_bse_tr): A+B matrix is not positive definit"
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

      !===========================================================!

      !===========================================================!
      ! Take the square root of (A-B)                             !
      ! Note: It is assumed to be positive definit.               !
      !===========================================================!
      write(unitout, '("Info(setup_bse_tr): Taking square root of RR-RA matrix")')
      call timesec(t0)
      ! rrmat -> (A-B)^1/2
      call sqrtzmat_hepd(rrmat)
      call timesec(t1)
      write(unitout, '("  RR-RA is positive definite")')
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
      !===========================================================!

      if(fsmat) then 
        !===========================================================!
        ! Construct S Matrix                                        !
        ! S = (A-B)^{1/2} (A+B) (A-B)^{1/2}                         !
        !===========================================================!

        write(unitout, '("Info(setup_bse_tr): Constructing S matrix")')
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
      write(unitout, '("Info(setup_bse_tr): Total time needed",f12.3,"s")') ts1-ts0

    end subroutine setup_bse_tr
    !EOC

    !BOP
    ! !ROUTINE: setup_bse_block
    ! !INTERFACE:
    subroutine setup_bse_block(ham, iqmt, fcoup)
    ! !USES:
      use modmpi
      use modinput, only: input
      use mod_constants, only: zzero, zone
      use modxs, only: unitout
      use modbse, only: de, nou_bse_max, hamsize, kmap_bse_rg,&
                      & kousize, nkkp_bse, nk_bse, ofac,&
                      & scclifbasename, exclifbasename,&
                      & scclicfbasename,&
                      & infofbasename,&
                      & vwdiffrr, vwdiffar
      use m_getunit
      use m_genfilname
      use m_putgetbsemat
      use m_writecmplxparts
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   integer(4) :: iqmt  ! Index of momentum transfer Q
    !   logical    :: fcoup ! If true, builds RA instead of RR block of BSE matrix 
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
      logical, intent(in) :: fcoup

      !! Local variables
      ! Indices
      integer(4) :: ikkp, ik, jk, iknr, jknr
      integer(4) :: inou, jnou
      integer(4) :: i1, i2, j1, j2
      ! Work arrays
      complex(8), dimension(nou_bse_max, nou_bse_max) :: excli_t, sccli_t

      ! Read in vars
      character(256) :: sfname, sinfofname
      character(256) :: efname, einfofname
      logical :: efcmpt, efid
      logical :: sfcmpt, sfid
      logical :: esfcmpt, esfid, ferror
      logical :: useexc, usescc

      ! Timings
      real(8) :: ts0, ts1
      
      logical :: fchibarq
      logical :: fmeasure

      ! Write out the coupling measure ? 
      fmeasure = input%xs%bse%measure

      ! Use truncated Coulomb potential
      fchibarq = input%xs%bse%chibarq

      call timesec(ts0)
      if(.not. fcoup) then 
        write(unitout, '("Info(setup_bse_block): Setting up RR part of hamiltonian")')
      else
        write(unitout, '("Info(setup_bse_block): Setting up RA part of hamiltonian")')
        write(unitout, '("Info(setup_bse_block):&
          & Using time reversal symmetry.")')
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

      ! Zero ham
      ham = zzero

      ! Select form which files W and V are to be read
      if(.not. fcoup) then
        call genfilname(basename=scclifbasename, iqmt=iqmt, filnam=sfname)
      else
        call genfilname(basename=scclicfbasename, iqmt=iqmt, filnam=sfname)
      end if
      if(.not. fchibarq) then
        call genfilname(basename=exclifbasename, iqmt=iqmt, filnam=efname)
      else
        call genfilname(basename=exclifbasename, bsetype="-BAR", iqmt=iqmt, filnam=efname)
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
      if(usescc) call getbseinfo(trim(sinfofname), iqmt,&
        & fcmpt=sfcmpt, fid=sfid)
      efcmpt=.false.
      efid=.false.
      if(useexc) call getbseinfo(trim(einfofname), iqmt,&
        & fcmpt=efcmpt, fid=efid)

      if(usescc .and. useexc) then 
        esfcmpt = efcmpt .and. sfcmpt
        esfid = efid .and. sfid
        ferror = .not. (esfcmpt .and. esfid)
        if(ferror) then 
          write(*, '("Error(setup_bse_block): Info files differ")')
          write(*, '("  efcmpt, efid, sfcmpt, sfcmpt", 4l)') efcmpt, efid, sfcmpt, sfcmpt
          write(*, '("  efcmpt, sfcmpt, efid, sfid, esfcmpt, esfid", 6l)') efcmpt, sfcmpt, efid, sfid, esfcmpt, esfid
          call terminate
        end if
      end if

      if(usescc) write(unitout, '("  Reading form W from ", a)') trim(sfname)
      if(usescc) write(unitout, '("  compatible:",l," identical:",l)') sfcmpt, sfid
      if(useexc) write(unitout, '("  Reading form V from ", a)') trim(efname)
      if(usescc) write(unitout, '("  compatible:",l," identical:",l)') efcmpt, efid

      ! Allocate measures
      if(fmeasure) then
        if(fcoup) then 
          if(allocated(vwdiffar)) deallocate(vwdiffar)
          allocate(vwdiffar(hamsize))
          vwdiffar=0.0d0
        else
          if(allocated(vwdiffrr)) deallocate(vwdiffrr)
          allocate(vwdiffrr(hamsize))
          vwdiffrr=0.0d0
        end if
      end if

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
          call getbsemat(trim(sfname), iqmt, ikkp,&
            & sccli_t(1:inou,1:jnou), check=.false., fcmpt=sfcmpt, fid=sfid)
        end if

        if(useexc) then 
          ! Read RR/RA part of exchange interaction v_{iuioik,jujojk}(qmt)
          call getbsemat(trim(efname), iqmt, ikkp,&
            & excli_t(1:inou,1:jnou), check=.false., fcmpt=efcmpt, fid=efid)
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
          call setint(ham(i1:i2,j1:j2),&
            & ofac(i1:i2), ofac(j1:j2),&
            & scc=sccli_t(1:inou,1:jnou), exc=excli_t(1:inou,1:jnou))
        else if(usescc) then
          call setint(ham(i1:i2,j1:j2),&
            & ofac(i1:i2), ofac(j1:j2),&
            & scc=sccli_t(1:inou,1:jnou))
        else if(useexc) then
          call setint(ham(i1:i2,j1:j2),&
            & ofac(i1:i2), ofac(j1:j2),&
            & exc=excli_t(1:inou,1:jnou))
        end if

        ! Maximum V-W for each row
        if(fmeasure) then 
          call makemeasure(i1, ham(i1:i2,j1:j2), fcoup)
        end if

      ! ikkp loop end
      end do

      !! RR only
      ! For blocks on the diagonal, add the KS transition
      ! energies to the diagonal of the block.
      if(.not. fcoup) then 
        call addkstransdiag(1, ham(:,:))
      end if

      ! Write lower triangular part in case it is explicitly needed
      do i1 = 1, hamsize
        ham(i1, i1) = cmplx(dble(ham(i1,i1)), 0.0d0, 8)
        do i2 = i1+1, hamsize
          ham(i2,i1) = conjg(ham(i1,i2))
        end do
      end do

      call timesec(ts1)
      write(unitout, '(" Matrix build.")')
      write(unitout, '("Timing (in seconds)	   :", f12.3)') ts1 - ts0

      contains

        subroutine setint(hamblock, oc1, oc2, scc, exc, w, v)
          complex(8), intent(out) :: hamblock(:,:)
          real(8), intent(in) :: oc1(:), oc2(:)
          complex(8), intent(in), optional :: scc(:,:)
          complex(8), intent(in), optional :: exc(:,:)
          complex(8), intent(out), optional :: w(:,:), v(:,:)
          ! local variables
          real (8) :: excfac 
          integer(4) :: i, j

          ! Set prefactor for exchange term in Hamiltonian
          if (input%xs%bse%xas) then
            if ((.not. input%groundstate%tevecsv) .and. (input%xs%bse%xasedge == 'K')) then
              excfac=2.0d0
            else
              excfac=1.0d0
            end if
          else
            if (.not. input%groundstate%tevecsv) excfac=2.0d0
            if (input%groundstate%tevecsv) excfac=1.0d0
          end if
          if(present(exc) .and. present(scc)) then 
            do j= 1, size(hamblock,2)
              do i= 1, size(hamblock,1)
                hamblock(i,j) = oc1(i)*oc2(j) * (excfac* exc(i,j) - scc(i,j))
              end do
            end do
          else if(present(scc)) then 
            do j= 1, size(hamblock,2)
              do i= 1, size(hamblock,1)
                hamblock(i,j) = -oc1(i)*oc2(j) * scc(i,j)
              end do
            end do
          else if(present(exc)) then 
            do j= 1, size(hamblock,2)
              do i= 1, size(hamblock,1)
                hamblock(i,j) = oc1(i)*oc2(j) *excfac * exc(i,j)
              end do
            end do
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
    ! !ROUTINE: setup_bse_tr_dist
    ! !INTERFACE:
    subroutine setup_bse_tr_dist(iqmt, binfo, smat, cmat, cpmat)
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
      character(*), parameter :: thisname = "setup_bse_tr_dist"
      type(dzmat) :: rrmat
      type(dzmat) :: ramat
      type(dzmat) :: auxmat

      logical :: fcheckpos, fsmat, fcmat, fcpmat

      real(8) :: ts0, ts1, t1, t0
      integer(4) :: i, j

      real(8), allocatable :: evals(:)

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

      write(unitout, '("Info(",a,"):&
        & Setting up distributed matrices for squared EVP")') trim(thisname)
      call timesec(ts0)

      !===========================================================!
      ! Getting main and coupling blocks of BSE hamiltonian       !
      !===========================================================!
      ! Get RR part of BSE Hamiltonian
      write(unitout, '("Info(",a,"):&
        & Setting up RR Block of orignial BSE")') trim(thisname)
      call timesec(t0)
      call new_dzmat(rrmat, hamsize, hamsize, binfo)
      call setup_bse_block_dist(rrmat, iqmt, .false., binfo)
      call timesec(t1)
      write(unitout, '("  Time needed",f12.3,"s")') t1-t0
      ! Get RA^{ti} part of BSE Hamiltonian
      write(unitout, '("Info(",a,"):&
        & Setting up RA^{ti} Block of orignial BSE")') trim(thisname)
      call timesec(t0)
      call new_dzmat(ramat, hamsize, hamsize, binfo)
      call setup_bse_block_dist(ramat, iqmt, .true., binfo)
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

    end subroutine setup_bse_tr_dist
    !EOC

    !BOP
    ! !ROUTINE: setup_bse_block_dist
    ! !INTERFACE:
    subroutine setup_bse_block_dist(ham, iqmt, fcoup, binfo)
    ! !USES:
      use modmpi
      use modscl
      use modinput, only: input
      use mod_constants, only: zzero, zone
      use modxs, only: unitout
      use modbse, only: nou_bse_max, kmap_bse_rg,&
                      & kousize, nkkp_bse, nk_bse, ofac,&
                      & scclifbasename, exclifbasename,&
                      & scclicfbasename,&
                      & infofbasename, &
                      & vwdiffrr, vwdiffar, hamsize
      use m_getunit
      use m_genfilname
      use m_putgetbsemat
      use m_writecmplxparts
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   integer(4) :: iqmt ! Index of momentum transfer Q
    !   logical :: fcoup   ! If true, builds RA instead of RR block of BSE matrix 
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
      logical, intent(in) :: fcoup
      type(blacsinfo), intent(in) :: binfo

      !! Local variables
      integer(4) :: ikkp, ik, jk, iknr, jknr
      integer(4) :: inou, jnou

      ! Arrays to read V and W from file
      complex(8), allocatable, dimension(:,:) :: excli_t, sccli_t
      ! V-W array
      complex(8), allocatable, dimension(:,:) :: vwdiff_t

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
      logical :: esfcmpt, esfid, ferror
      logical :: useexc, usescc

      logical :: fmeasure
      logical :: fchibarq

      ! Make coupling measures? 
      fmeasure = input%xs%bse%measure

      ! Use truncated Coulomb potential? 
      fchibarq = input%xs%bse%chibarq


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
            write(unitout, '("Info(setup_bse_block_dist):&
              & Using time reversal symmetry.")')
          end if

          ! Select form which files W and V are to be read
          if(.not. fcoup) then
            call genfilname(basename=scclifbasename, iqmt=iqmt, filnam=sfname)
          else
            call genfilname(basename=scclicfbasename, iqmt=iqmt, filnam=sfname)
          end if
          if(fchibarq) then 
            call genfilname(basename=exclifbasename, bsetype="-BAR", iqmt=iqmt, filnam=efname)
          else
            call genfilname(basename=exclifbasename, iqmt=iqmt, filnam=efname)
          end if

          sinfofname = trim(infofbasename)//'_'//trim(sfname)
          einfofname = trim(infofbasename)//'_'//trim(efname)

          if(usescc) write(unitout, '("  Reading info from ", a)')&
            & trim(sinfofname)
          if(useexc) write(unitout, '("  Reading info from ", a)')&
            & trim(einfofname)

          ! Check saved quantities for compatibility
          sfcmpt=.false.
          sfid=.false.
          if(usescc) call getbseinfo(trim(sinfofname), iqmt,&
            & fcmpt=sfcmpt, fid=sfid)
          efcmpt=.false.
          efid=.false.
          if(useexc) call getbseinfo(trim(einfofname), iqmt,&
            & fcmpt=efcmpt, fid=efid)

          if(usescc .and. useexc) then 
            esfcmpt = efcmpt .and. sfcmpt
            esfid = efid .and. sfid
            ferror = .not. (esfcmpt .and. esfid)
            if(ferror) then 
              write(*, '("Error(setup_bse_block_dist): Info files differ")')
              write(*, '("  efcmpt, efid, sfcmpt, sfcmpt")') efcmpt, efid, sfcmpt, sfcmpt
              call terminate
            end if
          end if

          if(usescc) write(unitout, '("  Reading W from ", a)') trim(sfname)
          if(usescc) write(unitout, '("  compatible:",l," identical:",l)') sfcmpt, sfid
          if(useexc) write(unitout, '("  Reading V from ", a)') trim(efname)
          if(usescc) write(unitout, '("  compatible:",l," identical:",l)') efcmpt, efid

          allocate(excli_t(nou_bse_max, nou_bse_max))
          allocate(sccli_t(nou_bse_max, nou_bse_max))

          ! Allocate measure arrays
          if(fmeasure) then
            allocate(vwdiff_t(nou_bse_max, nou_bse_max))
            if(fcoup) then 
              if(allocated(vwdiffar)) deallocate(vwdiffar)
              allocate(vwdiffar(hamsize))
              vwdiffar=0.0d0
            else
              if(allocated(vwdiffrr)) deallocate(vwdiffrr)
              allocate(vwdiffrr(hamsize))
              vwdiffrr=0.0d0
            end if
          end if

        ! is root
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
              call getbsemat(trim(sfname), iqmt, ikkp,&
                & sccli_t(1:inou,1:jnou), check=.false., fcmpt=sfcmpt, fid=sfid)
            end if

            ! Read in exchange interaction for ikkp
            if(useexc) then 
              ! Read RR/RA part of exchange interaction v_{iuioik,jujojk}(qmt)
              call getbsemat(trim(efname), iqmt, ikkp,&
                & excli_t(1:inou,1:jnou), check=.false., fcmpt=efcmpt, fid=efid)
            end if

            ! Make ik=jk blocks explicitly symmetric/hermitian
            if(iknr == jknr) then 
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

            ! Find maximal |V-W| per row
            if(fmeasure) then 

              vwdiff_t = 0.0d0
              if(usescc .and. useexc) then
                vwdiff_t(1:inou,1:jnou) = 2.0d0*excli_t(1:inou,1:jnou)&
                  & - sccli_t(1:inou,1:jnou)
              else if(usescc) then
                vwdiff_t(1:inou,1:jnou) = -sccli_t(1:inou,1:jnou)
              else if(useexc) then
                vwdiff_t(1:inou,1:jnou) = 2.0d0*excli_t(1:inou,1:jnou)
              end if

              call makemeasure(ii, vwdiff_t(1:inou,1:jnou), fcoup)

            end if

          ! is root
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
                  if(usescc) sbuff2(1:jb,1:ib) =&
                    & conjg(transpose(sccli_t(i:i+ib-1, j:j+jb-1)))
                  if(useexc) ebuff2(1:jb,1:ib) =&
                    & conjg(transpose(excli_t(i:i+ib-1, j:j+jb-1)))
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

            ! i while loop
            end do

            ! Next column block
            j = j + jb

          ! j while loop
          end do

        ! ikkp loop
        end do
#else
        write(*,*) "Error(setup_bse_block_dist): Scalapack needed."
        call terminate
#endif

      else

        call setup_bse_block(ham%za, iqmt, fcoup)

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
      complex(8) :: excfac 
      
      ! Set prefactor for exchange term in Hamiltonian
      if (input%xs%bse%xas) then
        if ((.not. input%groundstate%tevecsv) .and. (input%xs%bse%xasedge == 'K')) then
          excfac=ztwo
        else
          excfac=zone
        end if
      else
        if (.not. input%groundstate%tevecsv) excfac=ztwo
        if (input%groundstate%tevecsv) excfac=zone
      end if
      
      do c = 1, jb
        do r = 1, ib
          if(present(exc) .and. present(scc)) then
            ! Singlet case with exchange interaction
            hamblck(r, c) = occ1(r) * (excfac * exc(r, c) - scc(r, c)) * occ2(c)
          else if(present(scc)) then
            ! Triplet case without exchange interaction
            hamblck(r, c) = -occ1(r) * scc(r, c) * occ2(c)
          else if(present(exc)) then
            ! RPA
            hamblck(r, c) = occ1(r) * excfac * exc(r, c) * occ2(c)
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

    ! Finds the maximal component of |V-W| for each row
    subroutine makemeasure(ig, hamblock, fcoup)
      use modbse, only: vwdiffar, vwdiffrr
      integer(4), intent(in) :: ig
      complex(8), intent(in) :: hamblock(:,:)
      logical, intent(in) :: fcoup

      integer(4) :: i, m, n
      real(8) :: maxdiffrow

      m = size(hamblock,1)
      n = size(hamblock,2)

      ! Get maximal |V-W| for each row 
      do i = 1, m
        maxdiffrow = maxval(abs(hamblock(i,:)))
        if(fcoup) then 
          if(vwdiffar(ig+i-1) < maxdiffrow) then
            vwdiffar(ig+i-1) = maxdiffrow
          end if
        else
          if(vwdiffrr(ig+i-1) < maxdiffrow) then
            vwdiffrr(ig+i-1) = maxdiffrow
          end if
        end if
      end do
    end subroutine makemeasure


end module m_setup_bse
