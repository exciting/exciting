module m_setup_bse
  use modmpi
  use modscl
  use modinput, only: input
  use mod_constants, only: zzero, zone
  use mod_eigenvalue_occupancy, only: evalsv
  use modxs, only: evalsv0
  use modbse
  use m_getunit
  use m_genfilname
  use m_putgetbsemat
  use m_writecmplxparts
      
  implicit none

  private

  public :: setup_bse, setup_distributed_bse

  contains

    !BOP
    ! !ROUTINE: setup_bse
    ! !INTERFACE:
    subroutine setup_bse(ham, iqmt, fcoup, fti)
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
      if(mpiglobal%rank == 0) then 
        if(.not. fcoup) then 
          write(unitout, '("Info(setup_bse): Setting up RR part of hamiltonian")')
        else
          write(unitout, '("Info(setup_bse): Setting up RA part of hamiltonian")')
          if(fti) then 
            write(unitout, '("Info(setup_bse):&
              & Using time inverted anti-resonant basis")')
          end if
        end if
      end if

      usescc = .false.
      useexc = .false.
      select case(trim(input%xs%bse%bsetype))
        case('IP')
          usescc = .false.
          useexc = .false.
          if(mpiglobal%rank==0) then
            write(unitout, '("Info(setup_bse): Using IP")')
          end if
        case('singlet')
          usescc = .true.
          useexc = .true.
          if(mpiglobal%rank==0) then
            write(unitout, '("Info(setup_bse): Using singlet")')
          end if
        case('triplet')
          usescc = .true.
          useexc = .false.
          if(mpiglobal%rank==0) then
            write(unitout, '("Info(setup_bse): Using triplet")')
          end if
        case('RPA')
          usescc = .false.
          useexc = .true.
          if(mpiglobal%rank==0) then
            write(unitout, '("Info(setup_bse): Using RPA")')
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

      if(mpiglobal%rank == 0) then 
        if(usescc) write(unitout, '("  Reading form info from ", a)')&
          & trim(sinfofname)
        if(useexc) write(unitout, '("  Reading form info from ", a)')&
          & trim(einfofname)
      end if

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
        if(efcmpt /= sfcmpt .or. efid /= sfid) then 
          write(*, '("Error(setup_bse): Info files differ")')
          write(*, '("  efcmpt, efid, sfcmpt, sfcmpt")') efcmpt, efid, sfcmpt, sfcmpt
          call terminate
        end if
      end if

      if(mpiglobal%rank == 0) then 
        if(usescc) write(unitout, '("  Reading form W from ", a)') trim(sfname)
        if(usescc) write(unitout, '("  compatible:",l," identical:",l)') sfcmpt, sfid
        if(useexc) write(unitout, '("  Reading form V from ", a)') trim(efname)
        if(usescc) write(unitout, '("  compatible:",l," identical:",l)') efcmpt, efid
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
          call b_getbsemat(trim(sfname), iqmt, ikkp,&
            & sccli_t(1:inou,1:jnou), check=.false., fcmpt=sfcmpt, fid=sfid)
        end if

        if(useexc) then 
          ! Read RR/RA part of exchange interaction v_{iuioik,jujojk}(qmt)
          call b_getbsemat(trim(efname), iqmt, ikkp,&
            & excli_t(1:inou,1:jnou), check=.false., fcmpt=efcmpt, fid=efid)
        end if

        if(fwp) then 
          if(ikkp == 1 .and. fcoup == .false.) then 
            if(usescc) call writecmplxparts('ikkp1_read_Wrr',dble(sccli_t(1:inou,1:jnou)), aimag(sccli_t(1:inou,1:jnou)))
            if(useexc) call writecmplxparts('ikkp1_read_Vrr',dble(excli_t(1:inou,1:jnou)), aimag(excli_t(1:inou,1:jnou)))
          end if

          if(ikkp == 1 .and. fcoup == .true. .and. fti == .false.) then 
            if(usescc) call writecmplxparts('ikkp1_read_Wra',dble(sccli_t(1:inou,1:jnou)), aimag(sccli_t(1:inou,1:jnou)))
            if(useexc) call writecmplxparts('ikkp1_read_Vra',dble(excli_t(1:inou,1:jnou)), aimag(excli_t(1:inou,1:jnou)))
          end if

          if(ikkp == 1 .and. fcoup == .true. .and. fti == .true.) then 
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
            write(unitout, '("Info(setup_bse): Writing RR: H, W, V and E to file")')
          else
            if(fti) then 
              write(unitout, '("Info(setup_bse): Writing RA^ti: H, W, V to file")')
            else
              write(unitout, '("Info(setup_bse): Writing RA: H, W, V to file")')
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
            do j= 1, size(hamblock,2)
              do i= 1, size(hamblock,1)
                hamblock(i,j) = oc1(i)*oc2(j) * (2.0d0 * exc(i,j) - scc(i,j))
                !write(*,*) "i,j", i, j
                !write(*,*) "oc1,oc2", oc1(i), oc1(j)
                !write(*,*) "exc,scc", exc(i,j), scc(i,j)
                !write(*,*) "hamblock", hamblock(i,j)
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
                hamblock(i,j) = oc1(i)*oc2(j) * 2.0d0 * exc(i,j)
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

    end subroutine setup_bse
    !EOC

    !BOP
    ! !ROUTINE: setup_distributed_bse
    ! !INTERFACE:
    subroutine setup_distributed_bse(ham, iqmt, fcoup, fti, binfo)
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
              write(unitout, '("Info(setup_distributed_bse): Using IP")')
            end if
          case('singlet')
            usescc = .true.
            useexc = .true.
            if(mpiglobal%rank==0) then
              write(unitout, '("Info(setup_distributed_bse): Using singlet")')
            end if
          case('triplet')
            usescc = .true.
            useexc = .false.
            if(mpiglobal%rank==0) then
              write(unitout, '("Info(setup_distributed_bse): Using triplet")')
            end if
          case('RPA')
            usescc = .false.
            useexc = .true.
            if(mpiglobal%rank==0) then
              write(unitout, '("Info(setup_distributed_bse): Using RPA")')
            end if
        end select

        ! Check whether the saved data is usable and 
        ! allocate read in work arrays
        if(binfo%isroot) then 

          if(.not. fcoup) then 
            write(unitout, '("Info(setup_distributed_bse):&
              & Setting up RR part of hamiltonian")')
          else
            write(unitout, '("Info(setup_distributed_bse):&
              & Setting up RA part of hamiltonian")')
            if(fti) then 
              write(unitout, '("Info(setup_distributed_bse):&
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
            if(efcmpt /= sfcmpt .or. efid /= sfid) then 
              write(*, '("Error(setup_distributed_bse): Info files differ")')
              write(*, '("  efcmpt, efid, sfcmpt, sfcmpt")') efcmpt, efid, sfcmpt, sfcmpt
              call terminate
            end if
          end if

          if(usescc) write(unitout, '("  Reading form W from ", a)') trim(sfname)
          if(useexc) write(unitout, '("  Reading form V from ", a)') trim(efname)

          allocate(excli_t(nou_bse_max, nou_bse_max))
          allocate(sccli_t(nou_bse_max, nou_bse_max))

        end if

        ! Block sizes of global matrix
        iblck = ham%mblck
        jblck = ham%nblck 

!if(binfo%isroot) then 
!  write(*,*) "iblck, jblck", iblck, jblck
!end if

        ! Context
        context = ham%context
        if(context /= binfo%context) then 
          write(*, '("Error(setup_distributed_bse): Context mismatch:", 2i4)')&
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

!if(binfo%isroot) then 
!  write(*,*) "ikkp", ikkp
!end if
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

!if(binfo%isroot) then 
!  write(*,*) "inou x jnou", inou, jnou
!end if
!if(binfo%isroot) then 
!  write(*,*) "ii, jj", ii, jj
!end if

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
              if(ikkp == 1 .and. fcoup == .false.) then 
                if(usescc) then 
                  call writecmplxparts('ikkp1_read_Wrr',dble(sccli_t(1:inou,1:jnou)), aimag(sccli_t(1:inou,1:jnou)))
                end if
                if(useexc) then 
                  call writecmplxparts('ikkp1_read_Vrr',dble(excli_t(1:inou,1:jnou)), aimag(excli_t(1:inou,1:jnou)))
                end if
              end if

              if(ikkp == 1 .and. fcoup == .true. .and. fti == .false.) then 
                if(usescc) then 
                  call writecmplxparts('ikkp1_read_Wra',dble(sccli_t(1:inou,1:jnou)), aimag(sccli_t(1:inou,1:jnou)))
                end if
                if(useexc) then 
                  call writecmplxparts('ikkp1_read_Vra',dble(excli_t(1:inou,1:jnou)), aimag(excli_t(1:inou,1:jnou)))
                end if
              end if

              if(ikkp == 1 .and. fcoup == .true. .and. fti == .true.) then 
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

!if(binfo%isroot) then 
!  write(*,*) "j, jg", j, jg
!  write(*,*) "jb", jb
!end if
!if(binfo%isroot) then 
!  write(*,*) "pc", pc
!  write(*,*) "jl", jl
!end if

            ! Get corresponding quantities in the lower 
            ! triangular part of the matrix.
            if(iknr /= jknr) then 
              ! Get process grid row coordinate
              pr2 = indxg2p( jg, iblck, binfo%myprow, 0, binfo%nprows)
              ! Get row position of jb*ib block in local matrix
              il2 = indxg2l( jg, iblck, pr2, 0, binfo%nprows)

!if(binfo%isroot) then 
!  write(*,*) "pr2", pr2
!  write(*,*) "il2", il2
!end if
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

!if(binfo%isroot) then 
!  write(*,*) "i, ig", i, ig
!  write(*,*) "ib", ib
!end if
!if(binfo%isroot) then 
!  write(*,*) "pr", pr
!  write(*,*) "il", il
!end if

              ! Get corresponding quantities in the lower 
              ! triangular part of the matrix.
              if(iknr /= jknr) then 
                ! Get process grid column coordinate
                pc2 = indxg2p( ig, jblck, binfo%mypcol, 0, binfo%npcols)
                ! Get column position of jb*ib block in local matrix
                jl2 = indxg2l( ig, jblck, pc2, 0, binfo%npcols)

!if(binfo%isroot) then 
!  write(*,*) "pc2", pc2
!  write(*,*) "jl2", jl2
!end if
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

!write(*,'("(",i2,",",i2,"): Building upper")') binfo%myprow, binfo%mypcol

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

!write(*,'("(",i2,",",i2,"): Sending upper to: ("i2,",",i2,")")')&
! binfo%myprow, binfo%mypcol, pr, pc

                  if(usescc) call zgesd2d(context, ib, jb, sbuff, iblck, pr, pc)
                  if(useexc) call zgesd2d(context, ib, jb, ebuff, iblck, pr, pc)

                end if

                ! Lower triangular part
                if(iknr /= jknr) then 

                  if( pr2 == 0 .and. pc2 == 0) then 

!write(*,'("(",i2,",",i2,"): Building lower")') binfo%myprow, binfo%mypcol

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

!write(*,'("(",i2,",",i2,"): Sending lower to: ("i2,",",i2,")")')&
! binfo%myprow, binfo%mypcol, pr, pc

                    if(usescc) call zgesd2d(context, jb, ib, sbuff2, jblck, pr2, pc2)
                    if(useexc) call zgesd2d(context, jb, ib, ebuff2, jblck, pr2, pc2)

                  end if

                end if

              ! All others only receive
              else
               
                ! Upper triangular part
                if(pr == binfo%myprow  .and. pc == binfo%mypcol) then

!write(*,'("(",i2,",",i2,"): Receiving upper form: ("i2,",",i2,")")')&
! binfo%myprow, binfo%mypcol, 0 , 0

                  ! Receive block
                  if(usescc) call zgerv2d(context, ib, jb, sbuff, iblck, 0, 0)
                  if(useexc) call zgerv2d(context, ib, jb, ebuff, iblck, 0, 0)

!write(*,'("(",i2,",",i2,"): Building upper")') binfo%myprow, binfo%mypcol

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

!write(*,'("(",i2,",",i2,"): Receiving lower form: ("i2,",",i2,")")')&
! binfo%myprow, binfo%mypcol, 0 , 0

                    ! Receive block
                    if(usescc) call zgerv2d(context, jb, ib, sbuff2, jblck, 0, 0)
                    if(useexc) call zgerv2d(context, jb, ib, ebuff2, jblck, 0, 0)

!write(*,'("(",i2,",",i2,"): Building lower")') binfo%myprow, binfo%mypcol

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
        write(*,*) "Error(setup_distributed_bse): Scalapack needed."
        call terminate
#endif

      else

        call setup_bse(ham%za, iqmt, fcoup, fti)

      end if

    end subroutine setup_distributed_bse
    !EOC

    !BOP
    ! !ROUTINE: buildham
    ! !INTERFACE:
    subroutine buildham(fc, hamblck, ig, jg, ib, jb, occ1, occ2, scc, exc)
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
      
      do c = 1, jb
        do r = 1, ib
          if(present(exc) .and. present(scc)) then
            ! Singlet case with exchange interaction
            hamblck(r, c) = occ1(r) * (ztwo * exc(r, c) - scc(r, c)) * occ2(c)
          else if(present(scc)) then
            ! Triplet case without exchange interaction
            hamblck(r, c) = -occ1(r) * scc(r, c) * occ2(c)
          else if(present(exc)) then
            ! RPA
            hamblck(r, c) = occ1(r) * ztwo * exc(r, c) * occ2(c)
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
