module m_setup_pwmat
  use modmpi
  use modscl
  use modinput
  use mod_constants
  use modbse
  use mod_kpoint, only: vkl
  use mod_misc, only: filext
  use mod_APW_LO, only: lolmax
  use modxs, only: unitout, bcbs,&
                 & vkl0, usefilext0, filext0, iqmt0, iqmt1, iqmtgamma,&
                 & qvkloff, ngq, ivgigq, xsgnt
  use mod_Gvector, only: ivg
  use m_b_ematqk
  use m_genfilname
  use m_writecmplxparts
  use m_xsgauntgen
  use m_findgntn0

  implicit none

  contains

    !BOP
    ! !ROUTINE: setup_pwmat
    ! !INTERFACE:
    subroutine setup_pwmat(pwmat, iqmt, igqmt)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   integer(4) :: iqmt            ! Index of momentum transfer
    !   integer(4) :: igqmt           ! Index of G+qmt
    ! Out:
    !   complex(8) :: pwmat(hamsize)  ! Plane wave matrix elemens <ok|e^{-i(G+qmt)r}|uk+qmt>
    ! 
    ! !DESCRIPTION:
    !   The routine generates the plane wave matrix elements 
    !   $\tilde{M}_{\alpha}(G,qmt) = \sqrt{f_{v_\alpha, \vec{k}_\alpha}-f_{c_\alpha, \vec{k}_\alpha+\vec{q}_\text{mt}}}
    !    \langle v_\alpha \vec{k}_\alpha | e^{-i(G+qmt)r} | c_\alpha \vec{k}_\alpha+qmt \rangle
    !
    !   Alpha is the combined index used in the BSE Hamiltonian.
    !
    ! !REVISION HISTORY:
    !   Created. (Aurich)
    !EOP
    !BOC
      integer(4), intent(in) :: iqmt, igqmt
      complex(8), intent(out) :: pwmat(hamsize)

      real(8) :: t1, t0
      integer(4) :: io, iu, ik, iknr
      integer(4) :: ino, inu, ioabs1, iuabs1, ioabs2, iuabs2 
      integer(4) :: a1, numgq

      type(bcbs) :: ematbc
      character(256) :: fileext0_save, fileext_save
      complex(8), allocatable :: mou(:,:,:), moug(:,:,:)

      write(unitout, '("  Building Pwmat.")')
      call timesec(t0)

      ! Save file extension
      fileext_save = filext
      fileext0_save = filext0

      ! Set EVECFV_QMT001.OUT as bra state file
      usefilext0 = .true.
      iqmt0 = iqmtgamma
      call genfilname(iqmt=iqmt0, setfilext=.true.)
      filext0 = filext
      !write(*,*) "(setup_pwmat): filext0 =", trim(filext)

      if(iqmt /= 1) then 
        ! Set EVECFV_QMTXYZ.OUT as ket state file
        iqmt1 = iqmt
        call genfilname(iqmt=iqmt1, setfilext=.true.)
        !write(*,*) "(setup_pwmat): filext =", trim(filext)
      else
        ! Set EVECFV_QMT001.OUT as ket state file
        iqmt1 = iqmtgamma
        call genfilname(iqmt=iqmt1, setfilext=.true.)
        !write(*,*) "(setup_pwmat): filext =", trim(filext)
      end if

      ! Set vkl0 to k-grid
      ! Note: This needs init2 to be called form task 441 beforehand
      !write(*,*) "(setup_pwmat): iqmt=", iqmt
      !write(*,*) "(setup_pwmat): qvkloff(:,1) =", qvkloff(:,1)
      call init1offs(qvkloff(:,1))
      call xssave0
      ! Set vkl to k+qmt-grid
      !write(*,*) "(setup_pwmat): qvkloff(:,iqmt) =", qvkloff(:,iqmt)
      call init1offs(qvkloff(:,iqmt))

      ! Generate gaunt coefficients used in the construction of 
      ! the plane wave matrix elements in ematqk.
      call xsgauntgen(max(input%groundstate%lmaxapw, lolmax),&
        & input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))
      ! Find indices for non-zero gaunt coefficients
      call findgntn0(max(input%xs%lmaxapwwf, lolmax),&
        & max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)

      call ematqalloc

      ! Calculate radial integrals used in the construction 
      ! of the plane wave matrix elements for exponent (G+qmt)
      !write(*,*) "(setup_pwmat): calculating radial integrals"
      call ematrad(iqmt)

      numgq = ngq(iqmt)
      !write(*,*) "(setup_pwmat): numgq =", numgq

      allocate(mou(no_bse_max, nu_bse_max, numgq))
      allocate(moug(no_bse_max, nu_bse_max, nk_bse))
      mou=zzero
      moug=zzero

      ! Read in all needed momentum matrix elements
      do ik = 1, nk_bse

        !----------------------------------------------------------------!
        ! Calculate M_{io iu ik}(G, qmt) = <io ik|e^{-i(qmt+G)r}|iu ikp> !
        ! where ikp = ik+qmt                                             !
        !----------------------------------------------------------------!
        iknr = kmap_bse_rg(ik)

        !write(*,*) "(setup_pwmat): calculating pwmat for ik =", iknr

        iuabs1 = koulims(1,iknr)
        iuabs2 = koulims(2,iknr)
        ioabs1 = koulims(3,iknr)
        ioabs2 = koulims(4,iknr)
        inu = iuabs2 - iuabs1 + 1
        ino = ioabs2 - ioabs1 + 1

        ematbc%n1=ino
        ematbc%il1=ioabs1
        ematbc%iu1=ioabs2
        ematbc%n2=inu
        ematbc%il2=iuabs1
        ematbc%iu2=iuabs2

        ! Calculate M_{o1o2,G} at fixed (k, q)
        call b_ematqk(iqmt, iknr, mou(1:ino,1:inu,:), ematbc)

        !write(*,*) "passed ematqk"

        ! Save only selected G
        moug(1:ino, 1:inu, ik) = mou(1:ino,1:inu,igqmt)

      end do

      call ematqdealloc

      do a1 = 1, hamsize
        ! Relative indices
        iu = smap_rel(1,a1)
        io = smap_rel(2,a1)
        ik = smap_rel(3,a1)
        pwmat(a1) = ofac(a1)*moug(io, iu, ik)
      end do

      ! Restore file extension
      filext=fileext_save
      filext0=fileext0_save

      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
    end subroutine setup_pwmat
    !EOC


end module m_setup_pwmat
