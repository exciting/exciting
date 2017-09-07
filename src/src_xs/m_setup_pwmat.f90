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
                 & qvkloff, ngq, ivgigq, xsgnt, nqmt
  use mod_Gvector, only: ivg
  use mod_Gkvector, only: gkmax
  use mod_xsgrids
  use m_b_ematqk
  use m_genfilname
  use m_writecmplxparts
  use m_xsgauntgen
  use m_findgntn0
  use mod_variation, only: ematqk_sv

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
    !   complex(8) :: pwmat(hamsize)  ! Plane wave matrix elemens
    ! 
    ! !DESCRIPTION:
    !   The routine generates the plane wave matrix elements 
    !
    !   \begin{equation*}
    !   \tilde{M}_{\alpha}(\vec{G},\vec{q}_\text{mt}) = 
    !    \sqrt{
    !    f_{o_\alpha, \vec{k}_{\alpha+\vec{q}_\text{mt}/2}}
    !   -f_{u_\alpha, \vec{k}_{\alpha-\vec{q}_\text{mt}/2}}
    !    }
    !    \langle u_\alpha \vec{k}_{\alpha-\vec{q}_\text{mt}/2}
    !    | e^{-i(\vec{G}_\text{mt}+\vec{q}_\text{mt})r} |
    !    o_\alpha \vec{k}_{\alpha+\vec{q}_\text{mt}} \rangle
    !   \end{equation*}
    !
    !   Alpha is the combined index used in the BSE Hamiltonian.
    !
    ! !REVISION HISTORY:
    !   Created. (Aurich)
    !EOP
    !BOC
      integer(4), intent(in) :: iqmt, igqmt
      complex(8), intent(out) :: pwmat(hamsize)

      integer(4) :: io, iu, ik, iknr, ikmnr, ik1, ik2
      integer(4) :: ino, inu, ioabs1, iuabs1, ioabs2, iuabs2 
      integer(4) :: a1, numgq
      integer(4) :: i,j, igq

      type(bcbs) :: ematbc
      character(256) :: fileext0_save, fileext_save
      complex(8), allocatable :: muo(:,:,:), muog(:,:,:), mou_(:,:,:)
      integer(4), allocatable, target :: ikm2ikp_dummy(:,:)
      real(8), parameter :: epslat = 1.0d-8
      logical :: fsamekm, fsamekp

      !------------------------------------!
      ! Setup grids k, k-qmt/2 and k+qmt/2
      !------------------------------------!
      call xsgrids_init(totalqlmt(1:3, iqmt), gkmax) 

      ! Save the k,G arrays of the k-qmt/2 grid to modxs::vkl0 etc
      call init1offs(k_kqmtm%kqmtset%vkloff)
      call xssave0
      ! Save the k,G arrays of the k+qmt/2 grid to the default locations
      call init1offs(k_kqmtp%kqmtset%vkloff)

      ! Check whether k+-qmt/2 grids are identical to k grid
      if(all(abs(k_kqmtp%kqmtset%vkloff-k_kqmtp%kset%vkloff) < epslat)) then 
        write(unitout, '("Info(setup_pwmat):&
          & k+qmt/2-grid is identical for to iqmt=1 grid, iqmt=",i3)') iqmt
        fsamekp=.true.
      else
        fsamekp=.false.
      end if
      if(all(abs(k_kqmtm%kqmtset%vkloff-k_kqmtm%kset%vkloff) < epslat)) then 
        write(unitout, '("Info(setup_pwmat): k-qmt/2-grid is identical for to iqmt=1 grid, iqmt=",i3)') iqmt
        fsamekm=.true.
      else
        fsamekm=.false.
      end if

      ! Make link map connecting ik-qmt/2 to ik+qmt/2
      ! ( only the iqmt component will be used ...)
      if(allocated(ikm2ikp_dummy)) deallocate(ikm2ikp_dummy)
      allocate(ikm2ikp_dummy(k_kqmtp%kset%nkptnr, nqmt))
      ikm2ikp_dummy=0
      ikm2ikp_dummy(:, iqmt) = ikm2ikp(:)
      !------------------------------------!

      !----------------------------------!
      ! Specify form which files to read
      !----------------------------------!
      ! Save file extension
      fileext_save = filext
      fileext0_save = filext0

      usefilext0 = .true.
      if(fsamekm) then 
        ! Set EVECFV_QMT001.OUT as bra state file
        iqmt0 = iqmtgamma
        call genfilname(iqmt=iqmt0, setfilext=.true.)
      else
        ! Set EVECFV_QMTXYZ_m.OUT as bra state file
        iqmt0 = iqmt
        call genfilname(iqmt=iqmt0, auxtype="m", fileext=filext0)
      end if

      ! Set EVECFV_QMTXYZ.OUT as ket state file (k+qmt/2 grid)
      iqmt1 = iqmt
      ! Set to EVECFV_QMT001.OUT if it is the same grid
      if(fsamekp) iqmt1 = iqmtgamma
      call genfilname(iqmt=iqmt1, setfilext=.true.)
      !----------------------------------!

      ! Use normal computation of the plane wave matrix elements M
      emat_ccket=.false.

      ! Generate gaunt coefficients used in the construction of 
      ! the plane wave matrix elements in ematqk.
      call xsgauntgen(max(input%groundstate%lmaxapw, lolmax),&
        & input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))
      ! Find indices for non-zero gaunt coefficients
      call findgntn0(max(input%xs%lmaxapwwf, lolmax),&
        & max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)

      ! Allocate needed arrays (evecfv etc.)
      call ematqalloc

      ! Set up ikmapikq to link (jkm,iqmt) to (jkp)
      ikmapikq_ptr => ikm2ikp_dummy

      ! Set vkl0_ptr k-qmt/2-grid, vkl1_ptr, ... to k+qmt/2-grid
      ! Note: Needs to be called after ematqalloc
      call setptr01()

      ! Calculate radial integrals used in the construction 
      ! of the plane wave matrix elements for exponent (G+qmt)
      call ematrad(iqmt)

      numgq = ngq(iqmt)

      allocate(muo(nu_bse_max, no_bse_max, numgq))
      allocate(muog(nu_bse_max, no_bse_max, nk_bse))
      muo=zzero
      muog=zzero

      if(input%xs%bse%distribute) then 
        ik1=firstofset(mpiglobal%rank, nk_bse)
        ik2=lastofset(mpiglobal%rank, nk_bse)
      else
        ik1=1
        ik2=nk_bse
      end if

      do ik = ik1, ik2

        !------------------------------------------------------------------!
        ! Calculate M_{iu io ikm}(G, qmt) = <iu ikm|e^{-i(G+qmt)r}|io ikp> !
        ! where ikp = ik+qmt/2 and ikm = ik-qmt/2                          !
        !------------------------------------------------------------------!

        ! Get non reduced global k index 
        iknr = kmap_bse_rg(ik)

        ! Get selected occupation limits for that k point
        iuabs1 = koulims(1,iknr)
        iuabs2 = koulims(2,iknr)
        ioabs1 = koulims(3,iknr)
        ioabs2 = koulims(4,iknr)
        inu = iuabs2 - iuabs1 + 1
        ino = ioabs2 - ioabs1 + 1

        ! Get index of associated ik-qmt/2
        ikmnr = k_kqmtm%ik2ikqmt(iknr)

        ! Set ranges for the calculation of the plane wave matrix elements 
        ematbc%n1=inu
        ematbc%il1=iuabs1
        ematbc%iu1=iuabs2
        ematbc%n2=ino
        ematbc%il2=ioabs1
        ematbc%iu2=ioabs2

        ! Calculate M_{uo,G} at fixed (k, q)
        if(input%xs%bse%xas) then
          call xasgauntgen (input%xs%lmaxemat, Max(input%groundstate%lmaxapw, lolmax))
          call b_ematqk_core(iqmt, ikmnr, muo(1:inu,1:ino,:),ematbc,'uo')
        else
          if (.not. (input%groundstate%tevecsv)) then
            call b_ematqk(iqmt, ikmnr, muo(1:inu,1:ino,:), ematbc)
          else
            call ematqk_sv(iqmt, ikmnr, muo(1:inu,1:ino,:), ematbc)
          end if
        end if

        ! Save only selected G=G_mt
        muog(1:inu, 1:ino, ik) = muo(1:inu,1:ino,igqmt)
      end do

      ! Gather moug 
      if(input%xs%bse%distribute) then 
        call mpi_allgatherv_ifc(set=nk_bse, rlen=nu_bse_max*no_bse_max, zbuf=muog,&
          & inplace=.true., comm=mpiglobal)
      end if

      call ematqdealloc
      deallocate(muo)

      ! Return plane wave matrix elements weighted with 
      ! occupation factors in combined index notation.
      do a1 = 1, hamsize
        ! Relative indices
        iu = smap_rel(1,a1)
        io = smap_rel(2,a1)
        ik = smap_rel(3,a1)
        pwmat(a1) = ofac(a1)*muog(iu, io, ik)
      end do

      deallocate(muog)
      deallocate(ikm2ikp_dummy)
      call xsgrids_finalize()

      ! Restore file extension
      filext=fileext_save
      filext0=fileext0_save

    end subroutine setup_pwmat
    !EOC

    !BOP
    ! !ROUTINE: setup_pwmat_dist
    ! !INTERFACE:
    subroutine setup_pwmat_dist(dpwmat, iqmt, igqmt, binfo)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   integer(4) :: iqmt            ! Index of momentum transfer
    !   integer(4) :: igqmt           ! Index of G+qmt
    !   type(blacsinfo) :: binfo      ! Info type for BLACS grid
    ! Out:
    !   type(dzmat) :: dpwma          ! 2D block cyclic distributed plane wave
    !                                 ! matrix elements 
    ! 
    ! !DESCRIPTION:
    !   The routine generates the plane wave matrix elements 
    !   \begin{equation*}
    !     \tilde{M}_{\alpha}(G,qmt) =
    !     \sqrt{f_{v_\alpha, \vec{k}_\alpha}-f_{c_\alpha, \vec{k}_\alpha+\vec{q}_\text{mt}}}
    !     \langle v_\alpha \vec{k}_\alpha | e^{-i(G+qmt)r} | c_\alpha \vec{k}_\alpha+qmt \rangle
    !   \end{equation*}
    !
    !   Alpha is the combined index used in the BSE Hamiltonian.
    !
    ! !REVISION HISTORY:
    !   Created. (Aurich)
    !EOP
    !BOC
      integer(4), intent(in) :: iqmt, igqmt
      type(blacsinfo), intent(in) :: binfo
      type(dzmat), intent(inout) :: dpwmat

      complex(8) :: pwmat(hamsize,1)

      call setup_pwmat(pwmat(:,1), iqmt, igqmt)

      if(binfo%isactive) then 
        call dzmat_copy_global2local(pwmat, dpwmat, binfo)
      end if
      
    end subroutine setup_pwmat_dist
    !EOC

end module m_setup_pwmat
