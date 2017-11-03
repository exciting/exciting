! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
module m_b_getgrst
  contains
    
    !BOP
    ! !ROUTINE: b_getevecfv0
    ! !INTERFACE:
    subroutine b_getevecfv0(vpl, vgpl, evecfvt)
    ! !USES:
      use mod_kpoint, only: vkl_ptr
      use mod_Gkvector, only: ngkmax_ptr, vgkl_ptr, ngk_ptr
      use mod_eigensystem, only: nmatmax_ptr
      use mod_misc, only: filext
      use modxs, only: filext0
      use mod_spin, only: nspnfv
      use mod_eigenvalue_occupancy, only: nstfv
      use mod_ematptr
      use mod_constants, only: zzero
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    !   real(8) :: vpl(3)          ! k-point vector in lattice coordinates
    !   real(8) :: vgpl(3, ngkmax) ! G+k-vectors in lattice coordinates
    ! OUT:
    !   complex(8) :: evecfvt(nmatmax, nstfv, nspnfv) ! Eigenvectors at that k-point
    !
    ! !DESCRIPTION:
    !   This routine is a wrapper for {\tt getevecfv} that changes the $k$ and $G+k$ 
    !   quantities in {\tt mod_kpoint} and {\tt mod_Gkvector} to the corresponding
    !   quantities saved in {\tt modxs} (nmatmax0, vkl0,  ngk0, etc.), changes
    !   the file extension in {\tt mod_misc} accordingly, reads in the Eigenvector
    !   and finally restores the original state.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !
    !EOP
    !BOC
      implicit none

      ! Arguments
      real(8), intent(in) :: vpl(3)
      real(8), intent(in) :: vgpl(3,ngkmax0_ptr,nspnfv)
      complex(8), intent(out) :: evecfvt(nmatmax0_ptr,nstfv,nspnfv)

      ! Local variables
      integer, pointer :: nmatmax_ptr_save, ngkmax_ptr_save
      integer, pointer :: ngk_ptr_save(:,:)
      real(8), pointer :: vkl_ptr_save(:,:), vgkl_ptr_save(:,:,:,:)
      character(256) :: filext_save

      nullify(nmatmax_ptr_save)
      nullify(ngkmax_ptr_save)
      nullify(ngk_ptr_save)
      nullify(vkl_ptr_save)
      nullify(vgkl_ptr_save)

      ! Backup pointers to default locations
      nmatmax_ptr_save => nmatmax_ptr
      ngkmax_ptr_save => ngkmax_ptr
      ngk_ptr_save => ngk_ptr
      vkl_ptr_save => vkl_ptr
      vgkl_ptr_save => vgkl_ptr

      filext_save = filext

      ! Set default pointers to bra state quantities
      nmatmax_ptr => nmatmax0_ptr
      ngkmax_ptr => ngkmax0_ptr
      ngk_ptr => ngk0_ptr
      vkl_ptr => vkl0_ptr
      vgkl_ptr => vgkl0_ptr

      filext = filext0
      ! Call to getevecfv with changed default (G+)k-set pointers / matrix size
      call getevecfv(vpl, vgpl, evecfvt)

      ! Restore default pointers
      nmatmax_ptr => nmatmax_ptr_save
      ngkmax_ptr => ngkmax_ptr_save
      ngk_ptr => ngk_ptr_save
      vkl_ptr => vkl_ptr_save
      vgkl_ptr => vgkl_ptr_save

      filext = filext_save

    end subroutine b_getevecfv0
    !EOC
    
    !BOP
    ! !ROUTINE: b_getevecfv1
    ! !INTERFACE:
    subroutine b_getevecfv1(vpl, vgpl, evecfvt)
    ! !USES:
      use mod_kpoint, only: vkl_ptr
      use mod_Gkvector, only: ngkmax_ptr, vgkl_ptr, ngk_ptr
      use mod_eigensystem, only: nmatmax_ptr
      use mod_spin, only: nspnfv
      use mod_eigenvalue_occupancy, only: nstfv
      use mod_ematptr
      use mod_constants, only: zzero
      use mod_misc, only: filext
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    !   real(8) :: vpl(3)          ! k-point vector in lattice coordinates
    !   real(8) :: vgpl(3, ngkmax) ! G+k-vectors in lattice coordinates
    ! OUT:
    !   complex(8) :: evecfvt(nmatmax, nstfv, nspnfv) ! Eigenvectors at that k-point
    !
    ! !DESCRIPTION:
    !   This routine is a wrapper for {\tt getevecfv} that changes the $k$ and $G+k$ 
    !   quantities in {\tt mod_kpoint} and {\tt mod_Gkvector} to the corresponding
    !   quantities saved in {\tt modxs} (nmatmax0, vkl0,  ngk0, etc.), changes
    !   the file extension in {\tt mod_misc} accordingly, reads in the Eigenvector
    !   and finally restores the original state.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !
    !EOP
    !BOC
      implicit none
     ! Arguments
      real(8), intent(in) :: vpl(3)
      real(8), intent(in) :: vgpl(3,ngkmax1_ptr,nspnfv)
      complex(8), intent(out) :: evecfvt(nmatmax1_ptr,nstfv,nspnfv)
    
      ! Local variables
      integer, pointer :: nmatmax_ptr_save, ngkmax_ptr_save
      integer, pointer :: ngk_ptr_save(:,:)
      real(8), pointer :: vkl_ptr_save(:,:), vgkl_ptr_save(:,:,:,:)
      nullify(nmatmax_ptr_save)
      nullify(ngkmax_ptr_save)
      nullify(ngk_ptr_save)
      nullify(vkl_ptr_save)
      nullify(vgkl_ptr_save)

      ! Backup pointers to default locations
      nmatmax_ptr_save => nmatmax_ptr
      ngkmax_ptr_save => ngkmax_ptr
      ngk_ptr_save => ngk_ptr
      vkl_ptr_save => vkl_ptr
      vgkl_ptr_save => vgkl_ptr

      ! Set default pointers to ket state quantities
      nmatmax_ptr => nmatmax1_ptr
      ngkmax_ptr => ngkmax1_ptr
      ngk_ptr => ngk1_ptr
      vkl_ptr => vkl1_ptr
      vgkl_ptr => vgkl1_ptr

      ! Call to getevecfv with changed default (G+)k-set pointers / matrix size
      evecfvt(:,:,:)=zzero
      !write(*,*) 'filext in b_getevecsv:', filext
      call getevecfv(vpl, vgpl, evecfvt)

     ! Restore default pointers
      nmatmax_ptr => nmatmax_ptr_save
      ngkmax_ptr => ngkmax_ptr_save
      ngk_ptr => ngk_ptr_save
      vkl_ptr => vkl_ptr_save
      vgkl_ptr => vgkl_ptr_save

    end subroutine b_getevecfv1
    !EOC

    subroutine b_match1(ngp, gpc, tpgpc, sfacgp, apwalm)
      use mod_ematptr
      use mod_eigensystem, only: nmatmax_ptr 
      use mod_Gkvector, only: ngkmax_ptr, ngk_ptr, vgkl_ptr
      use mod_kpoint, only: vkl_ptr
      use mod_atoms, only: natmtot
      use mod_APW_LO, only: apwordmax
      use mod_muffin_tin, only: lmmaxapw
      implicit none
      ! arguments
      integer, intent (in) :: ngp
      real (8), intent (in) :: gpc (ngkmax1_ptr)
      real (8), intent (in) :: tpgpc (2, ngkmax1_ptr)
      complex (8), intent (in) :: sfacgp (ngkmax1_ptr, natmtot)
      complex (8), intent (out) :: apwalm (ngkmax1_ptr, apwordmax, lmmaxapw, &
     & natmtot)
      ! Local variables
      integer, pointer :: nmatmax_ptr_save, ngkmax_ptr_save
      integer, pointer :: ngk_ptr_save(:,:)
      real(8), pointer :: vkl_ptr_save(:,:), vgkl_ptr_save(:,:,:,:)
      nullify(nmatmax_ptr_save)
      nullify(ngkmax_ptr_save)
      nullify(ngk_ptr_save)
      nullify(vkl_ptr_save)
      nullify(vgkl_ptr_save)

      ! Backup pointers to default locations
      nmatmax_ptr_save => nmatmax_ptr
      ngkmax_ptr_save => ngkmax_ptr
      ngk_ptr_save => ngk_ptr
      vkl_ptr_save => vkl_ptr
      vgkl_ptr_save => vgkl_ptr

      ! Set default pointers to ket state quantities
      nmatmax_ptr => nmatmax1_ptr
      ngkmax_ptr => ngkmax1_ptr
      ngk_ptr => ngk1_ptr
      vkl_ptr => vkl1_ptr
      vgkl_ptr => vgkl1_ptr

      ! Call to match with changed default (G+)k-set pointers / matrix size
      call match(ngp, gpc, tpgpc, sfacgp, apwalm)

     ! Restore default pointers
      nmatmax_ptr => nmatmax_ptr_save
      ngkmax_ptr => ngkmax_ptr_save
      ngk_ptr => ngk_ptr_save
      vkl_ptr => vkl_ptr_save
      vgkl_ptr => vgkl_ptr_save
    
    end subroutine b_match1
    
    subroutine b_match0(ngp, gpc, tpgpc, sfacgp, apwalm)
      use mod_ematptr
      use mod_eigensystem, only: nmatmax_ptr 
      use mod_Gkvector, only: ngkmax_ptr, ngk_ptr, vgkl_ptr
      use mod_kpoint, only: vkl_ptr
      use mod_atoms, only: natmtot
      use mod_APW_LO, only: apwordmax
      use mod_muffin_tin, only: lmmaxapw
      use mod_misc, only: filext
      use modxs, only: filext0
      implicit none
      ! arguments
      integer, intent (in) :: ngp
      real (8), intent (in) :: gpc (ngkmax0_ptr)
      real (8), intent (in) :: tpgpc (2, ngkmax0_ptr)
      complex (8), intent (in) :: sfacgp (ngkmax0_ptr, natmtot)
      complex (8), intent (out) :: apwalm (ngkmax0_ptr, apwordmax, lmmaxapw, &
     & natmtot)
      ! Local variables
      character(256) :: filext_save
      integer, pointer :: nmatmax_ptr_save, ngkmax_ptr_save
      integer, pointer :: ngk_ptr_save(:,:)
      real(8), pointer :: vkl_ptr_save(:,:), vgkl_ptr_save(:,:,:,:)
      nullify(nmatmax_ptr_save)
      nullify(ngkmax_ptr_save)
      nullify(ngk_ptr_save)
      nullify(vkl_ptr_save)
      nullify(vgkl_ptr_save)

      ! Backup pointers to default locations
      nmatmax_ptr_save => nmatmax_ptr
      ngkmax_ptr_save => ngkmax_ptr
      ngk_ptr_save => ngk_ptr
      vkl_ptr_save => vkl_ptr
      vgkl_ptr_save => vgkl_ptr
      filext_save=filext
      ! Set default pointers to ket state quantities
      nmatmax_ptr => nmatmax0_ptr
      ngkmax_ptr => ngkmax0_ptr
      ngk_ptr => ngk0_ptr
      vkl_ptr => vkl0_ptr
      vgkl_ptr => vgkl0_ptr

      filext = filext0
      ! Call to match with changed default (G+)k-set pointers / matrix size
      call match(ngp, gpc, tpgpc, sfacgp, apwalm)
      
      filext = filext_save
     ! Restore default pointers
      nmatmax_ptr => nmatmax_ptr_save
      ngkmax_ptr => ngkmax_ptr_save
      ngk_ptr => ngk_ptr_save
      vkl_ptr => vkl_ptr_save
      vgkl_ptr => vgkl_ptr_save
      filext=filext_save
    end subroutine b_match0
    
    subroutine b_wavefmt1(lrstp, lmax, is, ia, ngp, apwalm, evecfv, ld, wfmt)
      use mod_ematptr, only: nmatmax1_ptr, ngkmax1_ptr, ngk1_ptr, vkl1_ptr, vgkl1_ptr
      use mod_eigensystem, only: nmatmax_ptr 
      use mod_Gkvector, only: ngkmax_ptr, ngk_ptr, vgkl_ptr
      use mod_kpoint, only: vkl_ptr
      use mod_atoms, only: natmtot
      use mod_APW_LO, only: apwordmax
      use mod_muffin_tin, only: lmmaxapw
      
      implicit none
      ! arguments
      integer, intent (in) :: lrstp
      integer, intent (in) :: lmax
      integer, intent (in) :: is
      integer, intent (in) :: ia
      integer, intent (in) :: ngp
      complex (8), intent (in) :: apwalm (ngkmax1_ptr, apwordmax, lmmaxapw, &
        & natmtot)
      complex (8), intent (in) :: evecfv (nmatmax1_ptr)
      integer, intent (in) :: ld
      complex (8), Intent (out) :: wfmt (ld,*)
      ! Local variables
      integer, pointer :: nmatmax_ptr_save, ngkmax_ptr_save
      integer, pointer :: ngk_ptr_save(:,:)
      real(8), pointer :: vkl_ptr_save(:,:), vgkl_ptr_save(:,:,:,:)
      nullify(nmatmax_ptr_save)
      nullify(ngkmax_ptr_save)
      nullify(ngk_ptr_save)
      nullify(vkl_ptr_save)
      nullify(vgkl_ptr_save)

      ! Backup pointers to default locations
      nmatmax_ptr_save => nmatmax_ptr
      ngkmax_ptr_save => ngkmax_ptr
      ngk_ptr_save => ngk_ptr
      vkl_ptr_save => vkl_ptr
      vgkl_ptr_save => vgkl_ptr

      ! Set default pointers to ket state quantities
      nmatmax_ptr => nmatmax1_ptr
      ngkmax_ptr => ngkmax1_ptr
      ngk_ptr => ngk1_ptr
      vkl_ptr => vkl1_ptr
      vgkl_ptr => vgkl1_ptr

      ! Call to match with changed default (G+)k-set pointers / matrix size
      call wavefmt(lrstp, lmax, is, ia, ngp, apwalm, evecfv, ld, wfmt)


     ! Restore default pointers
      nmatmax_ptr => nmatmax_ptr_save
      ngkmax_ptr => ngkmax_ptr_save
      ngk_ptr => ngk_ptr_save
      vkl_ptr => vkl_ptr_save
      vgkl_ptr => vgkl_ptr_save
    
    end subroutine b_wavefmt1

    subroutine b_wavefmt0(lrstp, lmax, is, ia, ngp, apwalm, evecfv, ld, wfmt)
      use mod_ematptr, only: nmatmax0_ptr, ngkmax0_ptr, ngk0_ptr, vkl0_ptr, vgkl0_ptr
      use mod_eigensystem, only: nmatmax_ptr 
      use mod_Gkvector, only: ngkmax_ptr, ngk_ptr, vgkl_ptr
      use mod_kpoint, only: vkl_ptr
      use mod_atoms, only: natmtot
      use mod_APW_LO, only: apwordmax
      use mod_muffin_tin, only: lmmaxapw
      
      implicit none
      ! arguments
      integer, intent (in) :: lrstp
      integer, intent (in) :: lmax
      integer, intent (in) :: is
      integer, intent (in) :: ia
      integer, intent (in) :: ngp
      complex (8), intent (in) :: apwalm (ngkmax0_ptr, apwordmax, lmmaxapw, &
        & natmtot)
      complex (8), intent (in) :: evecfv (nmatmax0_ptr)
      integer, intent (in) :: ld
      complex (8), Intent (out) :: wfmt (ld,*)
      ! Local variables
      integer, pointer :: nmatmax_ptr_save, ngkmax_ptr_save
      integer, pointer :: ngk_ptr_save(:,:)
      real(8), pointer :: vkl_ptr_save(:,:), vgkl_ptr_save(:,:,:,:)
      nullify(nmatmax_ptr_save)
      nullify(ngkmax_ptr_save)
      nullify(ngk_ptr_save)
      nullify(vkl_ptr_save)
      nullify(vgkl_ptr_save)

      ! Backup pointers to default locations
      nmatmax_ptr_save => nmatmax_ptr
      ngkmax_ptr_save => ngkmax_ptr
      ngk_ptr_save => ngk_ptr
      vkl_ptr_save => vkl_ptr
      vgkl_ptr_save => vgkl_ptr

      ! Set default pointers to ket state quantities
      nmatmax_ptr => nmatmax0_ptr
      ngkmax_ptr => ngkmax0_ptr
      ngk_ptr => ngk0_ptr
      vkl_ptr => vkl0_ptr
      vgkl_ptr => vgkl0_ptr

      ! Call to match with changed default (G+)k-set pointers / matrix size
      call wavefmt(lrstp, lmax, is, ia, ngp, apwalm, evecfv, ld, wfmt)


     ! Restore default pointers
      nmatmax_ptr => nmatmax_ptr_save
      ngkmax_ptr => ngkmax_ptr_save
      ngk_ptr => ngk_ptr_save
      vkl_ptr => vkl_ptr_save
      vgkl_ptr => vgkl_ptr_save
    
    end subroutine b_wavefmt0
    
    !BOP
    ! !ROUTINE: b_getevecfv0
    ! !INTERFACE:
    subroutine b_getevecsv0(ik, evecsvt)
    ! !USES:
      use mod_kpoint, only: vkl_ptr
      use mod_Gkvector, only: ngkmax_ptr, vgkl_ptr, ngk_ptr
      use mod_eigensystem, only: nmatmax_ptr
      use mod_misc, only: filext
      use modxs, only: filext0
      use mod_eigenvalue_occupancy, only: nstsv
      use mod_ematptr
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    !   integer :: ik        ! k-point index
    ! OUT:
    !   complex(8) :: evecsvt(nstsv, nstsv) ! Eigenvectors at that k-point
    !
    ! !DESCRIPTION:
    !   This routine is a wrapper for {\tt getevecfv} that changes the $k$ and $G+k$ 
    !   quantities in {\tt mod\_kpoint} and {\tt mod\_Gkvector} to the corresponding
    !   quantities saved in {\tt modxs} (nmatmax0, vkl0,  ngk0, etc.), changes
    !   the file extension in {\tt mod\_misc} accordingly, reads in the Eigenvector
    !   and finally restores the original state.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Extended to 2nd variation. (Vorwerk)
    !EOP
    !BOC
      implicit none

      ! Arguments
      integer, intent(in)     :: ik
      complex(8), intent(out) :: evecsvt(nstsv, nstsv)

      ! Local variables
      integer, pointer :: nmatmax_ptr_save, ngkmax_ptr_save
      integer, pointer :: ngk_ptr_save(:,:)
      real(8), pointer :: vkl_ptr_save(:,:), vgkl_ptr_save(:,:,:,:)
      character(256) :: filext_save

      nullify(nmatmax_ptr_save)
      nullify(ngkmax_ptr_save)
      nullify(ngk_ptr_save)
      nullify(vkl_ptr_save)
      nullify(vgkl_ptr_save)

      ! Backup pointers to default locations
      nmatmax_ptr_save => nmatmax_ptr
      ngkmax_ptr_save => ngkmax_ptr
      ngk_ptr_save => ngk_ptr
      vkl_ptr_save => vkl_ptr
      vgkl_ptr_save => vgkl_ptr

      filext_save = filext

      ! Set default pointers to bra state quantities
      nmatmax_ptr => nmatmax0_ptr
      ngkmax_ptr => ngkmax0_ptr
      ngk_ptr => ngk0_ptr
      vkl_ptr => vkl0_ptr
      vgkl_ptr => vgkl0_ptr

      filext = filext0
      ! Call to getevecfv with changed default (G+)k-set pointers / matrix size
      call getevecsv(vkl_ptr(1:3,ik), evecsvt)
      ! Restore default pointers
      nmatmax_ptr => nmatmax_ptr_save
      ngkmax_ptr => ngkmax_ptr_save
      ngk_ptr => ngk_ptr_save
      vkl_ptr => vkl_ptr_save
      vgkl_ptr => vgkl_ptr_save

      filext = filext_save

    end subroutine b_getevecsv0

    !BOP
    ! !ROUTINE: b_getevecsv1
    ! !INTERFACE:
    subroutine b_getevecsv1(ik,evecsvt)
    ! !USES:
      use mod_kpoint, only: vkl_ptr
      use mod_Gkvector, only: ngkmax_ptr, vgkl_ptr, ngk_ptr
      use mod_eigensystem, only: nmatmax_ptr
      use mod_eigenvalue_occupancy, only: nstsv
      use mod_ematptr
      ! !INPUT/OUTPUT PARAMETERS:
      ! IN:
      !   integer :: ik          ! k-point index
      ! OUT:
      !   complex(8) :: evecsvt(nstsv, nstsv) ! Eigenvectors at that k-point
      !
      ! !DESCRIPTION:
      !   This routine is a wrapper for {\tt getevecfv} that changes the $k$ and $G+k$ 
      !   quantities in {\tt mod\_kpoint} and {\tt mod\_Gkvector} to the corresponding
      !   quantities saved in {\tt modxs} (nmatmax0, vkl0,  ngk0, etc.), changes
      !   the file extension in {\tt mod\_misc} accordingly, reads in the Eigenvector
      !   and finally restores the original state.
      !
      ! !REVISION HISTORY:
      !   Added to documentation scheme. (Aurich)
      !   Extended to 2nd variation (Vorwerk)
      !
      !EOP
      !BOC
      implicit none

      ! Arguments
      integer, intent(in)     :: ik
      complex(8), intent(out) :: evecsvt(nstsv, nstsv)

      ! Local variables
      integer, pointer :: nmatmax_ptr_save, ngkmax_ptr_save
      integer, pointer :: ngk_ptr_save(:,:)
      real(8), pointer :: vkl_ptr_save(:,:), vgkl_ptr_save(:,:,:,:)
      nullify(nmatmax_ptr_save)
      nullify(ngkmax_ptr_save)
      nullify(ngk_ptr_save)
      nullify(vkl_ptr_save)
      nullify(vgkl_ptr_save)

      ! Backup pointers to default locations
      nmatmax_ptr_save => nmatmax_ptr
      ngkmax_ptr_save => ngkmax_ptr
      ngk_ptr_save => ngk_ptr
      vkl_ptr_save => vkl_ptr
      vgkl_ptr_save => vgkl_ptr

      ! Set default pointers to ket state quantities
      nmatmax_ptr => nmatmax1_ptr
      ngkmax_ptr => ngkmax1_ptr
      ngk_ptr => ngk1_ptr
      vkl_ptr => vkl1_ptr
      vgkl_ptr => vgkl1_ptr

      ! Call to getevecfv with changed default (G+)k-set pointers / matrix size
      call getevecsv(vkl_ptr(:,ik), evecsvt)

      ! Resorte default pointers
      nmatmax_ptr => nmatmax_ptr_save
      ngkmax_ptr => ngkmax_ptr_save
      ngk_ptr => ngk_ptr_save
      vkl_ptr => vkl_ptr_save
      vgkl_ptr => vgkl_ptr_save

    end subroutine b_getevecsv1
!--------------------------------------------------------------------------------
  subroutine b_wavefmtsv1(lrstp ,lmax ,is ,ia , ngp, ist1, apwalm, evecfv, evecsv, wfmtsv)

    use modinput
    use mod_muffin_tin, only: nrmt, lmmaxapw, nrmtmax
    use mod_atoms, only: idxas, natmtot
    use mod_eigenvalue_occupancy, only: nstfv, nstsv
    use mod_ematptr, only: ngkmax1_ptr
    use mod_APW_LO, only: apwordmax
    use mod_eigensystem, only: nmatmax
    use mod_spin, only: nspinor
    !use modmain
    implicit none

    ! arguments
    integer, intent (in) :: lrstp
    integer, intent (in) :: lmax
    integer, intent (in) :: is
    integer, intent (in) :: ia
    integer, intent (in) :: ngp
    integer,    intent(in)  :: ist1
    complex(8), intent(in)  :: apwalm(ngkmax1_ptr,apwordmax,lmmaxapw,natmtot)
    complex(8), intent(in)  :: evecfv(nmatmax,nstfv)
    complex(8), intent(in)  :: evecsv(nstsv,nstsv)
    complex(8), intent(out) :: wfmtsv(lmmaxapw,nrmtmax,nspinor)
    ! local variables
    integer    :: ispn, ias
    integer    :: i, j, n, ist, igp
    complex(8) :: zt1
    ! allocatable arrays
    logical,    allocatable :: done (:)
    complex(8), allocatable :: wfmt1(:,:,:)

    wfmtsv(:,:,:) = 0.d0

    allocate(done(nstfv))
    allocate (wfmt1(lmmaxapw,nrmtmax,nstsv))

    wfmt1(:,:,:)=0.0d0
    n = lmmaxapw*nrmt(is)
    ias = idxas(ia,is)
    done(:) = .false.
    ! generate spinor wavefunction from second-variational eigenvectors
    i = 0
    do ispn = 1, nspinor
      do ist = 1, nstfv
        i = i + 1
        zt1 = evecsv(i,ist1)
        if (abs(zt1)>1.d-8) then
          if (.not.done(ist)) then
            call b_wavefmt1(lrstp, lmax, is, ia, &
              & ngp, apwalm, evecfv(:,ist), lmmaxapw, wfmt1(:,:,ist))
            done(ist) = .true.
          end if
          ! add to spinor wavefunction
          call zaxpy(n, zt1, wfmt1(:,:,ist), 1, &
          &          wfmtsv(:,:,ispn), 1)
        end if
      end do ! ist
    end do ! ispn
    deallocate(done, wfmt1)
    return
 end subroutine b_wavefmtsv1
!--------------------------------------------------------------------------------
  subroutine b_wavefmtsv0(lrstp ,lmax ,is ,ia , ngp, ist1, apwalm, evecfv, evecsv, wfmtsv)

    use modinput
    use mod_muffin_tin, only: nrmt, lmmaxapw, nrmtmax
    use mod_atoms, only: idxas, natmtot
    use mod_eigenvalue_occupancy, only: nstfv, nstsv
    use mod_ematptr, only: ngkmax0_ptr, nmatmax0_ptr
    use mod_APW_LO, only: apwordmax
    use mod_spin, only: nspinor
    !use modmain
    implicit none

    ! arguments
    integer, intent (in) :: lrstp
    integer, intent (in) :: lmax
    integer, intent (in) :: is
    integer, intent (in) :: ia
    integer, intent (in) :: ngp
    integer,    intent(in)  :: ist1
    complex(8), intent(in)  :: apwalm(ngkmax0_ptr,apwordmax,lmmaxapw,natmtot)
    complex(8), intent(in)  :: evecfv(nmatmax0_ptr,nstfv)
    complex(8), intent(in)  :: evecsv(nstsv,nstsv)
    complex(8), intent(out) :: wfmtsv(lmmaxapw,nrmtmax,nspinor)
    ! local variables
    integer    :: ispn, ias
    integer    :: i, j, n, ist, igp
    complex(8) :: zt1
    ! allocatable arrays
    logical,    allocatable :: done (:)
    complex(8), allocatable :: wfmt1(:,:,:)

    wfmtsv(:,:,:) = 0.d0

    allocate(done(nstfv))
    allocate (wfmt1(lmmaxapw,nrmtmax,nstsv))

    wfmt1(:,:,:)=0.0d0
    n = lmmaxapw*nrmt(is)
    ias = idxas(ia,is)
    done(:) = .false.
    ! generate spinor wavefunction from second-variational eigenvectors
    i = 0
    do ispn = 1, nspinor
      do ist = 1, nstfv
        i = i + 1
        zt1 = evecsv(i,ist1)
        if (abs(zt1)>1.d-8) then
          if (.not.done(ist)) then
            call b_wavefmt0(lrstp, lmax, is, ia, &
              & ngp, apwalm, evecfv(:,ist), lmmaxapw, wfmt1(:,:,ist))
            done(ist) = .true.
          end if
          ! add to spinor wavefunction
          call zaxpy(n, zt1, wfmt1(:,:,ist), 1, &
          &          wfmtsv(:,:,ispn), 1)
        end if
      end do ! ist
    end do ! ispn
    deallocate(done, wfmt1)
    return
 end subroutine b_wavefmtsv0
!--------------------------------------------------------------------------------
  subroutine b_wavefmtsv(lrstp ,lmax ,is ,ia , ngp, ist1, apwalm, evecfv, evecsv, wfmtsv)

    use modinput
    use mod_muffin_tin, only: nrmt, lmmaxapw, nrmtmax
    use mod_atoms, only: idxas, natmtot
    use mod_eigenvalue_occupancy, only: nstfv, nstsv
    use mod_ematptr, only: ngkmax0_ptr
    use mod_APW_LO, only: apwordmax
    use mod_eigensystem, only: nmatmax
    use mod_spin, only: nspinor
    use mod_Gkvector, only: ngkmax
    !use modmain
    implicit none

    ! arguments
    integer, intent (in) :: lrstp
    integer, intent (in) :: lmax
    integer, intent (in) :: is
    integer, intent (in) :: ia
    integer, intent (in) :: ngp
    integer,    intent(in)  :: ist1
    complex(8), intent(in)  :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
    complex(8), intent(in)  :: evecfv(nmatmax,nstfv)
    complex(8), intent(in)  :: evecsv(nstsv,nstsv)
    complex(8), intent(out) :: wfmtsv(lmmaxapw,nrmtmax,nspinor)
    ! local variables
    integer    :: ispn, ias
    integer    :: i, j, n, ist, igp
    complex(8) :: zt1
    ! allocatable arrays
    logical,    allocatable :: done (:)
    complex(8), allocatable :: wfmt1(:,:,:)

    wfmtsv(:,:,:) = 0.d0

    allocate(done(nstfv))
    allocate (wfmt1(lmmaxapw,nrmtmax,nstsv))

    wfmt1(:,:,:)=0.0d0
    n = lmmaxapw*nrmt(is)
    ias = idxas(ia,is)
    done(:) = .false.
    ! generate spinor wavefunction from second-variational eigenvectors
    i = 0
    do ispn = 1, nspinor
      do ist = 1, nstfv
        i = i + 1
        zt1 = evecsv(i,ist1)
        if (abs(zt1)>1.d-8) then
          if (.not.done(ist)) then
            call wavefmt(lrstp, lmax, is, ia, &
              & ngp, apwalm, evecfv(:,ist), lmmaxapw, wfmt1(:,:,ist))
            done(ist) = .true.
          end if
          ! add to spinor wavefunction
          call zaxpy(n, zt1, wfmt1(:,:,ist), 1, &
          &          wfmtsv(:,:,ispn), 1)
        end if
      end do ! ist
    end do ! ispn
    deallocate(done, wfmt1)
    return
 end subroutine b_wavefmtsv
end module
