! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
module m_getgrst
  contains
    
    !BOP
    ! !ROUTINE: getevecfv0
    ! !INTERFACE:
    subroutine getevecfv0(vpl, vgpl, evecfvt)
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

    end subroutine getevecfv0
    !EOC
    
    !BOP
    ! !ROUTINE: getevecfv1
    ! !INTERFACE:
    subroutine getevecfv1(vpl, vgpl, evecfvt)
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
      !write(*,*) 'filext in getevecsv:', filext
      call getevecfv(vpl, vgpl, evecfvt)

     ! Restore default pointers
      nmatmax_ptr => nmatmax_ptr_save
      ngkmax_ptr => ngkmax_ptr_save
      ngk_ptr => ngk_ptr_save
      vkl_ptr => vkl_ptr_save
      vgkl_ptr => vgkl_ptr_save

    end subroutine getevecfv1
    !EOC

    subroutine match1(ngp, gpc, tpgpc, sfacgp, apwalm)
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
    
    end subroutine match1
    
    subroutine match0(ngp, gpc, tpgpc, sfacgp, apwalm)
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
    end subroutine match0
    
    subroutine wavefmt1(lrstp, lmax, is, ia, ngp, apwalm, evecfv, ld, wfmt)
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
    
    end subroutine wavefmt1

end module
