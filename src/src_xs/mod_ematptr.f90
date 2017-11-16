module mod_ematptr
  ! Default arrays k,G arrays
  use mod_kpoint, only: nkpt, vkl
  use mod_Gkvector, only: ngkmax, ngk, igkig, vgkl,&
                        & vgkc, gkc, tpgkc, sfacgk
  use mod_eigensystem, only: nmatmax, nmat
  ! modxs mirror images of k,G arrays
  use modxs, only: vkl0, ngk0, igkig0, vgkl0,&
                 & vgkc0, gkc0, tpgkc0, sfacgk0,&
                 & nmat0, nkpt0, ngkmax0, nmatmax0
  ! bra and ket eigenvectors
  use modxs, only: evecfv, evecfvb, evecfv0, evecfv0b,&
                 & apwcmt, apwcmt0

  implicit none

  ! Pointer for k and k' quantities
  integer(4), pointer :: nkpt0_ptr, nkpt1_ptr
  integer(4), pointer :: nmatmax0_ptr, nmatmax1_ptr
  integer(4), pointer :: ngkmax0_ptr, ngkmax1_ptr
  integer(4), pointer :: nmat0_ptr(:,:), nmat1_ptr(:,:)
  integer(4), pointer :: ngk0_ptr(:,:), ngk1_ptr(:,:)
  integer(4), pointer :: igkig0_ptr(:,:,:), igkig1_ptr(:,:,:)
  real(8), pointer :: vkl0_ptr(:,:), vkl1_ptr(:,:)
  real(8), pointer :: vgkl0_ptr(:,:,:,:), vgkl1_ptr(:,:,:,:)
  real(8), pointer :: vgkc0_ptr(:,:,:,:), vgkc1_ptr(:,:,:,:)
  real(8), pointer :: gkc0_ptr(:,:,:), gkc1_ptr(:,:,:)
  real(8), pointer :: tpgkc0_ptr(:,:,:,:), tpgkc1_ptr(:,:,:,:)
  complex(8), pointer :: sfacgk0_ptr(:,:,:,:), sfacgk1_ptr(:,:,:,:)

  ! Bra and ket eigenvectors
  complex(8), pointer :: evecfv0_ptr(:,:,:), evecfv1_ptr(:,:,:)
  complex(8), pointer :: apwcmt0_ptr(:,:,:,:), apwcmt1_ptr(:,:,:,:)

  ! Link map between k and k'=k+q
  integer(4), pointer :: ikmapikq_ptr(:, :)

  contains

    subroutine setptr01()

      ! Set 0 pointers (bra) to the modxs::*0 arrays
      nkpt0_ptr => nkpt0
      nmatmax0_ptr => nmatmax0
      ngkmax0_ptr => ngkmax0
      nmat0_ptr => nmat0
      ngk0_ptr => ngk0
      igkig0_ptr => igkig0
      vkl0_ptr => vkl0
      vgkl0_ptr => vgkl0
      vgkc0_ptr => vgkc0
      gkc0_ptr => gkc0
      tpgkc0_ptr => tpgkc0
      sfacgk0_ptr => sfacgk0

      evecfv0_ptr => evecfv0
      apwcmt0_ptr=> apwcmt0

      ! Set 1 pointers (ket) to the default arrays
      nkpt1_ptr => nkpt
      nmatmax1_ptr => nmatmax
      ngkmax1_ptr => ngkmax
      nmat1_ptr => nmat
      ngk1_ptr => ngk
      igkig1_ptr => igkig
      vkl1_ptr => vkl
      vgkl1_ptr => vgkl
      vgkc1_ptr => vgkc
      gkc1_ptr => gkc
      tpgkc1_ptr => tpgkc
      sfacgk1_ptr => sfacgk

      evecfv1_ptr => evecfv
      apwcmt1_ptr => apwcmt

    end subroutine setptr01

    subroutine setptr10()

      ! Set 0 pointers (bra) to the default arrays
      nkpt0_ptr => nkpt
      nmatmax0_ptr => nmatmax
      ngkmax0_ptr => ngkmax
      nmat0_ptr => nmat
      ngk0_ptr => ngk
      igkig0_ptr => igkig
      vkl0_ptr => vkl
      vgkl0_ptr => vgkl
      vgkc0_ptr => vgkc
      gkc0_ptr => gkc
      tpgkc0_ptr => tpgkc
      sfacgk0_ptr => sfacgk

      evecfv0_ptr => evecfv
      apwcmt0_ptr => apwcmt

      ! Set 1 pointers (ket) to the modxs::*0 arrays
      nkpt1_ptr => nkpt0
      nmatmax1_ptr => nmatmax0
      ngkmax1_ptr => ngkmax0
      nmat1_ptr => nmat0
      ngk1_ptr => ngk0
      igkig1_ptr => igkig0
      vkl1_ptr => vkl0
      vgkl1_ptr => vgkl0
      vgkc1_ptr => vgkc0
      gkc1_ptr => gkc0
      tpgkc1_ptr => tpgkc0
      sfacgk1_ptr => sfacgk0

      evecfv1_ptr => evecfv0
      apwcmt1_ptr => apwcmt0

    end subroutine setptr10

    subroutine setptr00()

      ! Set 0 pointers (bra) to the modxs::*0 arrays
      nkpt0_ptr => nkpt0
      nmatmax0_ptr => nmatmax0
      ngkmax0_ptr => ngkmax0
      nmat0_ptr => nmat0
      ngk0_ptr => ngk0
      igkig0_ptr => igkig0
      vkl0_ptr => vkl0
      vgkl0_ptr => vgkl0
      vgkc0_ptr => vgkc0
      gkc0_ptr => gkc0
      tpgkc0_ptr => tpgkc0
      sfacgk0_ptr => sfacgk0

      evecfv0_ptr => evecfv0
      apwcmt0_ptr => apwcmt0

      ! Set 1 pointers (ket) to the modxs::*0 arrays
      nkpt1_ptr => nkpt0
      nmatmax1_ptr => nmatmax0
      ngkmax1_ptr => ngkmax0
      nmat1_ptr => nmat0
      ngk1_ptr => ngk0
      igkig1_ptr => igkig0
      vkl1_ptr => vkl0
      vgkl1_ptr => vgkl0
      vgkc1_ptr => vgkc0
      gkc1_ptr => gkc0
      tpgkc1_ptr => tpgkc0
      sfacgk1_ptr => sfacgk0

      evecfv1_ptr => evecfv0b
      apwcmt1_ptr => apwcmt

    end subroutine setptr00

    subroutine setptr11()

      ! Set 0 pointers (bra) to the default arrays
      nkpt0_ptr => nkpt
      nmatmax0_ptr => nmatmax
      ngkmax0_ptr => ngkmax
      nmat0_ptr => nmat
      ngk0_ptr => ngk
      igkig0_ptr => igkig
      vkl0_ptr => vkl
      vgkl0_ptr => vgkl
      vgkc0_ptr => vgkc
      gkc0_ptr => gkc
      tpgkc0_ptr => tpgkc
      sfacgk0_ptr => sfacgk

      evecfv0_ptr => evecfv
      apwcmt0_ptr => apwcmt0

      ! Set 1 pointers (ket) to the default arrays
      nkpt1_ptr => nkpt
      nmatmax1_ptr => nmatmax
      ngkmax1_ptr => ngkmax
      nmat1_ptr => nmat
      ngk1_ptr => ngk
      igkig1_ptr => igkig
      vkl1_ptr => vkl
      vgkl1_ptr => vgkl
      vgkc1_ptr => vgkc
      gkc1_ptr => gkc
      tpgkc1_ptr => tpgkc
      sfacgk1_ptr => sfacgk

      evecfv1_ptr => evecfvb
      apwcmt1_ptr => apwcmt

    end subroutine setptr11

end module mod_ematptr
