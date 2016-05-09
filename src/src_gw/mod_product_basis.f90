
module mod_product_basis

    implicit none
    
    !----------------------------!
    !     mixed basis (general)  !
    !----------------------------!
    ! Size of the mixed basis
    integer(4) :: matsiz, matsizmax, mbsiz
 
    ! Matrix elements M^i_nm and \tilde{M}^i_nm
    complex(8), allocatable :: minmmat(:,:,:)
      
    ! Matrix elements M^i_cm and \tilde{M}^i_cm
    complex(8), allocatable :: micmmat(:,:,:)
      
    ! Matrix elements M^i_nc and \tilde{M}^i_nc
    complex(8), allocatable :: mincmat(:,:,:)
      
    ! Matrix elements between mixed functions and planewaves
    complex(8), allocatable :: mpwmix(:,:)
    
    !--------------------------------------!
    !     mixed basis (Muffin-Tins)        !
    !--------------------------------------!
    
    ! Upper size limit estimate of the number of possible radial function products
    integer(4) :: maxnup
    
    ! Actual number of radial function products (per atom)
    integer(4) :: nup
        
    ! Radial product functions (per atom)
    real(8), allocatable :: uprod(:,:)
    
    ! l,l' pairs of the product functions (per atom)
    integer(4), allocatable :: eles(:,:)
    
    ! Overlap matrix of the product functions (per atom)
    real(8), allocatable :: umat(:,:)
    
    ! Size of the local part of the mixed basis including LM combinations
    integer(4) :: locmatsiz
    
    ! maximum number of mixed functions per atom
    integer(4) :: lmixmax
    
    ! indexes of the local mixed basis functions
    integer(4), allocatable :: locmixind(:,:)
    
    ! Combined aNLM index of Mixed Product Basis functions
    ! mbindex(I,1) = is
    ! mbindex(I,2) = ia
    ! mbindex(I,3) = N 
    ! mbindex(I,4) = L
    ! mbindex(I,5) = M
    integer(4), allocatable :: mbindex(:,:)    
 
    !--------------

    ! Number of mixed radial functions per atom
    integer(4), allocatable :: nmix(:)
    
    ! Maximum number of radial functions per atom
    integer(4) :: maxnmix
    
    ! Maximum L of the mixed functions
    integer(4) :: maxbigl
    
    ! Radial mixed functions
    real(8), pointer :: umix(:,:,:)
    
    ! L quantum number of the mixed functions
    integer(4), pointer :: bigl(:,:)
    
    ! Maximum L of the mixed functions per atom
    integer(4), allocatable :: mbl(:)

    !--------------
      
    ! <umix(L)|ucore(l1)u(l2)> integrals
    real(8), allocatable :: bradketc(:,:,:,:,:,:)
    
    ! <umix(L)|u(l1)u(l2)> integrals
    real(8), allocatable :: bradketa(:,:,:,:,:,:,:)
    
    ! <umix(L)|ulo(l1)u(l2)> integrals
    real(8), allocatable :: bradketlo(:,:,:,:,:,:)
    
    ! <umix(L)|u_core> integrals
    real(8), allocatable :: umbucor(:,:,:)
    ! <umix(L)|u_apw> integrals
    real(8), allocatable :: umbuapw(:,:,:)
    ! <umix(L)|u_lo> integrals
    real(8), allocatable :: umbulor(:,:,:)
    
    ! the gaunt coefficients
    real(8), allocatable :: cgcoef(:)

    ! <umix(l)|r^(l+2)> integrals     
    real(8), allocatable :: rtl(:,:)
    
    ! <umix(l1)|r^(l1)/r^(l2+1)|umix(l2)> integrals
    real(8), allocatable :: rrint(:,:)
    
    !---------------------------------!
    !     mixed basis (Interstitial)  !
    !---------------------------------!
    
    !  number of G+q-vectors for the mixed basis
    integer(4), allocatable :: ngq(:)
    
    ! maximum number of G+q-vectors over all q-points
    integer ngqmax
    
    ! index from G+q-vectors to G-vectors
    integer, allocatable :: igqig(:,:)
    
    ! index from G-vectors to G+q-vectors
    integer, allocatable :: igigq(:,:)
    
    ! G+q-vectors in lattice coordinates
    real(8), allocatable :: vgql(:,:,:)
    
    ! G+q-vectors in Cartesian coordinates
    real(8), allocatable :: vgqc(:,:,:)
    
    ! Transformation matrix between IPW's and OIPW's
    complex(8), allocatable :: sgi(:,:)
    complex(8), allocatable :: sgi_fft(:,:)
    
    ! Matrix element between IPW's and PW's       
    complex(8), allocatable :: mpwipw(:,:)

    !---------------------------------------------------------------!
    ! Matrix representation of the symmetry operations in MB basis  !
    !---------------------------------------------------------------!
    complex(8), allocatable :: rotmat(:,:)

contains

    subroutine delete_product_basis()
        implicit none
        !------------------------------------
        if (allocated(nmix)) deallocate(nmix)
        if (associated(umix)) deallocate(umix)
        if (associated(bigl)) deallocate(bigl)
        if (allocated(mbl)) deallocate(mbl)
        !------------------------------------
        if (allocated(bradketc)) deallocate(bradketc)
        if (allocated(bradketa)) deallocate(bradketa)
        if (allocated(bradketlo)) deallocate(bradketlo)
        if (allocated(cgcoef)) deallocate(cgcoef)
        if (allocated(rtl)) deallocate(rtl)
        if (allocated(rrint)) deallocate(rrint)
        return
    end subroutine
    
end module

