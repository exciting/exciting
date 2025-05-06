!> This contains all elements related to the product basis
module mod_product_basis

    use gw_io, only: read_from_file, write_to_file, build_file_name
    use precision, only: dp, i32
#include "offload.fpp"

    implicit none
    
    !----------------------------!
    !     mixed basis (general)  !
    !----------------------------!
    ! Size of the mixed basis
    integer(i32) :: matsiz, matsizmax, mbsiz
 
    ! Matrix elements M^i_nm and \tilde{M}^i_nm
    complex(dp), allocatable :: minmmat(:,:,:)
      
    ! Matrix elements M^i_cm and \tilde{M}^i_cm
    complex(dp), allocatable :: micmmat(:,:,:)
      
    ! Matrix elements M^i_nc and \tilde{M}^i_nc
    complex(dp), allocatable :: mincmat(:,:,:)
      
    ! Matrix elements between mixed functions and planewaves
    complex(dp), allocatable :: mpwmix(:,:)
    
    !--------------------------------------!
    !     mixed basis (Muffin-Tins)        !
    !--------------------------------------!
    
    ! Upper size limit estimate of the number of possible radial function products
    integer(i32) :: maxnup
    
    ! Actual number of radial function products (per atom)
    integer(i32) :: nup
        
    ! Radial product functions (per atom)
    real(dp), allocatable :: uprod(:,:)
    
    ! l,l' pairs of the product functions (per atom)
    integer(i32), allocatable :: eles(:,:)
    
    ! Overlap matrix of the product functions (per atom)
    real(dp), allocatable :: umat(:,:)
    
    ! Size of the local part of the mixed basis including LM combinations
    integer(i32) :: locmatsiz
    
    ! maximum number of mixed functions per atom
    integer(i32) :: lmixmax
    
    ! indexes of the local mixed basis functions
    integer(i32), allocatable :: locmixind(:,:)
    
    ! Combined aNLM index of Mixed Product Basis functions
    ! mbindex(I,1) = is
    ! mbindex(I,2) = ia
    ! mbindex(I,3) = N 
    ! mbindex(I,4) = L
    ! mbindex(I,5) = M
    integer(i32), allocatable :: mbindex(:,:)
 
    !--------------

    ! Number of mixed radial functions per atom
    integer(i32), allocatable :: nmix(:)
    
    ! Maximum number of radial functions per atom
    integer(i32) :: maxnmix
    
    ! Maximum L of the mixed functions
    integer(i32) :: maxbigl
    
    ! Radial mixed functions
    real(dp), pointer :: umix(:,:,:)
    
    ! L quantum number of the mixed functions
    integer(i32), pointer :: bigl(:,:)
    
    ! Maximum L of the mixed functions per atom
    integer(i32), allocatable :: mbl(:)

    !--------------
      
    ! <umix(L)|ucore(l1)u(l2)> integrals
    real(dp), allocatable :: bradketc(:,:,:,:,:,:)
    
    ! <umix(L)|u(l1)u(l2)> integrals
    real(dp), allocatable :: bradketa(:,:,:,:,:,:,:)
    
    ! <umix(L)|ulo(l1)u(l2)> integrals
    real(dp), allocatable :: bradketlo(:,:,:,:,:,:)
    
    ! <umix(L)|u_core> integrals
    real(dp), allocatable :: umbucor(:,:,:)
    ! <umix(L)|u_apw> integrals
    real(dp), allocatable :: umbuapw(:,:,:)
    ! <umix(L)|u_lo> integrals
    real(dp), allocatable :: umbulor(:,:,:)
    
    ! the gaunt coefficients
    real(dp), allocatable :: cgcoef(:)

    ! <umix(l)|r^(l+2)> integrals     
    real(dp), allocatable :: rtl(:,:)
    
    ! <umix(l1)|r^(l1)/r^(l2+1)|umix(l2)> integrals
    real(dp), allocatable :: rrint(:,:)
    
    !---------------------------------!
    !     mixed basis (Interstitial)  !
    !---------------------------------!
    
    !  number of G+q-vectors for the mixed basis
    integer(i32), allocatable :: ngq(:)
    
    ! maximum number of G+q-vectors over all q-points
    integer ngqmax
    
    ! index from G+q-vectors to G-vectors
    integer(i32), allocatable :: igqig(:,:)
    
    ! index from G-vectors to G+q-vectors
    integer(i32), allocatable :: igigq(:,:)
    
    ! G+q-vectors in lattice coordinates
    real(dp), allocatable :: vgql(:,:,:)
    
    ! G+q-vectors in Cartesian coordinates
    real(dp), allocatable :: vgqc(:,:,:)
    
    ! Transformation matrix between IPW's and OIPW's
    complex(dp), allocatable :: sgi(:,:)
    complex(dp), allocatable :: sgi_fft(:,:)
    
    ! Matrix element between IPW's and PW's       
    complex(dp), allocatable :: mpwipw(:,:)

    !---------------------------------------------------------------!
    ! Matrix representation of the symmetry operations in MB basis  !
    !---------------------------------------------------------------!
    complex(dp), allocatable :: rotmat(:,:)

    !> Name of output file where sgi is written to
    character(len=*), parameter :: basename_sgi = 'SGI_Q'

contains

    subroutine delete_product_basis()
        implicit none
        !------------------------------------
        if (allocated(nmix)) deallocate(nmix)
        if (associated(umix)) deallocate(umix)
        if (associated(bigl)) deallocate(bigl)
        if (allocated(mbl)) deallocate(mbl)
        !------------------------------------
        if (allocated(cgcoef)) deallocate(cgcoef)
        if (allocated(rtl)) deallocate(rtl)
        if (allocated(rrint)) deallocate(rrint)
        !------------------------------------
        if (allocated(bradketc)) then
          OMP_OFFLOAD target exit data map(delete: bradketc)
          deallocate(bradketc)
        end if

        if (allocated(bradketa)) then
          OMP_OFFLOAD target exit data map(delete: bradketa)
          deallocate(bradketa)
        end if

        if (allocated(bradketlo)) then
          OMP_OFFLOAD target exit data map(delete: bradketlo)
          deallocate(bradketlo)
        end if

        if (allocated(mbindex)) then
          OMP_OFFLOAD target exit data map(delete: mbindex)
          deallocate(mbindex)
        end if

        return
    end subroutine

    !> Write `sgi` to an output file
    subroutine write_sgi_to_file( int, file_format )
      !> integer that will be added to `basename_sgi` to form the output name
      integer, intent(in) :: int
      !> Format of the output file
      character(len=*), intent(in) :: file_format
    
      character(len=30)   :: file_name
      
      call build_file_name( basename_sgi, int, file_name )
      call write_to_file( file_name, sgi, file_format )
    
    end subroutine


    !> Read `sgi` stored in a file
    subroutine read_sgi_from_file( int, file_format )
      !> Integer that will be added to `basename_sgi` to form the output name
      integer, intent(in) :: int
      !> Format of the file to be read from
      character(len=*), intent(in) :: file_format
    
      character(len=30)   :: file_name
      
      call build_file_name( basename_sgi, int, file_name )
      call read_from_file( file_name, sgi, file_format )
    
    end subroutine
    
end module

