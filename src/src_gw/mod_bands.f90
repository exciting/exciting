
!--------------------------------!
! DFT groundstate related stuff  !
!--------------------------------!

module mod_bands
    implicit none
    
! (spin-collinear case) Spin index of states
    integer, allocatable :: spindex(:,:)
    
! Eigenvectors at k
    complex(8), allocatable :: eveck(:,:)
      
! Eigenvectors at k'=k-q
    complex(8), allocatable :: eveckp(:,:)
      
! Spherical harmonic expansion coefficients at k
    complex(8), allocatable :: eveckalm(:,:,:,:)
      
! Spherical harmonic expansion coefficients at k'=k-q
    complex(8), allocatable :: eveckpalm(:,:,:,:)
    
! Position of Valence Band Maximum (VBM)      
    integer(4) :: nomax
    integer(4) :: ikvbm
    
! Position of Conduction Band Minimum (CBM)      
    integer(4) :: numin
    integer(4) :: ikcbm
    
! Position of the direct v->c (optical) gap
    integer(4) :: ikvcm

! Number of states used to calculate the dielectric function
    integer(4) :: nstdf

! Number of states used to calculate the self-energy
    integer(4) :: nstse
    
! Metallicity flag
    logical :: metallic

! Lower and upper indexes of the degenerated states
    integer, allocatable :: n12dgn(:,:,:)

!---------------------------------------------------------------
! To be used in the interpolation routine (band structure plot) 
!---------------------------------------------------------------

! Input 
    integer(4) :: nkp1
    real(8), allocatable :: kvecs1(:,:)
    real(8), allocatable :: eks1(:,:), eqp1(:,:)

! Output (interpolated)
    integer(4) :: nkp2
    real(8), allocatable :: kvecs2(:,:)
    real(8), allocatable :: eks2(:,:), eqp2(:,:)
    
contains

    subroutine delete_bands
        if (allocated(eveck)) deallocate(eveck)
        if (allocated(eveckp)) deallocate(eveckp)
        if (allocated(eveckalm)) deallocate(eveckalm)
        if (allocated(eveckpalm)) deallocate(eveckpalm)
        if (allocated(n12dgn)) deallocate(n12dgn)
    end subroutine
    
end module
