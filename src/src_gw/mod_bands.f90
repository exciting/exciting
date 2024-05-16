
!--------------------------------!
! DFT groundstate related stuff  !
!--------------------------------!

module mod_bands

    use precision, only: i32, dp

    implicit none

    private
    public :: evalfv, &
              occfv, &
              eveck, &
              eveckp, &
              eveckalm, &
              eveckpalm, &
              nomax, &
              ikvbm, &
              numin, &
              ikcbm, &
              ikvcm, &
              nstdf, &
              nstse, &
              metallic, &
              nkp1, &
              kvecs1, &
              eks1, &
              eqp1, &
              nkp2, &
              kvecs2, &
              eks2, &
              eqp2, &
              delete_bands

! First-variational eigenvalues
    real(dp), allocatable :: evalfv(:,:)

! First-variational occupations
    real(dp), allocatable :: occfv(:,:)

! Eigenvectors at k
    complex(dp), allocatable :: eveck(:,:)
      
! Eigenvectors at k'=k-q
    complex(dp), allocatable :: eveckp(:,:)
      
! Spherical harmonic expansion coefficients at k
    complex(dp), allocatable :: eveckalm(:,:,:,:)
      
! Spherical harmonic expansion coefficients at k'=k-q
    complex(dp), allocatable :: eveckpalm(:,:,:,:)
    
! Position of Valence Band Maximum (VBM)      
    integer(i32) :: nomax
    integer(i32) :: ikvbm
    
! Position of Conduction Band Minimum (CBM)      
    integer(i32) :: numin
    integer(i32) :: ikcbm
    
! Position of the direct v->c (optical) gap
    integer(i32) :: ikvcm

! Number of states used to calculate the dielectric function
    integer(i32) :: nstdf

! Number of states used to calculate the self-energy
    integer(i32) :: nstse
    
! Metallicity flag
    logical :: metallic

!---------------------------------------------------------------
! To be used in the interpolation routine (band structure plot) 
!---------------------------------------------------------------

! Input 
    integer(i32) :: nkp1
    real(dp), allocatable :: kvecs1(:,:)
    real(dp), allocatable :: eks1(:,:), eqp1(:,:)

! Output (interpolated)
    integer(i32) :: nkp2
    real(dp), allocatable :: kvecs2(:,:)
    real(dp), allocatable :: eks2(:,:), eqp2(:,:)
    
contains

    subroutine delete_bands
        if (allocated(eveck)) deallocate(eveck)
        if (allocated(eveckp)) deallocate(eveckp)
        if (allocated(eveckalm)) deallocate(eveckalm)
        if (allocated(eveckpalm)) deallocate(eveckpalm)
    end subroutine
    
end module
