!> This module contains procedures for the calculation of the response
!> of the macroscopic polarzation upon a generic perturbation.
module dfpt_polarization
  use dfpt_variables

  use precision, only: dp

  implicit none
  private

  !> radial integrals times Gaunt coeffiencients of plane wave operator
  complex(dp), allocatable, public :: PWmat_mt_basis(:,:,:)

  public :: dfpt_pol_free, dfpt_pol_init_p

contains
    
  !> Free memory from unneeded variables.
  subroutine dfpt_pol_free
    if( allocated( PWmat_mt_basis ) ) deallocate( PWmat_mt_basis )
  end subroutine dfpt_pol_free

  !> Initialize plane wave matrix elements for a given polarization direction / 
  !> displacement vector \({\bf b}\).
  !>
  !> Calculate corresponding radial integrals times Gaunt coefficients.
  subroutine dfpt_pol_init_p( vbc )
    use matrix_elements
    use math_utils, only: plane_wave_in_spherical_harmonics
    use constants, only: zzero, zone
    use mod_muffin_tin, only: nrmtmax, nrmt
    use mod_atoms, only: nspecies, natoms, idxas, spr, atposc
    !> displacement vector \({\bf b}\) along polarization direction in Cartesian coordinates
    real(dp), intent(in) :: vbc(3)

    integer :: is, ia, ias
    real(dp) :: dotp
    complex(dp) :: z1

    complex(dp), allocatable :: pwrfun(:,:)

    call me_mt_alloc( PWmat_mt_basis )

    do is = 1, nspecies
      call plane_wave_in_spherical_harmonics( vbc, spr(1:nrmt(is), is), dfpt_lmaxvr, pwrfun )
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        dotp = dot_product( vbc, atposc(:, ia, is) )
        z1 = cmplx( cos( dotp ), sin( dotp ), dp )
        pwrfun = pwrfun * z1
        call me_mt_prepare( is, ias, dfpt_lmaxvr, zone, pwrfun, zzero, PWmat_mt_basis(:, :, ias) )
        pwrfun = pwrfun / z1
      end do
    end do

    if( allocated( pwrfun ) ) deallocate( pwrfun )
  end subroutine dfpt_pol_init_p

end module dfpt_polarization
