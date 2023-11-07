!> This module carries global DFPT electric field variables that are used across multiple modules.
module efield_variables
  use dfpt_variables

  use precision, only: dp
  use modmpi
  use asserts, only: assert

  implicit none
  private

  ! ** auxiliary quantities
  !> spherical Bessel functions \(j_l(|{\bf G}|R_\alpha)\)
  real(dp), allocatable, public :: ef_jlgr(:,:,:)
  !> spherical harmonics \(Y_{lm}(\widehat{{\bf G}})\)
  complex(dp), allocatable, public :: ef_ylmg(:,:)
  !> structure factors \({\rm e}^{{\rm i} ({\bf G})\cdot{\bf \tau}_\alpha}\)
  complex(dp), allocatable, public :: ef_sfacg(:,:)

  public :: ef_var_init, ef_var_free

contains

  subroutine ef_var_init
    use mod_atoms, only: natmtot, nspecies, natoms, idxas, atposc
    use mod_muffin_tin, only: rmt
    use modinput

    integer :: is, ia, ias, ig, lmmax
    real(dp) :: v(3), gc, t1

    lmmax = (dfpt_lmaxmp + 1)**2
    allocate( ef_jlgr(0:(dfpt_lmaxmp+input%groundstate%npsden+1), dfpt_Gset%ngvec, nspecies) )
    allocate( ef_ylmg(lmmax, dfpt_Gset%ngvec) )
    allocate( ef_sfacg(dfpt_Gset%ngvec, natmtot) )
    
!$omp parallel default( shared ) private( ig, v, gc, is, ia, ias, t1 )
!$omp do
    do ig = 1, dfpt_Gset%ngvec
      v = dfpt_Gset%vgc(:, ig)
      gc = dfpt_Gset%gc(ig)
      call ylm( v, dfpt_lmaxmp, ef_ylmg(:, ig) )
      do is = 1, nspecies
        call sbessel( dfpt_lmaxmp+input%groundstate%npsden+1, gc*rmt(is), ef_jlgr(:, ig, is) )
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          t1 = dot_product( v, atposc(:, ia, is) )
          ef_sfacg(ig, ias) = cmplx( cos( t1 ), sin( t1 ), dp )
        end do
      end do
    end do
!$omp end do
!$omp end parallel
  end subroutine ef_var_init

  subroutine ef_var_free
    if( allocated( ef_jlgr ) ) deallocate( ef_jlgr )
    if( allocated( ef_ylmg ) ) deallocate( ef_ylmg )
    if( allocated( ef_sfacg ) ) deallocate( ef_sfacg )
  end subroutine ef_var_free

end module efield_variables
