!> This module contains the global variables for the RT-TDDFT implementation. 
!> For details, refer to the reference with DOI: 10.1088/2516-1075/ac7afc
module rttddft_GlobalMDVariables
  use mod_gkvector, only: ngk, ngkmax, gkc, tpgkc, sfacgk, vgkc
  use mod_gvector, only: ngvec, vgc, sfacg
  use mod_spin, only: nspnfv
  use precision, only: dp, i32
  use mod_gen_lo, only: genlofr

  implicit none

  private
  ! variables
  public :: mathcalH, mathcalB, B_time, B_past
  ! subroutines
  public :: update_exciting_globals_for_new_ions_positions

  !> `mathcalH` gives the impact of an ion displacement on the hamiltonian matrix
  !> \[ \left[ \left\langle 
  !> \frac{\partial \phi_{\mu'}^{\mathbf{k}}}{\partial \mathbf{R}_J}
  !> \Bigg|\hat{H}\Bigg|\phi_{\mu}^{\mathbf{k}}\right\rangle +
  !> \left\langle\phi_{\mu'}^{\mathbf{k}}\Bigg|\hat{H}\Bigg|\frac{\partial 
  !> \phi_{\mu}^{\mathbf{k}}}{\partial \mathbf{R}_J}\right\rangle \right] 
  !> \]
  complex(dp), allocatable  :: mathcalH(:,:,:,:,:)
  
  !> `mathcalB` measures how the ions displacements affect overlap elements
  !> \[ \mathcal{B}_{J\mu'\mu}^{\mathbf{k}} = \left \langle
  !> \phi_{\mu'}^{\mathbf{k}}\bigg| \frac{\partial}{\partial \mathbf{R}_J}
  !> \phi_{\mu}^{\mathbf{k}} \right\rangle \]
  complex(dp), allocatable  :: mathcalB(:,:,:,:,:)
  
  !> `B_time` quantifies the impact of the Ehrenfest molecular dynamics
  !> on the time evolution of the electonic wavefunctions.
  !> \[B_{\mu'\mu}^{\mathbf{k}} = \left \langle
  !> \phi_{\mu'}^{\mathbf{k}}\left|\frac{d}{d t}\right.
  !> \phi_{\mu}^{\mathbf{k}}\right\rangle =
  !> \sum_J \dot{\mathbf{R}}_J\cdot \mathcal{B}_{J\mu'\mu}^{\mathbf{k}} \]
  complex(dp), allocatable  :: B_time(:,:,:)
  !> Same as `B_time`, but at the previous time step: \(t-\Delta t\)
  complex(dp), allocatable  :: B_past(:,:,:)

contains

subroutine update_exciting_globals_for_new_ions_positions( first_kpt, apwalm )
  !> index of the first k-point
  integer(i32), intent(in) :: first_kpt
  !> Matching coefficients of the (L)APWs
  !> (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
  complex(dp), contiguous, intent(inout) :: apwalm(:, :, :, :, first_kpt :)

  integer(i32) :: ik, i_spin, last_kpt

  last_kpt = ubound( apwalm, 5 )

  call checkmt     ! check for overlapping muffin-tins
  call gencfun     ! generate the characteristic function
  call energynn    ! determine the nuclear-nuclear energy
  ! generate structure factors for G and G+k-vectors
  call gensfacgp (ngvec, vgc, ngvec, sfacg)
  do ik = first_kpt, last_kpt
    do i_spin = 1, nspnfv
      call gensfacgp (ngk(i_spin, ik), vgkc(:, :, i_spin, ik), ngkmax, sfacgk(:, :, i_spin, ik))
    end do
  end do
  call gencore( )       ! generate the core wavefunctions and densities
  call linengy( )       ! find the new linearization energies
  call genapwfr( )      ! generate the APW radial functions
  call genlofr( )       ! generate the local-orbital radial functions
  call olprad( )
  ! Matching coefficients (apwalm)
  do ik = first_kpt, last_kpt
    call match( ngk(1,ik), gkc(:,1,ik), tpgkc(:,:,1,ik), sfacgk(:,:,1,ik), apwalm(:,:,:,:,ik) )
  end do
end subroutine

end module rttddft_GlobalMDVariables
