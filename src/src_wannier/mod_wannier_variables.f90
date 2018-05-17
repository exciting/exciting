module mod_wannier_variables
  use mod_kpointset

  implicit none

  type wannier_group
    integer :: fst, lst, nst, nwf, fwf, lwf
    integer :: nprojused
    integer, allocatable :: projused(:)
    character(32) :: method
    real(8) :: win_i(2), win_o(2)
    integer, allocatable :: win_ii(:,:), win_io(:,:)
    integer, allocatable :: win_ni(:), win_no(:)
    complex(8), allocatable :: subspace(:,:,:)
  end type wannier_group

! module variables
  integer :: wf_nprojtot                                ! total number of projection functions
  integer :: wf_fst                                     ! lowest band used for generation
  integer :: wf_lst                                     ! highest band used for generation
  integer :: wf_nst                                     ! (maximum) number of bands used for generation
  integer :: wf_nwf                                     ! number of Wannier functions to generate
  integer :: wf_info                                    ! unit to write WANNIER_INFO.OUT
  logical :: wf_initialized = .false.
  real(8) :: wf_t0
  type( k_set) :: wf_kset
  type( k_set) :: wf_kset_red
  type( G_set) :: wf_Gset
  type( Gk_set) :: wf_Gkset
  character(265) :: wf_filename
  logical :: wf_fermizero

  type( wannier_group), allocatable :: wf_groups(:)
  integer :: wf_ngroups
  integer :: wf_group
  
  integer, allocatable :: wf_evecphase(:,:)             ! indices for phase correction of eigenvectors
  integer, allocatable :: wf_projst(:,:)                ! l, m, species, atom of local-orbitals for projection
  complex(8), allocatable :: wf_transform(:,:,:)        ! unitary transformation matrices
  complex(8), allocatable :: wf_projection(:,:,:)       ! overlap of wafefunctions and local-orbitals
  complex(8), allocatable :: wf_opf(:,:)                ! expansion coefficients for optimal projection functions
  complex(8), allocatable :: wf_m0(:,:,:,:)             ! original plane-wave matrix elements for neighboring k-points
  complex(8), allocatable :: wf_m(:,:,:,:)              ! transformed plane-wave matrix elements for neighboring k-points
  real(8), allocatable :: wf_centers(:,:)               ! centers of Wannier functions
  real(8), allocatable :: wf_omega(:)                   ! total localization functional for each Wannier function
  real(8), allocatable :: wf_omega_i(:)                 ! gauge independent localization functional for each Wannier function
  real(8), allocatable :: wf_omega_d(:)                 ! diagonal localization functional for each Wannier function
  real(8), allocatable :: wf_omega_od(:)                ! off-diagonal localization functional for each Wannier function
  integer, allocatable :: wf_sheet(:,:)

  ! disentanglement
  logical :: wf_disentangle = .false.
  real(8) :: wf_win_i(2), wf_win_o(2)                   ! boundaries of inner/outer window
  integer, allocatable :: wf_win_ii(:,:), wf_win_io(:,:)! bands inside inner/outer window per k-point
  integer, allocatable :: wf_win_ni(:), wf_win_no(:)    ! number of bands in inner/outer window per k-point

  ! geometry
  integer :: wf_n_nshells                               ! number of shells used for gradient calculation
  integer :: wf_n_ntot                                  ! total numbers of neighbors in used shells
  real(8), allocatable :: wf_n_dist(:)                  ! distance of neighbors
  integer, allocatable :: wf_n_ik(:,:)                  ! k-point index of given neighbor for given k-point
  integer, allocatable :: wf_n_ik2(:,:)                 ! k-point index of given neighbor for given k-point
  integer, allocatable :: wf_n_usedshells(:)            ! index of used shells
  real(8), allocatable :: wf_n_vl(:,:), wf_n_vc(:,:)    ! vector of given neighbor in lattice and cartesian coordinates
  real(8), allocatable :: wf_n_wgts(:)                  ! geometric weight of given shell
  real(8), allocatable :: wf_n_wgt(:)                   ! geometric weight of given neighbor
  real(8) :: wf_n_rot(3,3)                              ! rotation to xy-plane for 2D grids

! methods
  contains
    subroutine wannier_destroy
      if( allocated( wf_projst)) deallocate( wf_projst)
      if( allocated( wf_transform)) deallocate( wf_transform)
      if( allocated( wf_projection)) deallocate( wf_projection)
      if( allocated( wf_opf)) deallocate( wf_opf)
      if( allocated( wf_m0)) deallocate( wf_m0)
      if( allocated( wf_m)) deallocate( wf_m)
      if( allocated( wf_centers)) deallocate( wf_centers)
      if( allocated( wf_omega)) deallocate( wf_omega)
      if( allocated( wf_omega_i)) deallocate( wf_omega_i)
      if( allocated( wf_omega_d)) deallocate( wf_omega_d)
      if( allocated( wf_omega_od)) deallocate( wf_omega_od)
      if( allocated( wf_n_dist)) deallocate( wf_n_dist)
      if( allocated( wf_n_ik)) deallocate( wf_n_ik)
      if( allocated( wf_n_ik2)) deallocate( wf_n_ik2)
      if( allocated( wf_n_vl)) deallocate( wf_n_vl)
      if( allocated( wf_n_vc)) deallocate( wf_n_vc)
      if( allocated( wf_n_wgt)) deallocate( wf_n_wgt)
      if( allocated( wf_n_wgts)) deallocate( wf_n_wgts)
      if( allocated( wf_n_usedshells)) deallocate( wf_n_usedshells)
      wf_initialized = .false.
    end subroutine wannier_destroy
end module mod_wannier_variables
