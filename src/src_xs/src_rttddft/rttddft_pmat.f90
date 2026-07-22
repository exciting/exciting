!> Module for the momentum matrix in RT-TDDFT
module rttddft_pmat
#include "asserts.fpp"
  use constants, only: real_zero, zi, zone, zzero
  use exciting_mpi, only: mpiinfo
  use math_utils, only: integer_sqrt
  use mod_APW_LO, only: apword, lorbl, nlorb, nlotot, nlomax, lolmax
  use mod_atoms, only: idxas, natoms, natmtot, nspecies
  use mod_eigensystem, only: idxlo, nmat
  use mod_gkvector, only: gkc, igkig, ngk, vgkc
  use mod_gvector, only: cfunig, ivg, ivgig
  use modmpi, only: terminate_if_false
  use modxs, only: apwcmt, locmt, ripaa, ripalo, riploa, riplolo
  use precision, only: dp, i32
  use rttddft_arrays, only: generic_matrix_set, hermitian_matrix_set
  use rttddft_file_formats, only: file_handler
  use rttddft_io_unformatted, only: file_pmat_exists, file_pmat_mt_exists, &
    get_filename_pmat, get_filename_pmat_mt, read_pmat, read_pmat_mt, write_pmat, &
    write_pmat_mt
  use rttddft_Overlap, only: overlap_set
  use xlapack, only: matrix_multiply

  implicit none
  private 
  integer(i32), parameter :: n_cartesian = 3
  real(dp), parameter :: eps_scissor = 1.e-10_dp
  
  !> Type to encapsulate the x, y, and z components of the momentum matrix
  type, public :: pmat_set
  private
    !> `x`, `y`, and `z` component of the momentum matrix
    class(generic_matrix_set), allocatable, public :: components(:)
    !> MT part of the momentum matrix (x, y, and z components)
    complex(dp), allocatable, public :: MT(:, :, :, :, :)
    !> If `.true.`, the momentum matrix is built in the LAPW+lo basis set
    logical, private :: lapwlo_basis = .true.
  contains
  private
    procedure, private :: assert_all_allocated => pmat_set_assert_all_allocated
    procedure, private :: assert_consistency => pmat_set_assert_consistency
    procedure, public :: allocate => pmat_set_allocate
    procedure, private :: calculate_in_KSBasis => pmat_set_calculate_in_KSBasis
    procedure, private :: calculate_in_LAPWBasis => pmat_set_calculate_in_LAPWBasis
    procedure, public :: calculate => pmat_set_calculate
    procedure, public :: is_hermitian => pmat_set_is_hermitian !@Ronaldo I am not sure about this one. Should we 
    ! consider pmat to be non-hermitian when calculating current density with MD?
    procedure, public :: read_from_file => pmat_set_read
    procedure, public :: write_to_file => pmat_set_write
    procedure, public :: represented_in_lapwlo => pmat_set_represented_in_lapwlo
    procedure, public :: scale => pmat_scale
    final :: destructor
  end type

contains
  !> Allocate the set of momentum matrices
  subroutine pmat_set_allocate( pmat, m, ki, kf, is_hermitian, allocate_MT, is_LAPWLO_basis )
    class(pmat_set), intent(inout) :: pmat
    !> Dimension of the matrices for all \( \mathbf{k} \)-points
    integer(i32), intent(in) :: m
    !> First \( \mathbf{k} \)-point managed by MPI process
    integer(i32), intent(in) :: ki
    !> Last \( \mathbf{k} \)-point managed by MPI process
    integer(i32), intent(in) :: kf
    !> If .True., for each `ik`, force the x, y, and z components of `pmat` to be hermitian
    logical, intent(in) :: is_hermitian
    !> If .True., allocate the MT part of the momentum matrix
    logical, intent(in) :: allocate_MT
    !> if `.True`, LAPW+lo basis is used
    logical, intent(in) :: is_LAPWLO_basis

    integer(i32) :: i

    if ( allocate_MT ) then
      CALL_ASSERT( is_LAPWLO_basis, "LAPWLO basis must be used if MT momentum matrices are required")
    end if
    if( is_hermitian ) then
      allocate( hermitian_matrix_set :: pmat%components(n_cartesian) )
    else
      allocate( generic_matrix_set :: pmat%components(n_cartesian) )
    end if
    do i = 1, n_cartesian
      call pmat%components(i)%allocate_array( [1, 1, ki], [m, m, kf] )
    end do
    if ( allocate_MT ) allocate( pmat%MT(m, m, n_cartesian, natmtot, ki:kf) )
    pmat%lapwlo_basis = is_LAPWLO_basis
  end subroutine

  pure elemental logical function pmat_set_is_hermitian( pmat ) result( is_hermitian )
    class(pmat_set), intent(in) :: pmat

    is_hermitian = .false.
    if( allocated( pmat%components ) ) then
      associate( p_x => pmat%components(1) )
        select type( p_x )
          type is( hermitian_matrix_set )
            is_hermitian = .true.
        end select
      end associate
    end if
  end function

  impure elemental subroutine destructor( this )
    type(pmat_set), intent(inout) :: this

    integer(i32) :: i

    if( allocated( this%components ) ) then
      do i = 1, n_cartesian
        call this%components(i)%deallocate_if_allocated()
      end do
      deallocate( this%components )
    end if
    if( allocated( this%MT ) ) deallocate( this%MT )
  end subroutine

  !> Assert that all components of the momentum matrix are allocated
  subroutine pmat_set_assert_all_allocated( this )
    class(pmat_set), intent(in) :: this

    integer(i32) :: j

    CALL_ASSERT( allocated( this%components ), "pmat not allocated" )
    do j = 1, n_cartesian
      call this%components(j)%assert_allocated()
    end do
  end subroutine

  subroutine pmat_set_assert_consistency( this )
    class(pmat_set), intent(in) :: this

    call this%assert_all_allocated()
    CALL_ASSERT( all( shape( this%components(1)%array ) == shape( this%components(2)%array ) ), "pmat x and y: not consistent shapes" )
    CALL_ASSERT( all( shape( this%components(1)%array ) == shape( this%components(3)%array ) ), "pmat x and z: not consistent shapes" )
  end subroutine

  !> Wrapper for calling [[rttddft_io_unformatted(module):write_pmat]]
  subroutine pmat_set_write( this, mpi_env, handler, n_kpt )
    class(pmat_set), intent(inout) :: this
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(in) :: mpi_env
    !> File handler
    type(file_handler), optional, intent(in) :: handler
    !> Number of k-points
    integer(i32), optional, intent(in) :: n_kpt

    call this%assert_consistency()
    call write_pmat( lbound( this%components(1)%array, 3 ), this%components(1)%array, &
      this%components(2)%array, this%components(3)%array, mpi_env, handler, n_kpt )
    if( allocated( this%MT ) ) call write_pmat_MT( lbound( this%MT, 5 ), this%MT, mpi_env, handler, n_kpt )
  end subroutine

  !> Wrapper for calling [[rttddft_io_unformatted(module):read_pmat]]
  subroutine pmat_set_read( this, mpi_env, handler )
    class(pmat_set), intent(inout) :: this
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(in) :: mpi_env
    !> File handler
    type(file_handler), optional, intent(in) :: handler

    call this%assert_consistency()
    call terminate_if_false( file_pmat_exists( handler, mpi_env ), &
        'File:'//trim( get_filename_pmat() )//' not found')
    call read_pmat( lbound( this%components(1)%array, 3 ), this%components(1)%array, &
      this%components(2)%array, this%components(3)%array, mpi_env, handler )
    if( allocated( this%MT ) ) then
      call terminate_if_false( file_pmat_mt_exists( handler, mpi_env ), &
        'File:'//trim( get_filename_pmat_mt() )//' not found')
      call read_pmat_mt( lbound( this%MT, 5 ), this%MT, mpi_env, handler )
    end if
  end subroutine

  !> Choose between the LAPW+lo and KS basis for the calculation of the momentum matrix
  subroutine pmat_set_calculate( this, first_kpt, apwalm, psi_gnd_lapwlo, psi_gnd_second_variation )
    class(pmat_set), intent(inout) :: this
    !> The first k point
    integer(i32), intent(in) :: first_kpt
    !> Matching coefficients of the (L)APWs
    !> (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), intent(in) :: apwalm(:, :, :, :, first_kpt :)
    !> All initial KS states in the LAPW+lo basis (nmatmax, n_states, first_kpt : last_kpt)
    complex(dp), contiguous, optional, intent(in) :: psi_gnd_lapwlo(:, :, first_kpt :)
    !> Initial KS states: second-variational wavefunctions (n_states_sv, n_states_sv, first_kpt : last_kpt)
    complex(dp), contiguous, optional, intent(in) :: psi_gnd_second_variation(:, :, first_kpt:)

    if ( this%represented_in_lapwlo() ) then
      call pmat_set_calculate_in_LAPWBasis( this, first_kpt, apwalm )
    else
      CALL_ASSERT( present( psi_gnd_lapwlo ), "psi_gnd_lapwlo should be present if KS basis is used" )
      call pmat_set_calculate_in_KSBasis( this, first_kpt, apwalm, psi_gnd_lapwlo, psi_gnd_second_variation )
    end if
  end subroutine

  !> Calculate the momentum matrix elements in the unperturbed KS basis. This routine 
  !> is essentially a wrapper for `/src/src_xs/genpmatxs.f90`
  subroutine pmat_set_calculate_in_KSBasis( pmat, first_kpt, apwalm, psi_gnd_lapwlo, psi_gnd_second_variation )
    class(pmat_set), intent(inout) :: pmat
    !> The first k point
    integer(i32), intent(in) :: first_kpt
    !> Matching coefficients of the (L)APWs
    !> (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), intent(in) :: apwalm(:, :, :, :, first_kpt :)
    !> Initial KS states in the LAPW+lo basis (nmatmax, n_states, first_kpt : last_kpt)
    complex(dp), contiguous, intent(in) :: psi_gnd_lapwlo(:, :, first_kpt :)
    !> Initial KS states: second-variational wavefunctions (n_states_sv, n_states_sv, first_kpt : last_kpt)
    complex(dp), contiguous, optional, intent(in) :: psi_gnd_second_variation(:, :, first_kpt:)

    
    integer(i32) :: apwordmax, ik, j, jl, jr, last_kpt
    integer(i32) :: l_max_apw, lm_max_apw, n_atoms_tot, n_states_first_variation, n_states_second_variation, n_used_states
    complex(dp), allocatable :: pmat_tmp(:, :, :, :), fake_evecsv(:, :)

    apwordmax = size( apwalm, 2 )
    lm_max_apw = size( apwalm, 3 )
    n_atoms_tot = size( apwalm, 4 )
    last_kpt = ubound( apwalm, 5 )
    n_states_first_variation = size( psi_gnd_lapwlo, 2 )
    n_states_second_variation = n_states_first_variation
    if( present( psi_gnd_second_variation ) ) n_states_second_variation = size( psi_gnd_second_variation, 1 )

    l_max_apw = integer_sqrt( lm_max_apw ) - 1
    CALL_ASSERT( (l_max_apw + 1)**2 == lm_max_apw, "lm_max_apw is not a perfect square" )
    
    allocate( pmat_tmp(3, n_states_second_variation, n_states_second_variation, first_kpt : last_kpt), source = zzero )
    if( .not. present( psi_gnd_second_variation ) ) allocate( fake_evecsv(1, 1) )
    n_used_states = size( pmat%components(1)%array, 1 )
    ! we may not need pmat matrix elements for the lowest-lying states
    CALL_ASSERT( n_used_states <= n_states_second_variation, 'n_used_states > n_states_second_variation for pmat' )

    do j = 1, n_cartesian
      pmat%components(j)%array = zzero
    end do
    if ( allocated( apwcmt ) ) deallocate( apwcmt )
    allocate( apwcmt(n_states_first_variation, apwordmax, lm_max_apw, n_atoms_tot), source = zzero )
    if ( allocated( ripaa ) ) deallocate( ripaa )
    allocate( ripaa( apwordmax, lm_max_apw, apwordmax, lm_max_apw, n_atoms_tot, 3), source = real_zero )
    if( nlotot > 0 ) then
      if ( allocated ( locmt ) ) deallocate( locmt )
      allocate( locmt(n_states_first_variation, nlomax,-lolmax : lolmax, n_atoms_tot), source = zzero )
      if ( allocated ( ripalo ) ) deallocate( ripalo )
      allocate( ripalo(apwordmax, lm_max_apw, nlomax,-lolmax : lolmax, n_atoms_tot, 3), source = real_zero )
      if ( allocated ( riploa ) ) deallocate( riploa )
      allocate( riploa(nlomax,-lolmax : lolmax, apwordmax, lm_max_apw, n_atoms_tot, 3), source = real_zero )
      if ( allocated ( riplolo ) ) deallocate( riplolo )
      allocate( riplolo(nlomax,-lolmax : lolmax, nlomax,-lolmax : lolmax, n_atoms_tot, 3), source = real_zero )
    end if
    call pmatrad()

    do ik = first_kpt, last_kpt
      call genapwcmt( l_max_apw, ngk(1, ik), 1, n_states_first_variation, apwalm(:, :, :, :, ik), &
        psi_gnd_lapwlo(:, :, ik), apwcmt )
      if( nlotot > 0 ) call genlocmt( ngk(1, ik), 1, n_states_first_variation, psi_gnd_lapwlo(:, :, ik), &
        locmt )
      if( present( psi_gnd_second_variation ) ) then
        call genpmatxs( ngk(1, ik), igkig(:, 1, ik), vgkc(:, :, 1, ik), &
          psi_gnd_lapwlo(:, :, ik), psi_gnd_second_variation(:, :, ik), pmat_tmp(:, :, :, ik) )
      else
        call genpmatxs( ngk(1, ik), igkig(:, 1, ik), vgkc(:, :, 1, ik), &
          psi_gnd_lapwlo(:, :, ik), fake_evecsv, pmat_tmp(:, :, :, ik) )
      end if

      do j = 1, n_cartesian
        pmat%components(j)%array(:, :, ik) = pmat_tmp(j, :, :, ik)
      end do

      ! Forces the matrix to be hermitian
      if( pmat%is_hermitian() ) then
        do j = 1, n_cartesian
          do jl = 1, n_used_states
            pmat%components(j)%array(jl, jl, ik) = cmplx( real(pmat%components(j)%array(jl, jl, ik), dp), 0.0_dp, dp )
            do jr = jl + 1, n_used_states
              pmat%components(j)%array(jr, jl, ik) = conjg( pmat%components(j)%array(jl, jr, ik) )
            end do
          end do
        end do
      end if
    end do

    deallocate( apwcmt, ripaa )
    if ( nlotot > 0 ) deallocate( locmt, ripalo, riploa, riplolo )
    
  end subroutine

  !> Here, we calculate the momentum matrix elements considering as basis 
  !> (L)APW+lo. We copied most of the code from `/src/src_xs/genpmatxs.F90`, 
  !> but there the basis are the KS-wavefunctions
  subroutine pmat_set_calculate_in_LAPWBasis( pmat, first_kpt, apwalm )
    class(pmat_set), intent(inout) :: pmat
    !> The first k point
    integer(i32), intent(in) :: first_kpt
    !> Matching coefficients of the (L)APWs
    !> (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), intent(in) :: apwalm(:, :, :, :, first_kpt :)

    integer(i32) :: apwordmax, ik, j, last_kpt, lm_max_apw, n_atoms_tot

    do j = 1, n_cartesian
      pmat%components(j)%array = zzero
    end do
    if ( allocated( pmat%MT ) ) pmat%MT = zzero
    apwordmax = size( apwalm, 2 )
    lm_max_apw = size( apwalm, 3 )
    n_atoms_tot = size( apwalm, 4 )
    last_kpt = ubound( apwalm, 5 )

    if(allocated(ripaa)) deallocate(ripaa)
    allocate(ripaa(apwordmax, lm_max_apw, apwordmax, lm_max_apw, n_atoms_tot, 3))
    if(nlotot > 0) then
      if(allocated(ripalo)) deallocate(ripalo)
      allocate(ripalo(apwordmax, lm_max_apw, nlomax,-lolmax:lolmax, n_atoms_tot, 3))
      if(allocated(riploa)) deallocate(riploa)
      allocate(riploa(nlomax,-lolmax:lolmax, apwordmax, lm_max_apw, n_atoms_tot, 3))
      if(allocated(riplolo)) deallocate(riplolo)
      allocate(riplolo(nlomax,-lolmax:lolmax, nlomax,-lolmax:lolmax, n_atoms_tot, 3))
    end if

    ! Calculate gradient of radial functions times spherical harmonics
    call pmatrad
    
!$OMP PARALLEL DEFAULT(NONE), PRIVATE(ik) SHARED(pmat, first_kpt, last_kpt, apwalm)
!$OMP DO
    do ik = first_kpt, last_kpt
      call generate_pmat_ik( pmat, ik, apwalm(:, :, :, :, ik) )
    end do
!$OMP END DO
!$OMP END PARALLEL

    deallocate( ripaa )
    if(nlotot > 0) deallocate( ripalo, riploa, riplolo )
  end subroutine

  !> Obtain the momentum matrix elements for a given a `k-point`
  subroutine generate_pmat_ik(pmat, ik, apwalmk)
    class(pmat_set), intent(inout) :: pmat
    !> Index of the `k-point` considered
    integer, intent(in)       :: ik
    !> apwalmk:  The matching coefficients for the (L)APW's (coefficient that
    !>           smoothly connect LAPW's and their MT counterparts)
    !>           Dimensions: ngkmax,apwordmax,lmmaxapw,natmtot
    complex(dp), intent(in)   :: apwalmk(:, :, :, :)

    integer(i32) :: is, ia, ias, io, io1, io2, ilo, ilo1, ilo2
    integer(i32) :: ig1, igp1, ig2, igp2, ig, iv1 (3), iv (3)
    integer(i32) :: j, l1, m1, lm1, l2, m2, lm2, l_max_apw, lm_max_apw
    integer(i32) :: nmatp, ng
    complex(dp), allocatable :: zv (:), aux(:,:)

    nmatp = nmat(1,ik)
    ng = ngk(1,ik)

    allocate (zv(ng))
    allocate (aux(ng,ng))

    lm_max_apw = size( apwalmk, 3 )
    l_max_apw = integer_sqrt( lm_max_apw ) - 1
    CALL_ASSERT( (l_max_apw + 1)**2 == lm_max_apw, "lm_max_apw is not a perfect square" )

    ! loop over species and atoms
    Do is = 1, nspecies
      Do ia = 1, natoms (is)
        ias = idxas (ia, is)
        !---------------------------!
        !     APW-APW contribution  !
        !---------------------------!
        Do j = 1, n_cartesian ! (loop over the x,y,z components)
          aux = zzero
          lm1 = 0
          Do l1 = 0, l_max_apw
            Do m1 = - l1, l1
              lm1 = lm1 + 1
              Do io1 = 1, apword (l1, is)
                zv = zzero
                lm2 = 0
                Do l2 = 0, l_max_apw
                  Do m2 = - l2, l2
                    lm2 = lm2 + 1
                    Do io2 = 1, apword (l2, is)
                      zv(1:ng) = zv(1:ng) + ripaa(io1,lm1,io2,lm2,ias,j)*apwalmk(1:ng,io2,lm2,ias)
                    End Do
                  End Do
                End Do
                Call zoutpr(ng,ng,zone,apwalmk(1:ng,io1,lm1,ias),zv(1:ng),aux(1:ng,1:ng))
              End Do
            End Do
          End Do
          pmat%components(j)%array(1:ng, 1:ng, ik) = pmat%components(j)%array(1:ng, 1:ng, ik) + & 
            aux(1:ng, 1:ng)
          if ( allocated( pmat%MT ) ) pmat%MT(1:ng, 1:ng, j, ias, ik) = aux(1:ng, 1:ng)
        End Do !Do j = 1, n_cartesian
        if (nlotot > 0) then
        !--------------------------------------!
        !     APW-local-orbital contribution   !
        !--------------------------------------!
          Do j = 1, n_cartesian
            lm1 = 0
            Do l1 = 0, l_max_apw ! (loop over lapw)
              Do m1 = - l1, l1
                lm1 = lm1 + 1
                Do io = 1, apword(l1, is)
                  Do ilo = 1, nlorb(is) ! (loop over local orbitals)
                    l2 = lorbl(ilo, is)
                    lm2 = l2**2
                    Do m2 = - l2, l2
                      lm2 = lm2 + 1
                      zv(1:ng) = ripalo(io,lm1,ilo,m2,ias,j)*conjg(apwalmk(1:ng,io,lm1,ias))
                      pmat%components(j)%array(1:ng, ng+idxlo(lm2, ilo, ias), ik) = &
                        pmat%components(j)%array(1:ng, ng+idxlo(lm2, ilo, ias), ik) + &
                        zv(1:ng)
                      if( allocated( pmat%MT ) ) pmat%MT(1:ng, ng+idxlo(lm2, ilo, ias), j, ias, ik) = &
                          & pmat%MT(1:ng, ng+idxlo(lm2,ilo,ias), j, ias, ik) + zv(1:ng)
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
          !--------------------------------------!
          !     local-orbital-APW contribution   !
          !--------------------------------------!
          Do j = 1, n_cartesian
            Do ilo = 1, nlorb (is)
              l1 = lorbl(ilo, is)
              lm1 = l1**2
              Do m1 = - l1, l1
                lm1 = lm1 + 1
                lm2 = 0
                Do l2 = 0, l_max_apw
                  Do m2 = - l2, l2
                    lm2 = lm2 + 1
                    Do io = 1, apword(l2, is)
                      zv(1:ng) = riploa(ilo,m1,io,lm2,ias,j)*apwalmk(1:ng,io,lm2,ias)
                      pmat%components(j)%array(ng+idxlo(lm1,ilo,ias), 1:ng, ik) = &
                        pmat%components(j)%array(ng+idxlo(lm1,ilo,ias), 1:ng, ik) + zv(1:ng)
                      if ( allocated( pmat%MT ) ) pmat%MT(ng+idxlo(lm1,ilo,ias),1:ng,j,ias,ik) = &
                          & pmat%MT(ng+idxlo(lm1,ilo,ias),1:ng,j,ias,ik) + zv(1:ng)
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
          !------------------------------------------------!
          !     local-orbital-local-orbital contribution   !
          !------------------------------------------------!
          Do j = 1, n_cartesian
            Do ilo1 = 1, nlorb (is)
              l1 = lorbl(ilo1, is)
              lm1 = l1**2
              Do m1 = - l1, l1
                lm1 = lm1 + 1
                Do ilo2 = 1, nlorb (is)
                  l2 = lorbl(ilo2, is)
                  lm2 = l2**2
                  Do m2 = - l2, l2
                    lm2 = lm2 + 1
                    pmat%components(j)%array(ng+idxlo(lm1,ilo1,ias),ng+idxlo(lm2,ilo2,ias),ik) = &
                          & riplolo(ilo1,m1,ilo2,m2,ias,j)
                    if ( allocated( pmat%MT ) ) pmat%MT(ng+idxlo(lm1,ilo1,ias),ng+idxlo(lm2,ilo2,ias), j, ias, ik) = &
                        & riplolo(ilo1,m1,ilo2,m2,ias,j)
                  End Do
                End Do
              End Do
            End Do
          End Do
        ! end case of local orbitals
        End If
      ! end loop over atoms and species
      End Do
    End Do

    ! x and z components
    pmat%components(1)%array(:, :, ik) = -zi*pmat%components(1)%array(:, :, ik)
    pmat%components(3)%array(:, :, ik) = -zi*pmat%components(3)%array(:, :, ik)
    if ( allocated( pmat%MT ) ) then
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          pmat%MT(:, :, 1, ias, ik) = -zi*pmat%MT(:, :, 1, ias, ik)
          pmat%MT(:, :, 3, ias, ik) = -zi*pmat%MT(:, :, 3, ias, ik)
        end do
      end do
    end if

    !  calculate momentum matrix elements in the interstitial region
    do j = 1, n_cartesian
      do igp1 = 1, ng
        ig1 = igkig(igp1,1,ik)
        iv1 (:) = ivg (:, ig1)
        do igp2 = 1, ng
          ig2 = igkig(igp2,1,ik)
          iv (:) = iv1 (:) - ivg (:, ig2)
          ig = ivgig (iv(1), iv(2), iv(3))
          pmat%components(j)%array(igp1, igp2, ik) = pmat%components(j)%array(igp1, igp2, ik) + vgkc(j, igp2, 1, ik)*cfunig(ig)
        end do
      end do
    end do

    if( pmat%is_hermitian() ) then
      do j = 1, n_cartesian
        do ig1 = 1, nmatp
          pmat%components(j)%array(ig1, ig1, ik) = cmplx( real(pmat%components(j)%array(ig1, ig1, ik), dp), 0.0_dp, dp )
          do ig2 = ig1+1, nmatp
            pmat%components(j)%array(ig2, ig1, ik) = conjg( pmat%components(j)%array(ig1, ig2, ik) )
          end do
        end do
      end do
    end if

  end subroutine

  !> Scaled momentum matrix to account for the scissor shift correction \( \Delta E \).
  !> A scaling matrix $f$ in the Kohn-Sham basis is [see PRB 48, 11789 (1993)]
  !> \[
  !> f_{ij} = 1 + \frac{\Delta E}{|\epsilon_i - \epsilon_j| - \Delta E},
  !> \]
  !> where \( i \) and \( j \) lie in different (conduction or valence) bands,
  !> and \( \epsilon_i \) is a quasiparticel (scissor-shifted) energy.
  !> The momentum matrix is updated via an element-wise (Hadamard) product:
  !> \[
  !> p'_{KS} = f \circ p_{KS}.
  !> \]
  !> When LAPW+lo basis is used, the momentum matrix is first transformed into the KS basis, 
  !> then scaled using the element-wise product, and transformed back to the LAPW+lo basis.
  subroutine pmat_scale( this, ks_lapwlo_transition_matrix, initial_eigenvalues, overlap, scissor_shift, n_occupied )
    !> Type that encapsulates the momentum matrix (to be allocated in `array_allocation` block)
    class(pmat_set), intent(inout) :: this
    !> KS-LAPW+lo transition matrix (n_lapwlo, n_states, first_kpt : last_kpt)
    complex(dp), intent(in):: ks_lapwlo_transition_matrix(:, :, :)
    !> Scissor-shifted initial KS eigenvalues (n_states, first_kpt : last_kpt)
    real(dp), intent(in) :: initial_eigenvalues(:, :)
    !> Object that encapsulates the overlap matrix \( S \)
    class(overlap_set), intent(in) :: overlap
    !> Scissor correction
    real(dp), intent(in) :: scissor_shift
    !> Number of occupied states
    integer(i32), intent(in) :: n_occupied

    integer(i32) :: i, j, ik, n_ks_states, first_kpt, last_kpt, k_shift, n_lapwlo_functions
    complex(dp), allocatable :: scaling_matrix_ks_basis(:, :), aux(:, :), &
      overlap_times_psi_lapwlo(:, :), pmat_ks(:, :)

    if ( scissor_shift > eps_scissor ) then
      n_ks_states = size( initial_eigenvalues, 1 )
      n_lapwlo_functions = size( ks_lapwlo_transition_matrix, 1 )
      first_kpt = lbound( overlap%array, 3 )
      last_kpt = ubound( overlap%array, 3 )

      allocate( scaling_matrix_ks_basis(n_ks_states, n_ks_states), source = zone )
      if ( this%represented_in_lapwlo() ) then
        allocate( overlap_times_psi_lapwlo(n_lapwlo_functions, n_ks_states), source = zzero )
        allocate( aux, source = overlap_times_psi_lapwlo )
        allocate( pmat_ks(n_ks_states, n_ks_states), source = zzero )
      end if

      k_shift = 1 - first_kpt
      do ik = first_kpt, last_kpt
        scaling_matrix_ks_basis = zone
        do i = 1, n_occupied
          do j = n_occupied + 1, n_ks_states
            scaling_matrix_ks_basis(i, j) = zone + scissor_shift / ( abs( initial_eigenvalues(i, k_shift + ik) - &
              initial_eigenvalues(j, k_shift + ik) ) - scissor_shift )
            scaling_matrix_ks_basis(j, i) = scaling_matrix_ks_basis(i, j)
          end do
        end do

        if ( this%represented_in_lapwlo() ) then
          call matrix_multiply( overlap%array(:, :, ik), &
            ks_lapwlo_transition_matrix(:, :, ik + k_shift), overlap_times_psi_lapwlo )
          do j = 1, n_cartesian
            ! transform to KS basis first
            call matrix_multiply( this%components(j)%array(:, :, ik), ks_lapwlo_transition_matrix(:, :, ik + k_shift), aux )
            call matrix_multiply( ks_lapwlo_transition_matrix(:, :, ik + k_shift), aux, pmat_ks, trans_A = 'C' )
            ! elementwise product
            pmat_ks = scaling_matrix_ks_basis * pmat_ks
            ! transform back to LAPWlo basis
            call matrix_multiply( overlap_times_psi_lapwlo, pmat_ks, aux )
            call matrix_multiply( aux, overlap_times_psi_lapwlo, &
              this%components(j)%array(:, :, ik), trans_B = 'C' )
          end do
        else
          ! elementwise product
          do j = 1, n_cartesian
            this%components(j)%array(:, :, ik) = scaling_matrix_ks_basis * &
              this%components(j)%array(:, :, ik)
          end do
        end if
      end do
    end if
  end subroutine

  !> Return whether the momentum matrix is represented in the LAPW+lo basis
  pure logical function pmat_set_represented_in_lapwlo( this ) result( represented_in_lapwlo )
    class(pmat_set), intent(in) :: this

    represented_in_lapwlo = this%lapwlo_basis
  end function

end module rttddft_pmat
