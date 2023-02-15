!> This module contains procedures for the calculation of the 
!> Hamiltonian and overlap matrix and their response which are
!> common for all types of perturbations.
!>
!> For a generic perturbation operator \(\delta\), the variation/response
!> of the overlap matrix \({\bf S}\) reads
!> \[\begin{align*}
!>   \delta S_{mn} &= \delta \langle \phi_m | \phi_n \rangle \\
!>                 &= \left( \delta \langle \phi_m| \right) |\phi_n \rangle + 
!>                    \langle \phi_m| \left( \delta |\phi_n \rangle \right) \\
!>                 &= \delta S^L_{mn} + \delta S^R_{mn} \;.
!> \end{align*}\]
!> The variation/response of the Hamiltonian matrix \({\bf H}\) reads
!> \[\begin{align*}
!>   \delta H_{mn} &= \delta \langle \phi_m | V_{\rm eff} + {\bf \hat{T}} | \phi_n \rangle \\
!>                 &= \left( \delta \langle \phi_m| \right) V_{\rm eff} + {\bf \hat{T}} |\phi_n \rangle + 
!>                    \langle \phi_m | \delta V_{\rm eff} | \phi_n \rangle + 
!>                    \langle \phi_m| V_{\rm eff} + {\bf \hat{T}} \left( \delta |\phi_n \rangle \right) \\
!>                 &= \delta H^L_{mn} + \delta H^0_{mn} + \delta H^R_{mn} \;,
!> \end{align*}\]
!> where we assumed that the kinetic energy operator \({\bf \hat{T}}\) is invariant
!> under the perturbation \(\delta\).
module dfpt_eigensystem
  use dfpt_variables

  use precision, only: dp
  use matrix_elements
  use asserts, only: assert
  use modmpi, only: terminate_if_false

  implicit none
  private

  !> radial muffin-tin integrals times Gaunt coefficients for overlap matrix \(S\)
  complex(dp), allocatable, public :: Smat_mt_basis(:,:,:)
  !> radial muffin-tin integrals times Gaunt coefficients for Hamiltonian matrix \(H\)
  complex(dp), allocatable, public :: Hmat_mt_basis(:,:,:)
  !> characteristic function in reciprocal space
  complex(dp), allocatable :: cfun_ig(:)
  !> interstitial effective potential times characteristic function in reciprocal space
  complex(dp), allocatable :: pot_cfun_ig(:)
  !> interstitial (scalar relativistic) kinetic energy times characteristic function in reciprocal space
  complex(dp), allocatable :: kin_cfun_ig(:)

  public :: dfpt_eig_init, dfpt_eig_free
  public :: dfpt_eig_ks, dfpt_eig_geteval, dfpt_eig_getevec
  public :: dfpt_eig_gen_dHmat ,dfpt_eig_prepare_dHmat

  contains

    !> This subroutine initializes variables for the calculation of matrix elements
    !> that remain constant during the entire calculation.
    !>
    !> This includes: 
    !> 
    !> * generation of the muffin-tin basis functions
    !> * generation of intertistial characteristic function
    !> * initialization of the matrix elements module
    !> * calculation of radial muffin-tin integrals times Gaunt coefficients 
    !>   for Hamiltonian and overlap
    !> * calculation of reciprocal space representation of interstitial potential 
    !>   times characteristic function
    subroutine dfpt_eig_init
      use constants, only: zzero, zone
      use physical_constants, only: alpha
      use mod_potential_and_density, only: pot_mt => veffmt, pot_ir => veffir
      use modinput

      real(dp), allocatable :: kin_ir(:)

      ! initialize matrix elements module
      call me_init( mt_basis, dfpt_lmaxvr, dfpt_Gset )

      ! compute MT radial integrals times Gaunt coefficients
      call gen_overlap_hamiltonian_mt_basis( Smat_mt_basis, Hmat_mt_basis, pot_mt, lmax_apw=dfpt_lmaxapw, lmax_pot=dfpt_lmaxvr )

      ! compute interstitial representation of characteristic function, 
      ! effective potential and kinetic energy
      call me_ir_alloc( pot_cfun_ig )
      call me_ir_prepare( zone, pot_ir, zzero, pot_cfun_ig )
      allocate( cfun_ig(dfpt_Gset%ngvec) )
      allocate( kin_ir(dfpt_Gset%ngrtot) )
      call gencfunig( dfpt_Gset%ngvec, dfpt_Gset%gc, dfpt_Gset%vgc, cfun_ig )
      select case( input%groundstate%ValenceRelativity )
        case( 'iora*' )
          call terminate_if_false( .false., '(dfpt_eig_init) &
            DFPT in combination with `iora*` valence relativity not implemented.' )
        case( 'none' )
          kin_ir = 0.5_dp
        case default
          kin_ir = 0.5_dp / (1.0_dp - 0.5_dp * alpha**2 * pot_ir )
      end select
      call me_ir_alloc( kin_cfun_ig )
      call me_ir_prepare( zone, kin_ir, zzero, kin_cfun_ig )
    end subroutine dfpt_eig_init

    !> This subroutine frees memory from the module variables
    !> and cleans up the matrix elements module.
    subroutine dfpt_eig_free
      if( allocated( Smat_mt_basis ) ) deallocate( Smat_mt_basis )
      if( allocated( Hmat_mt_basis ) ) deallocate( Hmat_mt_basis )
      if( allocated( pot_cfun_ig ) ) deallocate( pot_cfun_ig )
      if( allocated( kin_cfun_ig ) ) deallocate( kin_cfun_ig )
    end subroutine dfpt_eig_free

    !> This subroutine solves the Kohn-Sham equations non self-consistently 
    !> for the given wavevector \({\bf p}\) and returns the respective
    !> eigenvalues and eigenvectors.
    !>
    !> If files with eigenvectors and eigenvalues are provided, it is tried
    !> to find the result in the file. If it could not be found, eigenvalues
    !> and eigenvectors are computed.
    subroutine dfpt_eig_ks( ip, pset, Gset, Gpset, nst, eval, evec, &
        p0set, Gp0set, fevec, feval )
      use constants, only: zzero, zone
      use m_linalg, only: zhegdiag, zhegauge
      use block_data_file, only: block_data_file_type
      use mod_APW_LO, only: nlotot
      use mod_kpointset, only: k_set, G_set, Gk_set
      use mod_Gkvector, only: ngkmax_ptr
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use mod_APW_LO, only: apwordmax
      use mod_muffin_tin, only: lmmaxapw
      !> index of the wavevector \({\bf p}\) in the set
      integer, intent(in) :: ip
      !> set of \({\bf p}\) vectors
      type(k_set), intent(in) :: pset
      !> set of \({\bf G}\) vectors
      type(G_set), intent(in) :: Gset
      !> set of \({\bf G+p}\) vectors
      type(Gk_set), intent(in) :: Gpset
      !> number of states to solve for (`nst` states with lowest eigenvalues)
      integer, intent(in) :: nst
      !> eigenvalues
      real(dp), intent(out) :: eval(nst)
      !> eigenvectors
      complex(dp), intent(out) :: evec(Gpset%ngkmax+nlotot, nst)
      !> reference set of \({\bf p}\) vectors the files correspond to
      type(k_set), optional, intent(in) :: p0set
      !> \({\bf G+p}\) vectors for the reference set of \({\bf p}\) vectors
      type(Gk_set), optional, intent(in) :: Gp0set
      !> reference files with eigenvectors and eigenvalues
      type(block_data_file_type), optional, intent(inout) :: fevec, feval

      integer :: ip0, isym, i
      integer :: ngp, nmatp, nmatmax, n
      integer, target :: ngpmax
      integer :: is, ia, ias

      complex(dp), allocatable :: apwalm(:,:,:,:), S(:,:), H(:,:)
      real(dp), allocatable :: rbuff(:)
      complex(dp), allocatable :: zbuff(:,:)

      eval = 0.0_dp
      evec = zzero

      ! set matrix sizes
      ngp = Gpset%ngk(1, ip)
      nmatp = ngp + nlotot
      ngpmax = Gpset%ngkmax
      nmatmax = ngpmax + nlotot
      ngkmax_ptr => ngpmax
      n = min( nst, nmatp )

      ! read result from file if possible
      if( present( feval ) .or. present( fevec ) ) then
        ! check input
        call assert( present( p0set ), &
          'When reference files are present so must be the reference p-point set.' )
        call assert( present( p0set ), &
          'When reference files are present so must be the reference G+p-point set.' )
        call assert( present( feval ) .and. present( fevec ), &
          'Both reference files must be present.' )
        ! check if requested p-point is in reference set
        call findkptinset( pset%vkl(:, ip), p0set, isym, ip0 )
        ! p-point is in reference set
        if( ip0 > 0 ) then
          call dfpt_eig_geteval( pset%vkl(:, ip), feval, p0set, [1, n], rbuff )
          call dfpt_eig_getevec( pset%vkl(:, ip), Gpset%vgkl(:, :, 1, ip), fevec, p0set, Gp0set, [1, n], zbuff )
          eval(1:n) = rbuff
          evec(:, 1:n) = zbuff
          if( allocated( rbuff ) ) deallocate( rbuff )
          if( allocated( zbuff ) ) deallocate( zbuff )
          return
        end if
      endif
     
      ! ** solve KS equation
      ! allocate local variables
      allocate( apwalm(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) )
      allocate( S(nmatp, nmatp), source=zzero )
      allocate( H(nmatp, nmatp), source=zzero )
      ! * set up overlap and Hamiltonian matrix
      ! muffin-tin contribution
      call match( ngp, Gpset%gkc(:, 1, ip), Gpset%tpgkc(:, :, 1, ip), Gpset%sfacgk(:, :, 1, ip), apwalm )
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          call me_mt_mat( is, ias, ngp, ngp, apwalm(:, :, :, ias), apwalm(:, :, :, ias), zone, Smat_mt_basis(:, :, ias), zone, S )
          call me_mt_mat( is, ias, ngp, ngp, apwalm(:, :, :, ias), apwalm(:, :, :, ias), zone, Hmat_mt_basis(:, :, ias), zone, H )
        end do
      end do
      ! interstitial contribution
      call me_ir_mat( Gpset, ip, Gpset, ip, zone, cfun_ig, zone, S, Gset_op=Gset )
      call me_ir_mat( Gpset, ip, Gpset, ip, zone, pot_cfun_ig, zone, H, Gset_op=Gset )
      do i = 1, 3
        call me_ir_mat( Gpset, ip, Gpset, ip, zone, kin_cfun_ig, zone, H, Gset_op=Gset, &
          left_gradient=i, right_gradient=i )
      end do
      !* solve eigensystem
      call zhegdiag( H, S, eval(1:n), evec=evec(1:nmatp, 1:n), irange=[1,n] )
      ! we fix a unige gauge of the eigenvectors
      call zhegauge( eval(1:n), evec(1:nmatp, 1:n), eps=1e-6_dp )
      ! deallocate local variables
      deallocate( apwalm, S, H )
    end subroutine dfpt_eig_ks

    !> For a given effective potential response \(\delta V_{\rm eff}({\bf r})\), 
    !> this subroutine computes the radial muffin-tin integrals times Gaunt coefficients
    !> and the reciprocal space representation of the effective potential response
    !> times the characteristic function.
    subroutine dfpt_eig_prepare_dHmat( pot_mt, pot_ir, dpot_mt, dpot_ir, dHmat_mt_basis, dpot_cfun_ig, dkin_cfun_ig, &
        Gset )
      use matrix_elements
      use constants, only: zzero, zone, y00
      use physical_constants, only: alpha
      use mod_kpointset, only: G_set
      use mod_atoms, only: nspecies, natoms, idxas
      use mod_muffin_tin, only: nrmtmax, nrmt
      use modinput
      !> muffin-tin effective potential as real spherical harmonics expansion
      real(dp), intent(in) :: pot_mt(:,:,:)
      !> interstitial effective potential on real space FFT grid
      real(dp), intent(in) :: pot_ir(:)
      !> muffin-tin effective potential response as complex spherical harmonics expansion
      complex(dp), intent(in) :: dpot_mt(:,:,:)
      !> interstitial effective potential response on real space FFT grid
      complex(dp), intent(in) :: dpot_ir(:)
      !> radial muffin-tin integrals times Gaunt coefficients
      complex(dp), intent(out) :: dHmat_mt_basis(:,:,:)
      !> interstitial effective potential response times characteristic function in reciprocal space
      complex(dp), intent(out) :: dpot_cfun_ig(:)
      !> interstitial (scalar relativistic) kinetic energy response times characteristic function in reciprocal space
      complex(dp), intent(out) :: dkin_cfun_ig(:)
      !> set of \({\bf G}\) vectors the interstitial potential response is expanded on (default: `dfpt_Gset`)
      type(G_set), optional, intent(in) :: Gset

      integer :: is, ia, ias, i

      complex(dp), allocatable :: rfun(:,:), dkin_ir(:)

      allocate( rfun(1, nrmtmax) )
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! potential
          call me_mt_prepare( is, ias, dfpt_lmaxvr, zone, dpot_mt(:, :, ias), zzero, dHmat_mt_basis(:, :, ias) )
          ! kinetic energy
          if( input%groundstate%ValenceRelativity /= 'none' ) then
            rfun(1, 1:nrmt(is)) = 0.5_dp * (0.5_dp * alpha**2 * dpot_mt(1, 1:nrmt(is), ias)) / (1.0_dp - 0.5_dp * alpha**2 * pot_mt(1, 1:nrmt(is), ias) * y00)**2
            do i = 1, 3
              call me_mt_prepare( is, ias, 0, zone, rfun, zone, dHmat_mt_basis(:, :, ias), &
                left_gradient=i, right_gradient=i )
            end do
          end if

        end do
      end do
      deallocate( rfun )

      if( present( Gset ) ) then
        ! potential
        call me_ir_prepare( zone, dpot_ir, zzero, dpot_cfun_ig, Gset_op=Gset )
        ! kinetic energy
        if( input%groundstate%ValenceRelativity == 'none' ) then
          dkin_cfun_ig = zzero
        else
          allocate( dkin_ir(Gset%ngrtot) )
          dkin_ir = 0.5_dp * (0.5_dp * alpha**2 * dpot_ir) / (1.0_dp - 0.5_dp * alpha**2 * pot_ir)**2
          call me_ir_prepare( zone, dkin_ir, zzero, dkin_cfun_ig, Gset_op=Gset )
          deallocate( dkin_ir )
        end if
      else
        ! potential
        call me_ir_prepare( zone, dpot_ir, zzero, dpot_cfun_ig )
        ! kinetic energy
        if( input%groundstate%ValenceRelativity == 'none' ) then
          dkin_cfun_ig = zzero
        else
          allocate( dkin_ir(dfpt_Gset%ngrtot) )
          dkin_ir = 0.5_dp * (0.5_dp * alpha**2 * dpot_ir) / (1.0_dp - 0.5_dp * alpha**2 * pot_ir)**2
          call me_ir_prepare( zone, dkin_ir, zzero, dkin_cfun_ig )
          deallocate( dkin_ir )
        end if
      end if
    end subroutine dfpt_eig_prepare_dHmat

    !> This subroutine calculates the contribution to the Hamiltonian response 
    !> coming from the effective potential response, i.e.,
    !> \[ \delta H^0_{mn} = \langle \psi_{m{\bf p'}} | \delta V_{\rm eff} | \psi_{n{\bf p}} \rangle \;,\]
    !> where \({\bf p'}\) might be different from \({\bf p}\) when the perturbation carries a non-zero
    !> wavevector (e.g. phonon-like perturbation).
    !>
    !> @note The result is added to the input matrix! @endnote
    subroutine dfpt_eig_gen_dHmat( ip, Gpset1, Gpset2, fst1, lst1, fst2, lst2, evec1, evec2, apwalm1, apwalm2, &
        dHmat_mt_basis, dpot_cfun_ig, dkin_cfun_ig, dHmat, &
        Gset, diagonal )
      use constants, only: zone
      use mod_kpointset, only: G_set, Gk_set
      use mod_atoms, only: nspecies, natoms, idxas
      !> index of the wavevector \({\bf p}\) in the set
      integer, intent(in) :: ip
      !> set of \({\bf G+p}\) vectors on the left and right
      type(Gk_set), intent(in) :: Gpset1, Gpset2
      !> first and last state on the left for which the matrix elements are calculated
      integer, intent(in) :: fst1, lst1
      !> first and last state on the left for which the matrix elements are calculated
      integer, intent(in) :: fst2, lst2
      !> eigenvectors on the left and right
      complex(dp), intent(in) :: evec1(:,:), evec2(:,:)
      !> (L)APW matching coefficients \(A^\alpha_{{\bf G+p},lm,\xi}\) on the left and right
      complex(dp), intent(in) :: apwalm1(:,:,:,:), apwalm2(:,:,:,:)
      !> radial muffin-tin integrals of effective potential response times Gaunt coefficients
      complex(dp), intent(in) :: dHmat_mt_basis(:,:,:)
      !> interstitial effective potential response times characteristic function in reciprocal space
      complex(dp), intent(in) :: dpot_cfun_ig(:)
      !> interstitial (scalar relativistic) kinetic energy response times characteristic function in reciprocal space
      complex(dp), intent(in) :: dkin_cfun_ig(:)
      !> Hamiltonian response
      complex(dp), intent(inout) :: dHmat(:,:)
      !> set of \({\bf G}\) vectors the interstitial potential response is expanded on (default: dfpt_Gset)
      type(G_set), optional, intent(in) :: Gset
      !> if true, only `fst2` to `lst2` diagonal elements are computed
      !> and stored in the first dimension of the output (default: `.false.`)
      logical, optional, intent(in) :: diagonal

      integer :: ngk1, ngk2, nst1, nst2
      integer :: is, ia, ias, i
      logical :: diag

      diag = .false.
      if( present( diagonal ) ) diag = diagonal

      ngk1 = Gpset1%ngk(1, ip)
      ngk2 = Gpset2%ngk(1, ip)
      nst1 = lst1 - fst1 + 1
      nst2 = lst2 - fst2 + 1

      ! ** muffin-tin part
      ! calculate muffin-tin integrals
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          call me_mt_mat( is, ias, ngk1, ngk2, apwalm1(:, :, :, ias), apwalm2(:, :, :, ias), &
                 evec1(:, fst1:lst1), evec2(:, fst2:lst2), &
                 zone, dHmat_mt_basis(:, :, ias), zone, dHmat, &
                 diagonal_only=diag )
        end do
      end do

      ! ** interstitial part
      if( present( Gset ) ) then
        call me_ir_mat( Gpset1, ip, Gpset2, ip, &
               evec1(:, fst1:lst1), evec2(:, fst2:lst2), &
               zone, dpot_cfun_ig, zone, dHmat, &
               Gset_op=Gset, diagonal_only=diag )
        do i = 1, 3
          call me_ir_mat( Gpset1, ip, Gpset2, ip, &
                 evec1(:, fst1:lst1), evec2(:, fst2:lst2), &
                 zone, dkin_cfun_ig, zone, dHmat, &
                 left_gradient=i, right_gradient=i, Gset_op=Gset, diagonal_only=diag )
        end do
      else
        call me_ir_mat( Gpset1, ip, Gpset2, ip, &
               evec1(:, fst1:lst1), evec2(:, fst2:lst2), &
               zone, dpot_cfun_ig, zone, dHmat, &
               diagonal_only=diag )
        do i = 1, 3
          call me_ir_mat( Gpset1, ip, Gpset2, ip, &
                 evec1(:, fst1:lst1), evec2(:, fst2:lst2), &
                 zone, dkin_cfun_ig, zone, dHmat, &
                 left_gradient=i, right_gradient=i, diagonal_only=diag )
        end do
      end if
    end subroutine dfpt_eig_gen_dHmat

    subroutine dfpt_eig_geteval( vpl, feval, pset, band_range, eval )
      use block_data_file, only: block_data_file_type
      use mod_kpointset, only: k_set
      !> wavevector \({\bf p}_0\) in lattice coordinates
      real(dp), intent(in) :: vpl(3)
      !> file with eigenvalues
      type(block_data_file_type), intent(inout) :: feval
      !> set of \({\bf p}\) vectors the file corresponds to
      type(k_set), intent(in) :: pset
      !> band range
      integer, intent(in) :: band_range(2)
      !> eigenvalues at \({\bf p}_0\)
      real(dp), allocatable, intent(out) :: eval(:)

      integer :: ip, isym, edim(1)

      real(dp), allocatable :: rdata(:)

      ! allocate eigenvectors
      edim = feval%get_block_shape()
      allocate( eval(band_range(1):band_range(2)) )
      allocate( rdata(edim(1)) )
      ! find p-point in set
      call findkptinset( vpl, pset, isym, ip )
      call terminate_if_false( ip > 0, '(dfpt_eig_geteval) &
        Requested point not found in set.' )
      ! read eigenvectors
      call feval%read( ip, rdata )
      eval = rdata(band_range(1):band_range(2))
      deallocate( rdata )
    end subroutine dfpt_eig_geteval

    subroutine dfpt_eig_getevec( vpl, vgpl, fevec, pset, Gpset, band_range, evec )
      use block_data_file, only: block_data_file_type
      use mod_kpointset, only: k_set, Gk_set
      !> wavevector \({\bf p}_0\) in lattice coordinates
      real(dp), intent(in) :: vpl(3)
      !> \({\bf G+p}_0\) vectors in lattice coordinates
      real(dp), intent(in) :: vgpl(3,*)
      !> file with eigenvectors
      type(block_data_file_type), intent(inout) :: fevec
      !> set of \({\bf p}\) vectors the file corresponds to
      type(k_set), intent(in) :: pset
      !> \({\bf G+p}\) vectors for the \({\bf p}\) vectors
      type(Gk_set), intent(in) :: Gpset
      !> band range
      integer, intent(in) :: band_range(2)
      !> eigenvectors at \({\bf p}_0\)
      complex(dp), allocatable, intent(out) :: evec(:,:)

      integer :: ip, isym, vdim(2)

      complex(dp), allocatable :: zdata(:,:)

      ! allocate eigenvectors
      vdim = fevec%get_block_shape()
      allocate( evec(vdim(1), band_range(1):band_range(2)) )
      allocate( zdata(vdim(1), vdim(2)) )
      ! find p-point in set
      call findkptinset( vpl, pset, isym, ip )
      call terminate_if_false( ip > 0, '(dfpt_eig_getevec) &
        Requested point not found in set.' )
      ! read eigenvectors
      call fevec%read( ip, zdata )
      evec = zdata(:, band_range(1):band_range(2))
      deallocate( zdata )
      if( all( abs( vpl - pset%vkl(:, ip) ) < 1e-6_dp ) ) return
      ! rotate eigenvectors
      call rotate_evecfv( isym, pset%vkl(:, ip), vpl, &
             Gpset%ngk(1, ip), Gpset%vgkl(:, :, 1, ip), vgpl, &
             evec, size( evec, dim=1 ), size( evec, dim=2 ) )
    end subroutine dfpt_eig_getevec

    !> This subroutine computes the radial muffin-tin integrals times Gaunt coefficients
    !> for the overlap and Hamiltonian matrix.
    subroutine gen_overlap_hamiltonian_mt_basis( Smat_mt_basis, Hmat_mt_basis, pot_mt, &
        lmax_apw, lmax_pot, kinetic, potential )
      use constants, only: y00, zzero, zone
      use physical_constants, only: alpha
      use mod_muffin_tin, only: nrmtmax, nrmt
      use mod_atoms, only: nspecies, natoms, idxas
      use modinput
      !> overlap radial muffin-tin integrals times Gaunt coefficients
      complex(dp), allocatable, intent(inout) :: Smat_mt_basis(:,:,:)
      !> Hamiltonian radial muffin-tin integrals times Gaunt coefficients
      complex(dp), allocatable, intent(inout) :: Hmat_mt_basis(:,:,:)
      !> muffin-tin effective potential as real spherical harmonics expansion
      real(dp), intent(in) :: pot_mt(:,:,:)
      !> maximum angular momentum \(l\) for APWs and potential expansion (default: from input file)
      integer, optional, intent(in) :: lmax_apw, lmax_pot
      !> include kinetic energy contribution (default: true)
      logical, optional, intent(in) :: kinetic
      !> include effective potential contribution (default: true)
      logical, optional, intent(in) :: potential

      integer :: lmaxapw, lmaxpot, is, ia, ias, i
      real(dp), allocatable :: rfun(:,:)
      logical :: kin, pot

      lmaxapw = dfpt_lmaxapw
      if( present( lmax_apw ) ) lmaxapw = lmax_apw
      lmaxpot = dfpt_lmaxvr
      if( present( lmax_pot ) ) lmaxpot = lmax_pot
      kin = .true.
      if( present( kinetic ) ) kin = kinetic
      pot = .true.
      if( present( potential ) ) pot = potential

      ! allocate integrals
      call me_mt_alloc( Smat_mt_basis )
      call me_mt_alloc( Hmat_mt_basis )

      allocate( rfun(1, nrmtmax) )

      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)

          ! overlap
          rfun = 1.0_dp / y00
          call me_mt_prepare( is, ias, 0, zone, rfun, zzero, Smat_mt_basis(:, :, ias) )

          ! Hamiltonian
          ! kinetic energy
          if( kin ) then
            rfun = 0.5_dp / y00
            if( input%groundstate%ValenceRelativity /= 'none' ) &
              rfun(1, 1:nrmt(is)) = rfun(1, 1:nrmt(is)) / (1.0_dp - 0.5_dp * alpha**2 * pot_mt(1, 1:nrmt(is), ias) * y00)
            do i = 1, 3
              call me_mt_prepare( is, ias, 0, zone, rfun, zone, Hmat_mt_basis(:, :, ias), &
                left_gradient=i, right_gradient=i )
            end do
          end if
          ! potential
          if( pot ) then
            call me_mt_prepare( is, ias, lmaxpot, zone, pot_mt(:, :, ias), zone, Hmat_mt_basis(:, :, ias) )
          end if

        end do
      end do
      deallocate( rfun )
    end subroutine gen_overlap_hamiltonian_mt_basis

end module dfpt_eigensystem
