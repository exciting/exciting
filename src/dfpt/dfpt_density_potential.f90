!> This module contains procedures for the calculation of the response
!> of the electronic density and the effective potential upon a generic
!> perturbation \(\delta\).
module dfpt_density_potential
  use dfpt_variables

  use precision, only: dp
  use asserts, only: assert
  use modmpi, only: terminate_if_false
  use mod_kpointset, only: G_set

  implicit none
  private

  ! LOCAL VARIABLES
  !> GGA or LDA xc-functional
  logical :: gga

  !> exchange-correlation kernel in muffin-tin spheres (in spherical coordinates)
  real(dp), allocatable :: xc_kernel_mt(:,:,:,:)
  !> exchange-correlation kernel in interstitial region
  real(dp), allocatable :: xc_kernel_ir(:,:)
  !> set of G-vectors for products of functions
  type(G_set) :: prod_Gset
  !> lmax and lmmax for products of functions
  integer :: prod_lmax, prod_lmmax
  !> forward and backward SHT matrix for products
  real(dp), allocatable :: rfsht(:,:), rbsht(:,:)

  ! for GGA only
  !> muffin-tin gradient of density on doubled real space SHT grid
  real(dp), public, allocatable :: grhomt2(:,:,:,:)
  !> interstitial gradient of density on doubled real space FFT grid
  real(dp), allocatable :: grhoir2(:,:)

  public :: dfpt_rhopot_init, dfpt_rhopot_free
  public :: dfpt_rhopot_gen_dpot
  public :: dfpt_rhopot_mixpack
  public :: dfpt_rhopot_gen_drho_mt, dfpt_rhopot_drho_k
  public :: apply_xckernel_mt, apply_xckernel_ir

  contains

    !> This subroutine initializes variables for the calculation of density
    !> and potential response that remain constant during the entire calculation.
    !>
    !> This includes:
    !> 
    !> * setting of module variables
    !> * calculation of the exchange-correlation kernel \(f_{\rm xc}({\bf r},{\bf r}')\)
    subroutine dfpt_rhopot_init
      use mod_potential_and_density, only: rhomt, rhoir
      use mod_SHT, only: gen_rshtmat

      ! G-set for products of functions
      prod_Gset = dfpt_2Gset

      ! lmax and lmmax for products of functions
      prod_lmax = 1 * dfpt_lmaxvr
      prod_lmmax = ( prod_lmax + 1 )**2

      ! SHT matrices
      call gen_rshtmat( prod_lmax, rfsht, rbsht )

      ! generate xc kernel
      call gen_xc_kernel( rhomt, rhoir, xc_kernel_mt, xc_kernel_ir )
    end subroutine dfpt_rhopot_init

    !> This subroutine frees memory from the module variables.
    subroutine dfpt_rhopot_free
      use mod_kpointset, only: delete_G_vectors
      if( allocated( xc_kernel_mt ) ) deallocate( xc_kernel_mt )
      if( allocated( xc_kernel_ir ) ) deallocate( xc_kernel_ir )
      if( allocated( grhomt2 ) ) deallocate( grhomt2 )
      if( allocated( grhoir2 ) ) deallocate( grhoir2 )
      if( allocated( rfsht ) ) deallocate( rfsht )
      if( allocated( rbsht ) ) deallocate( rbsht )
      call delete_G_vectors( prod_Gset )
    end subroutine dfpt_rhopot_free

    !> This subroutine calculates the effective potential response from a given
    !> density response using Weinert's method.
    !>
    !> For a generic perturbation operator \(\delta\), the variation/response
    !> of the effective potential is the sum of the Hartree potential response
    !> \[\delta V_{\rm H}({\bf r}) = \delta \int \frac{n({\bf r}')}{| {\bf r} - {\bf r}'|} {\rm d}{\bf r}'
    !>   = \int \frac{\delta n({\bf r}')}{|{\bf r} - {\bf r}'|} {\rm d}{\bf r}' \;,\]
    !> the external potential response
    !> \[\delta V_{\rm ext}({\bf r}) = -\delta \sum_{\alpha,{\bf R}} \frac{Z_\alpha}{|{\bf r} 
    !>   - {\bf \tau}_\alpha - {\bf R}|} \;,\]
    !> which together form the Coulomb potential response, and the exchange-correlation
    !> potential response
    !> \[\delta V_{\rm xc}({\bf r}) = \int f_{\rm xc}({\bf r},{\bf r}') \delta n({\bf r}') {\rm d}{\bf r}' \;,\]
    !> where \(f_{\rm xc}({\bf r},{\bf r}')\) is the exchange-correlation kernel.
    subroutine dfpt_rhopot_gen_dpot( drho_mt, drho_ir, dpot_mt, dpot_ir, Gpset, jlgpr, ylmgp, sfacgp, &
        coulomb, xc, mt_mp_add )
      use weinert
      use constants, only: zzero
      use mod_atoms, only: natmtot, nspecies, natoms, idxas, spr
      use mod_muffin_tin, only: idxlm, nrmtmax, nrmt, rmt
      use modinput
      !> muffin-tin density response as complex spherical harmonics expansion
      complex(dp), intent(in) :: drho_mt(:,:,:)
      !> interstitial density response on real space FFT grid
      complex(dp), intent(in) :: drho_ir(:)
      !> muffin-tin effective potential response as complex spherical harmonics expansion
      complex(dp), intent(out) :: dpot_mt(:,:,:)
      !> interstitial effective potential response on real space FFT grid
      complex(dp), intent(out) :: dpot_ir(:)
      !> set of \({\bf G+p}\) vectors the interstitial density / potential response is expanded in
      type(G_set), intent(in) :: Gpset
      !> spherical Bessel functions \(j_l(|{\bf G+p}| R_\alpha)\)
      real(dp), intent(in) :: jlgpr(0:,:,:)
      !> spherical harmonics \(Y_{lm}(\widehat{\bf G+p})\)
      complex(dp), intent(in) :: ylmgp(:,:)
      !> structure factors \({\rm e}^{{\rm i} ({\bf G+p}) \cdot {\bf \tau}_\alpha}\)
      complex(dp), intent(in) :: sfacgp(:,:)
      !> if true, the Coulomb / exchange-correlation potential response is calculated (default: `.true.`)
      logical, optional, intent(in) :: coulomb, xc
      !> additional muffin-tin multipole moments that do not come from the density response (default: none)
      complex(dp), optional, intent(in) :: mt_mp_add(dfpt_lmmaxvr, natmtot)

      integer :: is, ia, ias, l, m, lm, ir
      logical :: vcoul, vxc

      real(dp), allocatable :: rl(:,:), rfmt(:,:,:)
      complex(dp), allocatable :: zfmt(:,:), zfir(:)
      complex(dp), allocatable :: rho_mt_mp(:,:), rho_i_mp(:,:)
      complex(dp), allocatable :: drho_surf(:,:)

      vcoul = .true.
      if( present( coulomb ) ) vcoul = coulomb
      vxc = .true.
      if( present( xc ) ) vxc = xc

      dpot_mt = zzero
      dpot_ir = zzero

      ! **** response of Coulomb potential
      if( vcoul ) then
        ! ** initialize auxilliary variables
        allocate( rho_mt_mp(dfpt_lmmaxvr, natmtot), source=zzero )
        allocate( rho_i_mp(dfpt_lmmaxvr, natmtot), source=zzero )
        allocate( drho_surf(dfpt_lmmaxvr, natmtot), source=zzero )

        ! ** Fourier transform density response to reciprocal space
        allocate( zfir, source=drho_ir )
        call zfftifc( 3, dfpt_Gset%ngrid, -1, zfir )

        ! ** find interstitial multipole moments
        call multipoles_ir( dfpt_lmaxvr, Gpset%ngvec, Gpset%gc, Gpset%ivg, &
               jlgpr, ylmgp, sfacgp, &
               dfpt_Gset%intgv, dfpt_Gset%ivgig, dfpt_Gset%igfft, zfir, &
               rho_i_mp )

        ! ** find the muffin tin multipole moments 
        ! ** and solution to Poisson's equation in MT spheres
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            call poisson_and_multipoles_mt( dfpt_lmaxvr, nrmt(is), spr(1:nrmt(is), is), drho_mt(:, :, ias), &
                   dpot_mt(:, :, ias), rho_mt_mp(:, ias) )
          end do
        end do

        ! ** build total multipole moments
        ! difference between MT and IR multipoles
        rho_mt_mp = rho_mt_mp - rho_i_mp
        ! add additional multipoles if given
        if( present( mt_mp_add ) ) rho_mt_mp = rho_mt_mp + mt_mp_add

        ! ** find interstitial response of Coulomb potential
        ! find Fourier coefficient from multipole moments
        call poisson_ir( dfpt_lmaxvr, input%groundstate%npsden, Gpset%ngvec, Gpset%gc, Gpset%ivg, &
               jlgpr, ylmgp, sfacgp, dfpt_Gset%intgv, dfpt_Gset%ivgig, dfpt_Gset%igfft, &
               zfir, rho_mt_mp, dpot_ir )
        ! find interstitial potential response on MT sphere boundaries
        call surface_ir( dfpt_lmaxvr, Gpset%ngvec, Gpset%gc, Gpset%ivg, &
               jlgpr, ylmgp, sfacgp, &
               dfpt_Gset%intgv, dfpt_Gset%ivgig, dfpt_Gset%igfft, dpot_ir, &
               drho_surf )
        ! transform potential response back to real space
        call zfftifc( 3, dfpt_Gset%ngrid, 1, dpot_ir )

        ! ** add homogeneous solution of Poisson's equation in MT spheres
        ! ** in order to make potential response smooth at sphere boundaries
        allocate( rl(nrmtmax, 0:dfpt_lmaxvr) )
        do is = 1, nspecies
          rl(:, 0) = 1.0_dp
          rl(1:nrmt(is), 1) = spr(1:nrmt(is), is) / rmt(is)
          do l = 2, dfpt_lmaxvr
            rl(:, l) = rl(:, l-1) * rl(:, 1)
          end do

          do ia = 1, natoms(is)
            ias = idxas(ia, is)
!$omp parallel default( shared ) private( l, m, lm, ir )
!$omp do
            do l = 0, dfpt_lmaxvr
              do m = -l, l
                lm = idxlm(l, m)
                do ir = 1, nrmt(is)
                  dpot_mt( lm, ir, ias ) = dpot_mt(lm, ir, ias) + rl(ir, l) * drho_surf(lm, ias)
                end do
              end do
            end do
!$omp end do
!$omp end parallel

          end do
        end do

        ! delete auxilliary variables
        deallocate( rho_mt_mp, rho_i_mp, drho_surf, zfir, rl )
      end if

      ! **** response of exchange-correlation potential
      if( vxc ) then
        ! ** muffin-tin contribution
        ! initialize auxilliary variables
        allocate( rfmt(dfpt_lmmaxvr, nrmtmax, 4) )
        allocate( zfmt(dfpt_lmmaxvr, nrmtmax) )
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            ! decompose density response into real and imaginary part
            do ir = 1, nrmt(is)
              call decompose_zflm( dfpt_lmaxvr, drho_mt(:, ir, ias), rfmt(:, ir, 1), rfmt(:, ir, 2) )
            end do
            ! apply xc-kernel
            call apply_xckernel_mt( is, ia, rfmt(:, :, 1), rfmt(:, :, 3) )
            call apply_xckernel_mt( is, ia, rfmt(:, :, 2), rfmt(:, :, 4) )
            do ir = 1, nrmt(is)
              ! join real and imaginary part to spherical-harmonic expansion
              call compose_zflm( dfpt_lmaxvr, rfmt(:, ir, 3), rfmt(:, ir, 4), zfmt(:, ir) )
              ! add result to potential response
              dpot_mt(:, ir, ias) = dpot_mt(:, ir, ias) + zfmt(:, ir)
            end do
          end do
        end do
        deallocate( rfmt, zfmt )
  
        ! ** interstitial contribution
        ! apply xc-kernel
        call apply_xckernel_ir( drho_ir, Gpset, dpot_ir )
      end if
    end subroutine dfpt_rhopot_gen_dpot

    !> This subroutine calculates the contribution to the electron density response
    !> coming from a single \({\bf k}\) point.
    !> 
    !> For a generic perturbation operator \(\delta\), this is given by
    !> \[\delta n_{\bf k}({\bf r}) = 2 \sum_n f_{n{\bf k}} \, \psi_{n{\bf k}}^\ast({\bf r}) \delta\psi_{n{\bf k}}({\bf r})
    !>   + \sum_n \delta f_{n{\bf k}} \, |\psi_{n{\bf k}}({\bf r})|^2 \;,\]
    !> where \(f_{n{\bf k}}\) are the occupation numbers.
    !>
    !> In the muffin-tins, we assume that the wavefunctions (and their response) 
    !> are given as a spherical harmonics expansion
    !> \[\psi_{n{\bf k}}({\bf r}_\alpha) = \sum_\lambda C^{n{\bf k},\alpha}_\lambda \, \varphi_\lambda(r_\alpha) \, 
    !>   Y_{l_\lambda m_\lambda}(\hat{\bf r}_\alpha)\]
    !> and the muffin-tin density response is described in terms of the density matrix response
    !> \[\delta D^{\bf k,\alpha}_{\lambda \lambda'} = 2 \sum_n f_{n{\bf k}} \, {C^{n{\bf k},\alpha}_\lambda}^\ast \, 
    !>   \delta C^{n{\bf k},\alpha}_{\lambda'}
    !>   + \sum_n \delta f_{n{\bf k}} \, {C^{n{\bf k},\alpha}_\lambda}^\ast \, C^{n{\bf k},\alpha}_{\lambda'} \;.\]
    !>
    !> @note The result is added to the input variables.
    !> 
    !> @warning This subroutine does not account for contributions coming from the variation of the
    !>          basis functions upon the perturbation and only accounts for the variation of the 
    !>          expansion coefficients of the wavefunctions.
    subroutine dfpt_rhopot_drho_k( ik, kset, Gkset1, Gkset2, fst, lst, occk, eveck, deveck, apwalmk1, apwalmk2, drho_mat, drho_ir, &
        docck )
      use constants, only: zzero, zone
      use mod_kpointset, only: k_set, Gk_set
      use mod_atoms, only: nspecies, natoms, idxas
      use mod_APW_LO, only: nlotot
      use mod_eigensystem, only: idxlo
      use mod_lattice, only: omega
      use modinput
      !> index of the wavevector \({\bf k}\) in the set
      integer, intent(in) :: ik
      !> set of \({\bf k}\) vectors
      type(k_set), intent(in) :: kset
      !> set of \({\bf G+k}\) vectors for the wavefunctions
      type(Gk_set), intent(in) :: Gkset1
      !> set of \({\bf G+k}\) vectors for the wavefunction response
      type(Gk_set), intent(in) :: Gkset2
      !> first and last state to consider
      integer, intent(in) :: fst, lst
      !> occupation numbers at \({\bf k}\)
      real(dp), intent(in) :: occk(:)
      !> eigenvectors at \({\bf k}\)
      complex(dp), intent(in) :: eveck(:,:)
      !> eigenvector response at \({\bf k}\)
      complex(dp), intent(in) :: deveck(:,:)
      !> (L)APW matching coefficients \(A^\alpha_{{\bf G+p},lm,\xi}\) for wavefunctions
      complex(dp), intent(in) :: apwalmk1(:,:,:,:)
      !> (L)APW matching coefficients \(A^\alpha_{{\bf G+p},lm,\xi}\) for wavefunction response
      complex(dp), intent(in) :: apwalmk2(:,:,:,:)
      !> muffin-tin density response matrix
      complex(dp), intent(inout) :: drho_mat(:,:,:)
      !> interstitial density response on real space FFT grid
      complex(dp), intent(inout) :: drho_ir(:)
      !> occupation number response at \({\bf k}\)
      real(dp), optional, intent( in) :: docck(:)

      integer :: ngk1, ngk2, i, ist, nst
      integer :: is, ia, ias
      real(dp) :: t1, t2

      integer, allocatable :: bands(:)
      real(dp), allocatable :: wgt1(:), wgt2(:)
      complex(dp), allocatable :: evecmt1(:,:), evecmt2(:,:), devecmt(:,:)
      complex(dp), allocatable :: zfft(:,:)

      ngk1 = Gkset1%ngk(1, ik)
      ngk2 = Gkset2%ngk(1, ik)

      ! find k-point weights for bands and bands that contribute
      allocate( wgt1(fst:lst), wgt2(fst:lst) )
      wgt1 = [(2*kset%wkpt(ik)*occk(i), i=fst, lst)]
      if( present( docck ) ) then
        wgt2 = [(kset%wkpt(ik)*docck(i), i=fst, lst)]
        bands = pack( [(i, i=fst, lst)], &
                      [(wgt1(i) > input%groundstate%epsocc .or. &
                        wgt2(i) > input%groundstate%epsocc, i=fst, lst)] )
      else
        bands = pack( [(i, i=fst, lst)], &
                      [(wgt1(i) > input%groundstate%epsocc, i=fst, lst)] )
      end if
      nst = size( bands )
      if( nst == 0 ) return

      ! **** muffin-tin density response matrix
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! generate shortened MT eigenvector response at current k-point
          call mt_basis%transform_evec( is, ngk2, nlotot, idxlo(:, :, ias), apwalmk2(:, :, :, ias), deveck(:, bands), devecmt )
          ! generate shortened MT eigenvector at current k-point
          call mt_basis%transform_evec( is, ngk1, nlotot, idxlo(:, :, ias), apwalmk1(:, :, :, ias), eveck(:, bands), evecmt1 )
          ! contribution from occupation response
          if( present( docck ) ) then
            allocate( evecmt2, source=evecmt1 )
            do i = 1, nst
              ist = bands(i)
              evecmt2(:, i) = wgt2(ist) * evecmt1(:, i)
            end do
            call zgemm( 'n', 'c', mt_basis%n_basis_fun(is), mt_basis%n_basis_fun(is), nst, zone, &
                   evecmt2, size(evecmt2, dim=1), &
                   evecmt1, size(evecmt1, dim=1), zone, &
                   drho_mat(:, :, ias), size(drho_mat, dim=1) )
            deallocate( evecmt2 )
          end if
          ! sum up MT eigenvectors over occupied states
          do i = 1, nst
            ist = bands(i)
            evecmt1(:, i) = wgt1(ist) * evecmt1(:, i)
          end do
          call zgemm( 'n', 'c', mt_basis%n_basis_fun(is), mt_basis%n_basis_fun(is), nst, zone, &
                 devecmt, size(devecmt, dim=1), &
                 evecmt1, size(evecmt1, dim=1), zone, &
                 drho_mat(:, :, ias), size(drho_mat, dim=1) )
        end do
      end do
      if( allocated( evecmt1 ) ) deallocate( evecmt1 )
      if( allocated( devecmt ) ) deallocate( devecmt )

      ! **** interstitial density response
      ! sum over occupied states
      allocate( zfft(dfpt_Gset%ngrtot, 2) )
      do i = 1, nst
        ist = bands(i)
        t1 = wgt1(ist) / omega
        t2 = wgt2(ist) / omega
        zfft = zzero
        ! Fourier transform wavefunction to real space
        zfft(dfpt_Gset%igfft(Gkset1%igkig(1:ngk1, 1, ik)), 1) = eveck(1:ngk1, ist)
        call zfftifc( 3, dfpt_Gset%ngrid, 1, zfft(:, 1) )
        ! Fourier transform wavefunction response to real space
        zfft(dfpt_Gset%igfft(Gkset2%igkig(1:ngk2, 1, ik)), 2) = deveck(1:ngk2, ist)
        call zfftifc( 3, dfpt_Gset%ngrid, 1, zfft(:, 2) )
        ! multiply wavefunction and wavefunction response and add to density response
        if( present( docck ) ) then
          drho_ir = drho_ir + conjg( zfft(:, 1) ) * (t1 * zfft(:, 2) + t2 * zfft(:, 1))
        else
          drho_ir = drho_ir + t1 * conjg( zfft(:, 1) ) * zfft(:, 2)
        end if
      end do
      deallocate( zfft )
      deallocate( wgt1, wgt2 )
    end subroutine dfpt_rhopot_drho_k

    !> This subroutine calculates the muffin-tin density response 
    !> from the given density response matrix \(\delta{\bf D}^\alpha\).
    !>
    !> This is given by
    !> \[\delta n({\bf r}_\alpha) = \sum_{l,m} \delta n_{lm}(r_\alpha) \, Y_{lm}(\hat{\bf r}_\alpha) \;,\]
    !> with
    !> \[\delta n_{lm}(r_\alpha) = \sum_{\lambda,\lambda'} \delta D^\alpha_{\lambda \lambda'} \,
    !>   \varphi_{\lambda}(r_\alpha) \, \varphi_{\lambda'}^\ast(r_\alpha) \, 
    !>   \langle Y_{l_\lambda m_\lambda} | Y_{lm} | Y_{l_{\lambda'} m_{\lambda'}} \rangle \;.\]
    subroutine dfpt_rhopot_gen_drho_mt( drho_mat, drho_mt )
      use constants, only: zzero
      use mod_muffin_tin, only: nrmtmax
      use mod_atoms, only: nspecies, natoms, idxas
      use mod_muffin_tin, only: nrmt
      use gaunt
      !> muffin-tin density response matrix
      complex(dp), intent(in) :: drho_mat(:,:,:)
      !> (soft) muffin-tin density response as complex spherical harmonics expansion
      complex(dp), intent(out) :: drho_mt(:,:,:)

      integer :: is, ia, ias
      integer :: l1, l2, m1, m2, lm1, lm2, lm3, lam1, lam2, idx1, idx2, i
      complex(dp) :: z1
      type(non_zero_gaunt_real), pointer :: yyy

      real(dp), allocatable :: f1(:), f2(:)
      complex(dp), allocatable :: drhotmp(:,:)

      drho_mt = zzero

      yyy => gaunt_coeff_yyy
      allocate( drhotmp(nrmtmax, dfpt_lmmaxvr) )
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          drhotmp = zzero

!$omp parallel default( shared ) private( l1, l2, m1, m2, lm1, lm2, lm3, lam1, lam2, idx1, idx2, f1, f2, z1, i ) reduction( +:drhotmp )
!$omp do collapse(2)
          do l1 = 0, dfpt_lmaxapw
            do l2 = 0, dfpt_lmaxapw

              do lam1 = 1, mt_basis%n_rad_fun(l1, is)
                f1 = mt_basis%get_rad_fun( l1, is, ias, lam1 )
                do lam2 = 1, mt_basis%n_rad_fun(l2, is)
                  f2 = mt_basis%get_rad_fun( l2, is, ias, lam2 )
                  f2 = f1 * f2

                  do m1 = -l1, l1
                    lm1 = l1 * (l1 + 1) + m1 + 1
                    idx1 = mt_basis%idx_basis_fun(lm1, lam1, is)
                    do m2 = -l2, l2
                      lm2 = l2 * (l2 + 1) + m2 + 1
                      idx2 = mt_basis%idx_basis_fun(lm2, lam2, is)

                      do i = 1, yyy%num(lm2, lm1)
                        lm3 = yyy%lm2(i, lm2, lm1)
                        if( lm3 > dfpt_lmmaxvr ) exit
                        z1 = yyy%val(i, lm2, lm1) * drho_mat(idx2, idx1, ias)
                        drhotmp(1:nrmt(is), lm3) = drhotmp(1:nrmt(is), lm3) + z1 * f2(1:nrmt(is))
                      end do

                    end do
                  end do

                end do
              end do
            
            end do
          end do
!$omp end do
!$omp end parallel
          
          drho_mt(:, :, ias) = transpose( drhotmp )

        end do
      end do
      if( allocated( f1 ) ) deallocate( f1 )
      if( allocated( f2 ) ) deallocate( f2 )
      deallocate( drhotmp )
    end subroutine dfpt_rhopot_gen_drho_mt

    !> This subroutine (un)packs the density or potential response for the mixing.
    subroutine dfpt_rhopot_mixpack( drho_mt, drho_ir, dpot_mt, dpot_ir, tpack, nf, n, nu )
      use modinput
      !> muffin-tin density / potential response as complex spherical harmonics expansion
      complex(dp), intent(inout) :: drho_mt(:,:,:,:), dpot_mt(:,:,:,:)
      !> interstitial density / potential response on real space FFT grid
      complex(dp), intent(inout) :: drho_ir(:,:), dpot_ir(:,:)
      !> `.true.` for packing and `.false.` for unpacking
      logical, intent(in) :: tpack
      !> number of functions
      integer, intent(in) :: nf
      !> number of elements per function
      integer, intent(out) :: n
      !> (un)packed function
      real(dp), intent(inout) :: nu(*)

      n = 0
      ! density for mixing
      if( input%groundstate%mixerswitch == 2 ) then
        call zfpack( tpack, n, dfpt_lmmaxvr, 1, nf, drho_mt, drho_ir, nu )
      ! potential for mixing
      else
        call zfpack( tpack, n, dfpt_lmmaxvr, 1, nf, dpot_mt, dpot_ir, nu )
      end if
    end subroutine dfpt_rhopot_mixpack

    !> This subroutine calculates the exchange-correlation kernel \(f_{\rm xc}({\bf r},{\bf r}')\)
    !> from a given electronic density using the `libxc` library.
    !>
    !> The xc-kernel is defined as the second functional derivative of the xc-energy
    !> \[ f_{\rm xc}({\bf r},{\bf r}') = \frac{ \delta^2 E_{\rm xc}[n({\bf r})]}{\delta n({\bf r}) \, \delta n({\bf r}')} \;.\]
    !> For both LDA and GGA xc-functionals, the xc-kernel is local, i.e., 
    !> \(f_{\rm xc}({\bf r},{\bf r}') = f_{\rm xc}({\bf r}) \, \delta({\bf r} - {\bf r}')\).
    !> The quantities `libxc` provides are the derivatives of the xc-energy density
    !> per unit volume, \(\epsilon_{\rm xc}[n({\bf r})]\), w.r.t. to electron density \(n({\bf r})\)
    !> and its gradient \(\nabla n({\bf r})\) in the form of
    !> \[ \frac{\partial \epsilon_{\rm xc}[n({\bf r})]}{\partial n({\bf r})} \;,
    !>    \frac{\partial \epsilon_{\rm xc}[n({\bf r})]}{\partial \sigma({\bf r})} \;,
    !>    \frac{\partial^2 \epsilon_{\rm xc}[n({\bf r})]}{\partial n({\bf r})^2} \;,
    !>    \frac{\partial^2 \epsilon_{\rm xc}[n({\bf r})]}{\partial n({\bf r})\, \partial \sigma({\bf r})} \;,
    !>    \frac{\partial^2 \epsilon_{\rm xc}[n({\bf r})]}{\partial \sigma({\bf r})^2} \;, \]
    !> where \(\sigma({\bf r}) = (\nabla n({\bf r}))^2\).
    !>
    !> In the case of LDA, the xc-kernel is returned in a single part 
    !> \[ f^{\rm LDA}_{\rm xc,1}({\bf r}) = \frac{\partial^2 \epsilon_{\rm xc}[n({\bf r})]}{\partial n({\bf r})^2} \,.\]
    !> In the case of GGA, the xc-kernel is returned in four parts
    !> \[\begin{align*}
    !>   f^{\rm GGA}_{\rm xc,1}({\bf r}) &= \frac{\partial^2 \epsilon_{\rm xc}[n({\bf r})]}{\partial n({\bf r})^2} 
    !>     - 2 \nabla \cdot \left( \frac{\partial^2 \epsilon_{\rm xc}[n({\bf r})]}{\partial n({\bf r})\, \partial \sigma({\bf r})} 
    !>       \nabla n({\bf r}) \right) \\
    !>   f^{\rm GGA}_{\rm xc,2}({\bf r}) &= \frac{\partial^2 \epsilon_{\rm xc}[n({\bf r})]}{\partial n({\bf r})\, \partial \sigma({\bf r})} \\
    !>   f^{\rm GGA}_{\rm xc,3}({\bf r}) &= \frac{\partial^2 \epsilon_{\rm xc}[n({\bf r})]}{\partial \sigma({\bf r})^2} \\
    !>   f^{\rm GGA}_{\rm xc,4}({\bf r}) &= \frac{\partial \epsilon_{\rm xc}[n({\bf r})]}{\partial \sigma({\bf r})} \;.
    !> \end{align*}\]
    !> See subroutines [[apply_xckernel_mt(subroutine)]] and [[apply_xckernel_ir(subroutine)]] for more details
    !> and the reason of the exact form \(f^{\rm GGA}_{\rm xc,1}({\bf r})\).
    !> @note Exchange-correlation types other than LDA and GGA are not yet supported by this subroutine.
    subroutine gen_xc_kernel( rhomt, rhoir, xc_kernel_mt, xc_kernel_ir )
      use constants, only: zzero
      use mod_potential_and_density, only: xctype, xcdescr
      use mod_muffin_tin, only: nrmtmax, nrmt
      use mod_atoms, only: natmtot, nspecies, natoms, idxas, spr
      use mod_symmetry, only: nsymcrys, symmetrize_real_mt, symmetrize_real_ir
#ifdef LIBXC
      use xc_f90_lib_m
#endif
      !> muffin-tin density as real spherical harmonics expansion
      real(dp), intent(in) :: rhomt(:,:,:)
      !> interstitial density on real space FFT grid
      real(dp), intent(in) :: rhoir(:)
      !> (parts of) muffin-tin xc kernel as real spherical harmonics expansion
      real(dp), allocatable, intent(out) :: xc_kernel_mt(:,:,:,:)
      !> (parts of) interstitial xc kernel on real space FFT grid
      real(dp), allocatable, intent(out) :: xc_kernel_ir(:,:)

      integer :: i, j, k, np, nr
      integer :: is, ia, ias, ig, ifg
      integer :: xctype_libxc(2), xcf

      real(dp), allocatable :: grhomt(:,:,:), rhomt2(:,:), rhoir2(:)
      real(dp), allocatable :: grho2(:), fxc(:,:), gfxc(:,:,:)
      complex(dp), allocatable :: zfft(:,:)
      
#ifdef LIBXC
      type(xc_f90_pointer_t) :: xct
      type(xc_f90_pointer_t) :: info
#else
      call terminate_if_false( .false., '(gen_xc_kernel) &
        LIBXC library not available.' )
#endif

      ! NOTE:
      ! We want to use the LIBXC functions
      !   xc_lda_fxc and xc_gga_fxc
      ! to get the xc-kernel. Therefore, we first translate the exciting internal
      ! xc-code to the corresponding one of LIBX.
#ifdef LIBXC
      xctype_libxc = 0
      select case( abs( xctype(1) ) )
        case( 100 )  ! already LIBXC type
          xctype_libxc = xctype(2:3)
        case( 1 )    ! no density derived xc energy
          xctype_libxc = 0
        case( 2 )    ! LDA Perdew-Zunger
          xctype_libxc(1) = XC_LDA_X
          xctype_libxc(2) = XC_LDA_C_PZ
        case( 3 )    ! LDA Perdew-Wang
          xctype_libxc(1) = XC_LDA_X
          xctype_libxc(2) = XC_LDA_C_PW
        case( 4 )    ! LDA Slater X-alpha
          xctype_libxc(1) = XC_LDA_X
          xctype_libxc(2) = XC_LDA_C_XALPHA
        case( 5 )    ! LDA Barth-Hedin
          xctype_libxc(1) = XC_LDA_X
          xctype_libxc(2) = XC_LDA_C_VBH
        case( 20 )   ! GGA PBE
          xctype_libxc(1) = XC_GGA_X_PBE
          xctype_libxc(2) = XC_GGA_C_PBE
        case( 21 )   ! GGA revised PBE, Zhang-Yang
          xctype_libxc(1) = XC_GGA_X_PBE_R
          xctype_libxc(1) = XC_GGA_C_PBE
        case( 22 )   ! GGA PBEsol
          xctype_libxc(1) = XC_GGA_X_PBE_SOL
          xctype_libxc(2) = XC_GGA_C_PBE_SOL
        case( 26 )   ! GGA Wo-Cohen
          xctype_libxc(1) = XC_GGA_X_WC
          xctype_libxc(2) = XC_GGA_C_PBE
        case( 30 )   ! GGA Armiento-Mattsson
          xctype_libxc(1) = XC_GGA_X_AM05
          xctype_libxc(2) = XC_GGA_C_AM05
        ! this would be available in a newer version of LIBXC
        !case( 300 )  ! GGA acPBE
        !  xctype_libxc(1) = XC_GGA_X_BCGP
        !  xctype_libxc(2) = XC_GGA_C_ACGGA
        case default
          call terminate_if_false( .false., '(gen_xc_kernel) &
            xc kernel not available for this xc-type. '// &
            trim( adjustl( xcdescr ) ) )
      end select

      ! check, if we need density gradients
      k = 1
      do i = 1, 2
        if( xctype_libxc(i) == 0 ) cycle
        ! get the xc functional family (LDA or GGA)
        xcf = xc_f90_family_from_id( xctype_libxc(i) )
        select case( xcf )
          case( XC_FAMILY_LDA )
            k = 1 ! only one part for LDA
          case( XC_FAMILY_GGA, XC_FAMILY_HYB_GGA )
            gga = .true.
            k = 4 ! four parts for GGA
          case default
            call terminate_if_false( .false., '(gen_xc_kernel) &
              Unsupported LIBX functional family.' )
        end select
      end do

      ! allocate kernel (and gradient) arrays
      allocate( grhomt2(prod_lmmax, nrmtmax, 3, natmtot), source=0._dp )
      allocate( grhoir2(prod_Gset%ngrtot, 3) )
      allocate( xc_kernel_mt(prod_lmmax, nrmtmax, natmtot, k) )
      allocate( xc_kernel_ir(prod_Gset%ngrtot, k) )

      ! **** GET EXCHANGE-CORRELATION KERNEL
      ! ** kernel in muffin-tin spheres
      xc_kernel_mt = 0.0_dp
      np = prod_lmmax * nrmtmax
      allocate( grhomt(dfpt_lmmaxvr, nrmtmax, 3), source=0.0_dp )
      allocate( rhomt2(prod_lmmax, nrmtmax) )
      allocate( grho2(np) )
      allocate( fxc(np, k) )

      do is = 1, nspecies
        nr = nrmt(is)
        np = prod_lmmax * nr
        do ia = 1, natoms(is)
          ias = idxas(ia, is)

          ! get density on doubled grid
          call dgemm( 'n', 'n', prod_lmmax, nr, dfpt_lmmaxvr, 1.0_dp, &
                 rbsht, prod_lmmax, &
                 rhomt(:, :, ias), dfpt_lmmaxvr, 0.0_dp, &
                 rhomt2, prod_lmmax )

          ! get density gradient and gradient^2 on doubled grid
          call gradrfmt( dfpt_lmaxvr, nr, spr(1:nr, is), dfpt_lmmaxvr, nrmtmax, rhomt(:, :, ias), grhomt )
          grho2 = 0.0_dp
          do i = 1, 3
            grhomt2(1:dfpt_lmmaxvr, :, i, ias) = grhomt(:, :, i)
            call dgemm( 'n', 'n', prod_lmmax, nr, dfpt_lmmaxvr, 1.0_dp, &
                   rbsht, prod_lmmax, &
                   grhomt(:, :, i), dfpt_lmmaxvr, 0.0_dp, &
                   fxc(:, 1), prod_lmmax )
            grho2 = grho2 + fxc(:, 1)**2
          end do
          where( grho2 < 1e-12_dp ) grho2 = 1e-12_dp

          ! get xc-kernel
          do i = 1, 2
            if( xctype_libxc(i) == 0 ) cycle
            ! initialize the xc-functional
            call xc_f90_func_init( xct, info, xctype_libxc(i), XC_UNPOLARIZED )
            ! get the xc functional family (LDA or GGA)
            xcf = xc_f90_family_from_id(xctype_libxc(i))
            select case( xcf )
              case( XC_FAMILY_LDA )
                call xc_f90_lda_fxc( xct, np, rhomt2(1, 1), fxc(1, 1) )
                where( isnan( fxc(:, 1) ) ) fxc(:, 1) = 1.0_dp
                xc_kernel_mt(:, :, ias, 1) = xc_kernel_mt(:, :, ias, 1) + reshape( fxc(:, 1), [prod_lmmax, nrmtmax] )
              case( XC_FAMILY_GGA, XC_FAMILY_HYB_GGA )
                call xc_f90_gga_vxc( xct, np, rhomt2(1, 1), grho2(1), fxc(1, 1), fxc(1, 2) )
                where( isnan( fxc(:, 2) ) ) fxc(:, 2) = 0.0_dp
                xc_kernel_mt(:, :, ias, 4) = xc_kernel_mt(:, :, ias, 4) + reshape( fxc(:, 2), [prod_lmmax, nrmtmax] )
                call xc_f90_gga_fxc( xct, np, rhomt2(1, 1), grho2(1), fxc(1, 1), fxc(1, 2), fxc(1, 3) )
                where( isnan( fxc(:, 1) ) ) fxc(:, 1) = 1.0_dp
                xc_kernel_mt(:, :, ias, 1) = xc_kernel_mt(:, :, ias, 1) + reshape( fxc(:, 1), [prod_lmmax, nrmtmax] )
                where( isnan( fxc(:, 2) ) ) fxc(:, 2) = 0.0_dp
                xc_kernel_mt(:, :, ias, 2) = xc_kernel_mt(:, :, ias, 2) + reshape( fxc(:, 2), [prod_lmmax, nrmtmax] )
                where( isnan( fxc(:, 3) ) ) fxc(:, 3) = 0.0_dp
                xc_kernel_mt(:, :, ias, 3) = xc_kernel_mt(:, :, ias, 3) + reshape( fxc(:, 3), [prod_lmmax, nrmtmax] )
              case default
            end select
          end do

          ! transform xc-kernel to spherical harmonics
          do i = 1, k
            fxc(:, i) = reshape( xc_kernel_mt(:, :, ias, i), [prod_lmmax * nrmtmax] )
            call dgemm( 'n', 'n', prod_lmmax, nr, prod_lmmax, 1.0_dp, &
                   rfsht, prod_lmmax, &
                   fxc(:, i), prod_lmmax, 0.0_dp, &
                   xc_kernel_mt(:, :, ias, i), prod_lmmax )
          end do

        end do
      end do

      deallocate( grhomt, rhomt2, grho2, fxc )

      ! ** kernel in interstitial region
      xc_kernel_ir = 0.0_dp
      np = prod_Gset%ngrtot
      allocate( rhoir2(prod_Gset%ngrtot) )
      allocate( grho2(np) )
      allocate( fxc(np, k) )
      allocate( zfft(prod_Gset%ngrtot, 0:3) )

      ! get density on doubled grid
      zfft = zzero
      zfft(1:dfpt_Gset%ngrtot, 1) = cmplx( rhoir, 0.0_dp, dp )
      call zfftifc( 3, dfpt_Gset%ngrid, -1, zfft(:, 1) )
      call dfpt_Gset%igfft2ig( zfft(:, 1), zfft(:, 2) )
      call dfpt_Gset%change_set( prod_Gset, zfft(:, 2), zfft(:, 0), 'pull', gmax=dfpt_Gset%gmaxvr )
      call prod_Gset%ig2igfft( zfft(:, 0), zfft(:, 3), gmax=dfpt_Gset%gmaxvr )
      call zfftifc( 3, prod_Gset%ngrid, 1, zfft(:, 3) )
      rhoir2 = dble( zfft(:, 3) )

      ! get density gradient and gradient^2 on doubled grid
      grho2 = 0._dp
      do i = 1, 3
        zfft(:, 1) = zzero
        do ig = 1, prod_Gset%ngvec
          if( prod_Gset%gc(ig) > dfpt_Gset%gmaxvr ) exit
          zfft(prod_Gset%igfft(ig), 1) = prod_Gset%vgc(i, ig) * cmplx( -aimag( zfft(ig, 0) ), dble( zfft(ig, 0) ), dp )
        end do
        call zfftifc( 3, prod_Gset%ngrid, 1, zfft(:, 1) )
        grhoir2(:, i) = dble( zfft(:, 1) )
        grho2 = grho2 + grhoir2(:, i)**2
      end do
      where( grho2 < 1e-12_dp ) grho2 = 1e-12_dp

      ! get xc-kernel
      do i = 1, 2
        if( xctype_libxc(i) == 0 ) cycle
        ! initialize the xc-functional
        call xc_f90_func_init( xct, info, xctype_libxc(i), XC_UNPOLARIZED )
        ! get the xc functional family (LDA or GGA)
        xcf = xc_f90_family_from_id( xctype_libxc(i) )
        select case( xcf )
          case( XC_FAMILY_LDA )
            call xc_f90_lda_fxc( xct, np, rhoir2(1), fxc(1, 1) )
            where( isnan( fxc(:, 1) ) ) fxc(:, 1) = 1.0_dp
            xc_kernel_ir(:, 1) = xc_kernel_ir(:, 1) + fxc(:, 1)
          case( XC_FAMILY_GGA, XC_FAMILY_HYB_GGA )
            call xc_f90_gga_vxc( xct, np, rhoir2(1), grho2(1), fxc(1, 1), fxc(1, 2) )
            where( isnan( fxc(:, 2) ) ) fxc(:, 2) = 0.0_dp
            xc_kernel_ir(:, 4) = xc_kernel_ir(:, 4) + fxc(:, 2)
            call xc_f90_gga_fxc( xct, np, rhoir2(1), grho2(1), fxc(1, 1), fxc(1, 2), fxc(1, 3) )
            where( isnan( fxc(:, 1) ) ) fxc(:, 1) = 1.0_dp
            xc_kernel_ir(:, 1) = xc_kernel_ir(:, 1) + fxc(:, 1)
            where( isnan( fxc(:, 2) ) ) fxc(:, 2) = 0.0_dp
            xc_kernel_ir(:, 2) = xc_kernel_ir(:, 2) + fxc(:, 2)
            where( isnan( fxc(:, 3) ) ) fxc(:, 3) = 0.0_dp
            xc_kernel_ir(:, 3) = xc_kernel_ir(:, 3) + fxc(:, 3)
          case default
        end select
      end do
      deallocate( grho2, rhoir2, fxc )

      ! symmetrize xc-kernel
      do i = 1, k
        call symmetrize_real_mt( xc_kernel_mt(:, :, :, i), prod_lmax, nrmt, [(j,j=1,nsymcrys)] )
        call symmetrize_real_ir( xc_kernel_ir(:, i), &
               prod_Gset%ngrtot, prod_Gset%ivg, prod_Gset%intgv, prod_Gset%ivgig, prod_Gset%igfft, &
               [(j,j=1,nsymcrys)] )
      end do

      ! **** PREPARE CONSTANT PART OF KERNEL APPLICATION
      if( gga ) then
        ! ** interstitial region
        ! d^2e/dndsigma * grad n in reciprocal space
        do i = 1, 3
          zfft(:, i) = cmplx( xc_kernel_ir(:, 2) * grhoir2(:, i), 0.0_dp, dp )
          call zfftifc( 3, prod_Gset%ngrid, -1, zfft(:, i) )
        end do
        zfft(:, 0) = zzero
        ! grad.(d^2e/dndsigma * grad n)
        do ig = 1, prod_Gset%ngvec
          if( prod_Gset%gc(ig) > dfpt_Gset%gmaxvr ) cycle
          ifg = prod_Gset%igfft(ig)
          do i = 1, 3
            zfft(ifg, 0) = zfft(ifg, 0) + prod_Gset%vgc(i, ig) * cmplx( -aimag( zfft(ifg, i) ), dble( zfft(ifg, i) ), dp )
          end do
        end do
        call zfftifc( 3, prod_Gset%ngrid, 1, zfft(:, 0) )
        ! d^2/dn^2 - 2 * grad.(d^2/dndsigma * grad n)
        xc_kernel_ir(:, 1) = xc_kernel_ir(:, 1) - 2._dp * dble( zfft(:, 0) )
        
        ! ** muffin-tin spheres
        allocate( fxc(prod_lmmax, nrmtmax) )
        allocate( gfxc(prod_lmmax, nrmtmax, 3) )
        do is = 1, nspecies
          nr = nrmt(is)
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            do i = 1, 3
              ! d^2e/dndsigma * grad n
              call dshmul( dfpt_lmaxvr, dfpt_lmaxvr, prod_lmax, nr, 1.0_dp, &
                     xc_kernel_mt(:, :, ias, 2), prod_lmmax, &
                     grhomt2(:, :, i, ias), prod_lmmax, 0.0_dp, &
                     fxc, prod_lmmax )
              ! grad.(d^2e/dndsigma * grad n)
              call gradrfmt( prod_lmax, nr, spr(1:nr, is), prod_lmmax, nrmtmax, fxc, gfxc )
              ! d^2/dn^2 - 2 * grad.(d^2/dndsigma * grad n)
              xc_kernel_mt(:, :, ias, 1) = xc_kernel_mt(:, :, ias, 1) - 2.0_dp * gfxc(:, :, i)
            end do
          end do
        end do
        deallocate( fxc, gfxc )
      end if

      deallocate( zfft )
#endif
    end subroutine gen_xc_kernel

    !> This subroutine applies the exchange-correlation kernel 
    !> to a given real valued test function \(\phi({\bf r})\) in the muffin-tin spheres.
    !>
    !> The application of the xc-kernel to a testfunction \(\phi({\bf r})\), i.e.,
    !> \[ \int f_{\rm xc}({\bf r},{\bf r}') \, \phi({\bf r}') \, {\rm d}{\bf r}' \]
    !> is given by
    !> \[ \int f^{\rm LDA}_{\rm xc}({\bf r},{\bf r}') \, \phi({\bf r}') \, {\rm d}{\bf r}' 
    !>   = \frac{\partial^2 \epsilon_{\rm xc}[n({\bf r})]}{\partial n({\bf r})^2} \phi({\bf r}) 
    !>   = f^{\rm LDA}_{\rm xc,1}({\bf r}) \, \phi({\bf r}) \]
    !> and
    !> \[\begin{align*}
    !>   \int f^{\rm GGA}_{\rm xc}({\bf r},{\bf r}') \, \phi({\bf r}') \, {\rm d}{\bf r}' 
    !>   &= \left[ \frac{\partial^2 \epsilon_{\rm xc}[n({\bf r})]}{\partial n({\bf r})^2} 
    !>      -2 \nabla \cdot \left( \frac{\partial^2 \epsilon_{\rm xc}[n({\bf r})]}
    !>      {\partial n({\bf r}) \, \partial \sigma({\bf r})} \nabla n({\bf r}) \right) \right] \phi({\bf r}) \\
    !>      &\phantom{=} -\nabla \cdot \left[ 2 \frac{\partial \epsilon_{\rm xc}[n({\bf r})]}{\partial \sigma({\bf r})}
    !>      \nabla \phi({\bf r}) + 4\frac{\partial^2 \epsilon_{\rm xc}[n({\bf r})]}{\partial \sigma({\bf r})^2}
    !>      \left( \nabla n({\bf r}) \cdot \nabla \phi({\bf r}) \right) \nabla n({\bf r}) \right] \\
    !>   &= f^{\rm LDA}_{\rm xc,1}({\bf r}) \, \phi({\bf r})
    !>      -\nabla \cdot \left[ 2 f^{\rm LDA}_{\rm xc,4}({\bf r}) \, \nabla \phi({\bf r}) 
    !>      +4 f^{\rm LDA}_{\rm xc,3}({\bf r}) \left( \nabla n({\bf r}) \cdot \nabla \phi({\bf r}) \right) \nabla n({\bf r}) \right]
    !> \end{align*}\]
    !> for LDA and GGA, respectively.
    !> See also [[gen_xc_kernel(subroutine)]].
    subroutine apply_xckernel_mt( is, ia, rfmt, res )
      use mod_atoms, only: idxas, spr
      use mod_muffin_tin, only: nrmt
      !> species and atom index
      integer, intent(in) :: is, ia
      !> real muffin-tin function as real spherical harmonics expansion
      real(dp), intent(in) :: rfmt(dfpt_lmmaxvr, *)
      !> xc-kernel applied to the function as real spherical harmonics expansion
      real(dp), intent(out) :: res(dfpt_lmmaxvr, *)

      integer :: ias, nr, i
      real(dp), allocatable :: gfmt(:,:,:), gfmt2(:,:,:), ggprod(:,:), gfxc(:,:,:)

      ias = idxas(ia, is)
      nr = nrmt(is)

      call dshmul( dfpt_lmaxvr, dfpt_lmaxvr, dfpt_lmaxvr, nr, 1.0_dp, &
             xc_kernel_mt(:, :, ias, 1), prod_lmmax, &
             rfmt, dfpt_lmmaxvr, 0.0_dp, &
             res, dfpt_lmmaxvr )

      if( gga ) then
        allocate( gfmt(dfpt_lmmaxvr, nr, 3) )
        allocate( gfmt2(prod_lmmax, nr, 0:3) )
        allocate( ggprod(prod_lmmax, nr) )
        allocate( gfxc(prod_lmmax, nr, 3) )
        ! get gradient of function
        call gradrfmt( dfpt_lmaxvr, nr, spr(1:nr, is), dfpt_lmmaxvr, nr, rfmt, gfmt )
        ! compute grad f . grad n
        ggprod = 0.0_dp
        do i = 1, 3
          call dshmul( dfpt_lmaxvr, dfpt_lmaxvr, prod_lmax, nr, 1.0_dp, &
                 gfmt(:, :, i), dfpt_lmmaxvr, &
                 grhomt2(:, :, i, ias), prod_lmmax, 1.0_dp, &
                 ggprod, prod_lmmax )
        end do
        ! compute 2 * de/dsigma * grad f + 4 * d^2e/dsigma^2 * (grad f . grad n) * grad n
        call dshmul( dfpt_lmaxvr, prod_lmax, prod_lmax, nr, 1.0_dp, &
               xc_kernel_mt(:, :, ias, 3), prod_lmmax, &
               ggprod, prod_lmmax, 0.0_dp, &
               gfmt2(:, :, 0), prod_lmmax )
        do i = 1, 3
          call dshmul( dfpt_lmaxvr, dfpt_lmaxvr, prod_lmax, nr, 2.0_dp, &
                 xc_kernel_mt(:, :, ias, 4), prod_lmmax, &
                 gfmt(:, :, i), dfpt_lmmaxvr, 0.0_dp, &
                 gfmt2(:, :, i), prod_lmmax )
          call dshmul( prod_lmax, dfpt_lmaxvr, prod_lmax, nr, 4.0_dp, &
                 gfmt2(:, :, 0), prod_lmmax, &
                 grhomt2(:, :, i, ias), prod_lmmax, 1.0_dp, &
                 gfmt2(:, :, i), prod_lmmax )
        end do
        ! compute divergence of all
        ggprod = 0.0_dp
        do i = 1, 3
          call gradrfmt( dfpt_lmaxvr, nr, spr(1:nr, is), prod_lmmax, nr, gfmt2(:, :, i), gfxc )
          ggprod = ggprod + gfxc(:, :, i)
        end do
        ! add to result
        res(:, 1:nr) = res(:, 1:nr) - ggprod(1:dfpt_lmmaxvr, :)

        deallocate( gfmt, gfmt2, ggprod, gfxc )
      end if
    end subroutine apply_xckernel_mt

    !> This subroutine applies the exchange-correlation kernel 
    !> to a given complex valued test function \(\phi({\bf r})\) in the interstitial region.
    !>
    !> See [[apply_xckernel_mt(subroutine)]] for further details.
    subroutine apply_xckernel_ir( zfir, Gpset, res )
      use constants, only: zzero
      use mod_kpointset, only: G_set
      !> complex interstitial function (possibly with Bloch wave vector \({\bf p}\))
      complex(dp), intent(in) :: zfir(dfpt_Gset%ngrtot)
      !> set of \({\bf G+p}\) vectors in which the function is expanded
      type(G_set), intent(in) :: Gpset
      !> xc-kernel applied to the function
      complex(dp), intent(inout) :: res(dfpt_Gset%ngrtot)
      
      integer :: i, ig, igp, ifg
      complex(dp), allocatable :: res2(:), zfft(:,:), zfir2(:), gfir(:,:), ggprod(:)

      allocate( res2(prod_Gset%ngrtot) )
      allocate( zfft(prod_Gset%ngrtot, 3), source=zzero )

      ! get input function on doubled grid
      allocate( zfir2(prod_Gset%ngrtot), source=zzero )
      zfft(1:dfpt_Gset%ngrtot, 1) = zfir
      call zfftifc( 3, dfpt_Gset%ngrid, -1, zfft(:, 1) )
      call dfpt_Gset%igfft2ig( zfft(:, 1), zfft(:, 2) )
      zfft(:, 1) = zzero
      call dfpt_Gset%change_set( Gpset, zfft(:, 2), zfft(:, 1), 'pull', ng=Gpset%ngvec )
      call Gpset%change_set( prod_Gset, zfft(:, 1), zfft(:, 3), 'push', ng=Gpset%ngvec )
      call prod_Gset%ig2igfft( zfft(:, 3), zfir2, gmax=Gpset%gmaxvr )
      call zfftifc( 3, prod_Gset%ngrid, 1, zfir2 )

      ! apply multiplicative part of kernel
      res2 = xc_kernel_ir(:, 1) * zfir2

      deallocate( zfir2 )

      if( gga) then
        allocate( gfir(prod_Gset%ngrtot, 3), ggprod(prod_Gset%ngrtot) )
        ! get gradient of function on doubled grid
        gfir = zzero
        do i = 1, 3
          do igp = 1, Gpset%ngvec
            ig = prod_Gset%ivgig( Gpset%ivg(1, igp), Gpset%ivg(2, igp), Gpset%ivg(3, igp) )
            ifg = prod_Gset%igfft(ig)
            gfir(ifg, i) = Gpset%vgc(i, igp) * cmplx( -aimag( zfft(igp, 1) ), dble( zfft(igp, 1) ), dp )
          end do
          call zfftifc( 3, prod_Gset%ngrid, 1, gfir(:, i) )
        end do
        ! compute grad f . grad n
        ggprod = zzero
        do i = 1, 3
          ggprod = ggprod + gfir(:, i) * grhoir2(:, i)
        end do
        ! compute 2 * de/dsigma * grad f + 4 * d^2e/dsigma^2 * (grad f . grad n) * grad n
        do i = 1, 3
          zfft(:, i) = 2.0_dp * xc_kernel_ir(:, 4) * gfir(:, i) + &
                       4.0_dp * xc_kernel_ir(:, 3) * ggprod * grhoir2(:, i)
        end do
        ! compute divergence of all
        ggprod = zzero
        do i = 1, 3
          call zfftifc( 3, prod_Gset%ngrid, -1, zfft(:, i) )
          do igp = 1, Gpset%ngvec
            ig = prod_Gset%ivgig( Gpset%ivg(1, igp), Gpset%ivg(2, igp), Gpset%ivg(3, igp) )
            ifg = prod_Gset%igfft(ig)
            ggprod(ifg) = ggprod(ifg) + Gpset%vgc(i, igp) * cmplx( -aimag( zfft(ifg, i) ), dble( zfft(ifg, i) ), dp )
          end do
        end do
        call zfftifc( 3, prod_Gset%ngrid, 1, ggprod )
        ! add to result
        res2 = res2 - ggprod
        deallocate( gfir, ggprod )
      end if

      ! get result on original grid
      zfft = zzero
      call zfftifc( 3, prod_Gset%ngrid, -1, res2 )
      call prod_Gset%igfft2ig( res2, zfft(:, 1) )
      call prod_Gset%change_set( Gpset, zfft(:, 1), zfft(:, 2), 'pull', ng=Gpset%ngvec )
      call Gpset%change_set( dfpt_Gset, zfft(:, 2), zfft(:, 3), 'push', ng=Gpset%ngvec )
      res2 = zzero
      call dfpt_Gset%ig2igfft( zfft(:, 3), res2 )
      call zfftifc( 3, dfpt_Gset%ngrid, 1, res2 )

      ! add result to output
      res = res + res2(1:dfpt_Gset%ngrtot)

      deallocate( res2, zfft )
    end subroutine apply_xckernel_ir

end module dfpt_density_potential
