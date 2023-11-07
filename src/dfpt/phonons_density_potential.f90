!> This module contains procedures for the calculation of the response
!> of the electronic density and the effective potential upon a phonon-
!> like perturbation \(\delta^{\rm q}_{I \mu}\).
module phonons_density_potential
  use dfpt_variables
  use dfpt_density_potential
  use phonons_variables

  use precision, only: dp
  use asserts, only: assert

  implicit none
  private

  ! LOCAL VARIABLES
  !> gradient of muffin-tin density
  real(dp), allocatable :: grho_mt(:,:,:,:)
  !> multipoles from gradient of electronic density
  complex(dp), allocatable :: grho_mt_mp(:,:,:)
  !> solution to Poisson's equation in muffin-tin spheres for density gradient
  complex(dp), allocatable :: gpot_coul_sph(:,:,:,:)
  !> gradient of Coulomb potential on muffin-tin sphere boundaries
  complex(dp), allocatable :: gpot_coul_surf(:,:,:)
  !> gradient of exchange correlation potential in muffin-tin spheres
  real(dp), allocatable :: gpot_xc_mt(:,:,:,:)
  !> surface integral of density discontinuity
  real(dp), allocatable :: drho_surf_int(:,:)

  public :: ph_rhopot_init, ph_rhopot_free
  public :: ph_rhopot_gen_dpot
  public :: ph_rhopot_init_drho, ph_rhopot_gen_drho_k, ph_rhopot_gen_drho_mt
  public :: ph_rhopot_symmetrize
  public :: ph_rhopot_rotate_q_canonical

  contains

    !> This subroutine initializes variables for the calculation of density
    !> and potential response that remain constant during the entire 
    !> phonon calculation.
    !>
    !> This includes:
    !> 
    !> * calculation of Coulomb potential gradient
    !> * calculation of density gradient multipole moments
    !> * calculation of exchange correlation potential gradient
    subroutine ph_rhopot_init
      use mod_potential_and_density, only: rho_mt => rhomt
      use mod_atoms, only: natmtot, nspecies, natoms, idxas, spr
      use mod_muffin_tin, only: nrmt, nrmtmax
      integer :: is, ia, ias

      ! get gradient of muffin-tin density
      allocate( grho_mt(dfpt_lmmaxvr, nrmtmax, 3, natmtot) )
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          call gradrfmt( dfpt_lmaxvr, nrmt(is), spr(1:nrmt(is), is), dfpt_lmmaxvr, nrmtmax, rho_mt(:, :, ias), grho_mt(:, :, :, ias) )
        end do
      end do

      ! generate gradient of Coulomb potential and density gradient multipole moments
      call gen_gpot_coul

      ! generate gradient of xc potential and xc kernel
      call gen_gpot_xc
    end subroutine ph_rhopot_init

    !> This subroutine frees memory from the module variables.
    subroutine ph_rhopot_free
      if( allocated( grho_mt ) ) deallocate( grho_mt )
      if( allocated( grho_mt_mp ) ) deallocate( grho_mt_mp )
      if( allocated( gpot_coul_sph ) ) deallocate( gpot_coul_sph )
      if( allocated( gpot_coul_surf ) ) deallocate( gpot_coul_surf )
      if( allocated( gpot_xc_mt ) ) deallocate( gpot_xc_mt )
      if( allocated( drho_surf_int ) ) deallocate( drho_surf_int )
    end subroutine ph_rhopot_free

    !> This subroutine calculates the effective potential response from a given
    !> soft density response using Weinert's method. See also [[dfpt_rhopot_gen_dpot(subroutine)]].
    !>
    !> This subroutine expects the soft density response \(\breve{\delta}\phantom{}^{\bf q}_{I \mu} n({\bf r})\)
    !> as input (i.e., without the gradient). Accordingly, in the muffin-tins 
    !> additional multipole moments coming from the gradient arise and must be passed 
    !> to the call of [[dfpt_rhopot_gen_dpot(subroutine)]]. Further, the homogeneous solution
    !> from the Coulomb potential gradient has to be added. See also [[gen_gpot_coul(subroutine)]].
    !>
    !> If `soft=.true.`, then only the soft part of the potential response 
    !> \(\breve{\delta}\phantom{}^{\bf q}_{I \mu} V_{\rm eff}\) is calculated.
    !> Otherwise, the full potential response
    !> \[ \delta^{\bf q}_{I \mu} V_{\rm eff}({\bf r}) 
    !>    = \breve{\delta}\phantom{}^{\bf q}_{I \mu} V_{\rm eff}({\bf r})
    !>    - \sum\limits_{\kappa,\alpha} p^{I \mu}_{\kappa\alpha}({\bf q})\, 
    !>      \nabla_\alpha \left. V_{\rm eff}({\bf r}) \right|_{{\rm MT}\kappa} \]
    !> is returned.
    !>
    !> If `coulomb=.true.`, then the Coulomb potential response is included.
    !> If `xc=.true.`, then the exchange correlation potential response is included.
    subroutine ph_rhopot_gen_dpot( pat, drho_mt, drho_ir, dpot_mt, dpot_ir, &
        soft, coulomb, xc )
      use constants, only: zzero
      use mod_atoms, only: natmtot, nspecies, natoms, idxas, spr
      use mod_muffin_tin, only: nrmtmax, rmt, nrmt
      !> displacement patterns \(p^{I \mu}_{\kappa\alpha}({\bf q})\)
      complex(dp), intent(in) :: pat(3, natmtot)
      !> soft muffin-tin density response as complex spherical harmonics expansion
      complex(dp), intent(in) :: drho_mt(:,:,:)
      !> interstitial density response on real space FFT grid
      complex(dp), intent(in) :: drho_ir(:)
      !> muffin-tin effective potential response as complex spherical harmonics expansion
      complex(dp), intent(out) :: dpot_mt(:,:,:)
      !> interstitial effective potential response on real space FFT grid
      complex(dp), intent(out) :: dpot_ir(:)
      !> compute only the soft effective potential response (default: `.true.`)
      logical, optional, intent(in) :: soft
      !> include the Coulomb / exchange correlation potential response (default: `.true.`)
      logical, optional, intent(in) :: coulomb, xc

      integer :: is, ia, ias, ip, l, m, lm, ir
      real(dp) :: rl(nrmtmax, 0:dfpt_lmaxvr)
      logical :: dosoft, vcoul, vxc

      complex(dp), allocatable :: mt_mp(:,:)

      dosoft = .true.
      if( present( soft ) ) dosoft = soft
      vcoul = .true.
      if( present( coulomb ) ) vcoul = coulomb
      vxc = .true.
      if( present( xc ) ) vxc = xc

      ! set additional multipole moments coming from density gradient
      allocate( mt_mp(dfpt_lmmaxvr, natmtot), source=zzero )
      do ias = 1, natmtot
        do ip = 1, 3
          mt_mp(:, ias) = mt_mp(:, ias) - pat(ip, ias) * grho_mt_mp(:, ias, ip)
        end do
      end do

      call dfpt_rhopot_gen_dpot( drho_mt, drho_ir, dpot_mt, dpot_ir, &
             ph_Gqset, ph_jlgqr, ph_ylmgq, ph_sfacgq, &
             coulomb=vcoul, xc=vxc, mt_mp_add=mt_mp )

      ! add homogeneous solution of Poisson's equation in MT spheres
      do is = 1, nspecies
        rl(:, 0) = 1.0_dp
        rl(1:nrmt(is), 1) = spr(1:nrmt(is), is) / rmt(is)
        do l = 2, dfpt_lmaxvr
          rl(:, l) = rl(:, l-1) * rl(:, 1)
        end do

        do ia = 1, natoms(is)
          ias = idxas(ia, is)
!$omp parallel default( shared ) private( ir, ip, l, m, lm, mt_mp )
!$omp do
          do ir = 1, nrmt(is)
            do ip = 1, 3
              ! for soft response add homogeneous solution from potential gradient for displaced atom
              if( dosoft ) then
                if( vcoul ) then
                  do l = 0, dfpt_lmaxvr
                    do m = -l, l
                      lm = l * (l + 1) + m + 1
                      dpot_mt(lm, ir, ias) = dpot_mt(lm, ir, ias) + pat(ip, ias) * rl(ir, l) * gpot_coul_surf(lm, ias, ip)
                    end do
                  end do
                end if
              ! for full response subtract gradient of potential for displaced atom
              else
                if( vcoul ) then
                  do l = 0, dfpt_lmaxvr
                    do m = -l, l
                      lm = l * (l + 1) + m + 1
                      dpot_mt(lm, ir, ias) = dpot_mt(lm, ir, ias) - pat(ip, ias) * gpot_coul_sph(lm, ir, ias, ip)
                    end do
                  end do
                end if
                if( vxc ) then
                  call rtozflm( dfpt_lmaxvr, gpot_xc_mt(:, ir, ip, ias), mt_mp(:, 1) )
                  dpot_mt(:, ir, ias) = dpot_mt(:, ir, ias) - pat(ip, ias) * mt_mp(:, 1)
                end if
              end if
            end do
          end do
!$omp end do
!$omp end parallel

        end do
      end do

      deallocate( mt_mp )
    end subroutine ph_rhopot_gen_dpot

    !> This subroutine initializes the density response for the self-consistency
    !> cycle.
    !> 
    !> Therefore, the density response is set to the density gradient, i.e.,
    !> in the muffin-tin spheres it is set to
    !> \[ \delta^{\bf q}_{I \mu} n({\bf r}_\kappa) 
    !>    = \sum_\alpha p^{I \mu}_{\kappa\alpha}({\bf q})\, \nabla_\alpha n({\bf r}_\kappa) \]
    !> and in the interstitial region it is set to
    !> \[ \delta^{\bf q}_{I \mu} n({\bf r}) 
    !>    = \frac{1}{N_{\rm at}}\sum\limits_{\kappa,\alpha} 
    !>      p^{I \mu}_{\kappa\alpha}({\bf q})\, \nabla_\alpha n({\bf r}) \;. \]
    !>
    !> If `soft=.true.`, then only the soft density response 
    !> \(\breve{\delta}\phantom{}^{\bf q}_{I \mu} n\) is initialized. It differs from
    !> \(\delta^{\bf q}_{I \mu} n\) in that it is zero inside the muffin-tin spheres.
    subroutine ph_rhopot_init_drho( pat, drho_mt, drho_ir, &
        soft )
      use constants, only: zzero
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use mod_muffin_tin, only: nrmt
      !> displacement patterns \(p^{I \mu}_{\kappa\alpha}({\bf q})\)
      complex(dp), intent(in) :: pat(3, natmtot)
      !> (soft) muffin-tin density response as complex spherical harmonics expansion
      complex(dp), intent(out) :: drho_mt(:,:,:)
      !> interstitial density response on real space FFT grid
      complex(dp), intent(out) :: drho_ir(:)
      !> compute only the soft density response (default: `.true.`)
      logical, optional, intent(in) :: soft

      integer :: is, ia, ias, ip, ir
      logical :: dosoft

      complex(dp), allocatable :: zfft(:)

      dosoft = .true.
      if( present( soft ) ) dosoft = soft

      ! initialize interstitial density response to density gradient
      drho_ir = zzero

      ! initialize muffin-tin density response
      drho_mt = zzero
      if( .not. dosoft ) then
        allocate( zfft(dfpt_lmmaxvr) )
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            do ip = 1, 3
              do ir = 1, nrmt(is)
                call rtozflm( dfpt_lmaxvr, grho_mt(:, ir, ip, ias), zfft )
                drho_mt(:, ir, ias) = drho_mt(:, ir, ias) - pat(ip, ias) * zfft
              end do
            end do
          end do
        end do
        deallocate( zfft )
      end if
    end subroutine ph_rhopot_init_drho

    !> This subroutine calculates the contribution to the electron density response
    !> coming from a single \({\bf k}\) point. See also [[dfpt_rhopot_drho_k(subroutine)]].
    !> 
    !> It differs from [[dfpt_rhopot_drho_k(subroutine)]] in that it accounts for the
    !> response of the basis functions coming from the perturbed matching coefficients
    !> \(\delta^{\bf q}_{I \mu} A^\kappa_{{\bf G+k},lm,\xi}\) (see [[gen_dapwalm(subroutine)]]).
    subroutine ph_rhopot_gen_drho_k( ik, kset, Gkset, Gkqset, fst, lst, occk, docck, eveck, deveck, apwalmk, apwalmkq, pat, gamma, drho_mat, drho_ir )
      use phonons_eigensystem, only: gen_dapwalm

      use constants, only: zzero, zone
      use mod_kpointset, only: k_set, Gk_set
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use mod_APW_LO, only: nlotot
      use mod_eigensystem, only: idxlo
      use mod_lattice, only: omega
      use modinput
      !> index of the wavevector \({\bf k}\) in the set
      integer, intent(in) :: ik
      !> set of \({\bf k}\) vectors
      type(k_set), intent(in) :: kset
      !> set of \({\bf G+k}\) vectors for the wavefunctions
      type(Gk_set), intent(in) :: Gkset
      !> set of \({\bf G+k+q}\) vectors for the wavefunction response
      type(Gk_set), intent(in) :: Gkqset
      !> first and last state to consider
      integer, intent(in) :: fst, lst
      !> occupation numbers at \({\bf k}\)
      real(dp), intent(in) :: occk(:)
      !> occupation number response at \({\bf k}\)
      real(dp), intent( in) :: docck(:)
      !> eigenvectors at \({\bf k}\)
      complex(dp), intent(in) :: eveck(:,:)
      !> eigenvector response at \({\bf k}\)
      complex(dp), intent(in) :: deveck(:,:)
      !> (L)APW matching coefficients \(A^\kappa_{{\bf G+p},lm,\xi}\) at \({\bf k}\) and \({\bf k+q}\)
      complex(dp), intent(in) :: apwalmk(:,:,:,:), apwalmkq(:,:,:,:)
      !> displacement patterns \(p^{I \mu}_{\kappa\alpha}({\bf q})\)
      complex(dp), intent(in) :: pat(3, natmtot)
      !> `.true.` for Gamma point phonons
      logical, intent(in) :: gamma
      !> muffin-tin density response matrix
      complex(dp), intent(inout) :: drho_mat(:,:,:)
      !> interstitial density response on real space FFT grid
      complex(dp), intent(inout) :: drho_ir(:)

      integer :: ngk, ngkq, i, ist, nst
      integer :: is, ia, ias
      real(dp) :: t1, t2
      logical :: unpert

      integer, allocatable :: bands(:)
      real(dp), allocatable :: wgt1(:), wgt2(:)
      complex(dp), allocatable :: dapwalmk(:,:,:)
      complex(dp), allocatable :: evecsum(:,:), evecmt1(:,:), evecmt2(:,:), devecmt(:,:)
      complex(dp), allocatable :: zfft(:,:)

      ngk = Gkset%ngk(1, ik)
      ngkq = Gkqset%ngk(1, ik)
      unpert = (sum( abs( pat ) ) < 1e-12_dp)

      ! find k-point weights for bands and bands that contribute
      allocate( wgt1(fst:lst), wgt2(fst:lst) )
      wgt1 = [(2*kset%wkpt(ik)*occk(i), i=fst, lst)]
      wgt2 = [(kset%wkpt(ik)*docck(i), i=fst, lst)]
      if( gamma ) then
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
      allocate( dapwalmk, source=apwalmk(:, :, :, 1) )
      allocate( evecsum(mt_basis%n_basis_fun_max, mt_basis%n_basis_fun_max) )
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          evecsum = zzero
          ! generate shortened MT eigenvector response at current k-point
          call mt_basis%transform_evec( is, ngkq, nlotot, idxlo(:, :, ias), apwalmkq(:, :, :, ias), deveck(:, bands), devecmt )
          ! generate shortened MT eigenvector at current k-point
          call mt_basis%transform_evec( is, ngk, nlotot, idxlo(:, :, ias), apwalmk(:, :, :, ias), eveck(:, bands), evecmt1 )
          ! contribution from occupation response
          if( gamma ) then
            allocate( evecmt2, source=evecmt1 )
            do i = 1, nst
              ist = bands(i)
              evecmt2(:, i) = wgt2(ist) * evecmt1(:, i)
            end do
            call zgemm( 'n', 'c', mt_basis%n_basis_fun(is), mt_basis%n_basis_fun(is), nst, zone, &
                   evecmt2, size(evecmt2, dim=1), &
                   evecmt1, size(evecmt1, dim=1), zone, &
                   evecsum, mt_basis%n_basis_fun_max )
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
                 evecsum, mt_basis%n_basis_fun_max )
          if( sum( abs( pat(:, ias) ) ) > 1e-12_dp ) then
            ! get perturbed matching coefficients at k
            call gen_dapwalm( ngk, Gkset%vgkc(:, :, 1, ik), pat(:, ias), apwalmk(:, :, :, ias), dapwalmk )
            call mt_basis%transform_evec( is, ngk, nlotot, idxlo(:, :, ias), dapwalmk, eveck(:, bands), evecmt2, &
                   use_local_orbitals=.false. )
            call zgemm( 'n', 'c', mt_basis%n_basis_fun(is), mt_basis%n_basis_fun(is), nst, zone, &
                   evecmt2, size(evecmt2, dim=1), &
                   evecmt1, size(evecmt1, dim=1), zone, &
                   evecsum, mt_basis%n_basis_fun_max )
            if( allocated( evecmt2 ) ) deallocate( evecmt2 )
          end if
          drho_mat(:, :, ias) = drho_mat(:, :, ias) + evecsum(:, :)
        end do
      end do
      if( allocated( evecmt1 ) ) deallocate( evecmt1 )
      if( allocated( devecmt ) ) deallocate( devecmt )
      deallocate( dapwalmk, evecsum )

      ! **** interstitial density response
      ! sum over occupied states
      allocate( zfft(dfpt_Gset%ngrtot, 2) )
      do i = 1, nst
        ist = bands(i)
        t1 = wgt1(ist) / omega
        t2 = wgt2(ist) / omega
        zfft = zzero
        ! Fourier transform wavefunction to real space
        zfft(dfpt_Gset%igfft(Gkset%igkig(1:ngk, 1, ik)), 1) = eveck(1:ngk, ist)
        call zfftifc( 3, dfpt_Gset%ngrid, 1, zfft(:, 1) )
        ! Fourier transform wavefunction response to real space
        zfft(dfpt_Gset%igfft(Gkqset%igkig(1:ngkq, 1, ik)), 2) = deveck(1:ngkq, ist)
        call zfftifc( 3, dfpt_Gset%ngrid, 1, zfft(:, 2) )
        ! multiply wavefunction and wavefunction response and add to density response
        if( gamma ) then
          drho_ir = drho_ir + conjg( zfft(:, 1) ) * (t1 * zfft(:, 2) + t2 * zfft(:, 1))
        else
          drho_ir = drho_ir + t1 * conjg( zfft(:, 1) ) * zfft(:, 2)
        end if
      end do
      deallocate( zfft )
      deallocate( wgt1, wgt2 )
    end subroutine ph_rhopot_gen_drho_k

    !> This subroutine calculates the muffin-tin density response 
    !> from the given density response matrix \(\delta^{\bf q}_{I \mu}{\bf D}^\kappa\).
    !> See also [[dfpt_rhopot_gen_drho_mt(subroutine)]].
    !>
    !> If `soft=.true.`, then only the soft part of the density response 
    !> \(\breve{\delta}\phantom{}^{\bf q}_{I \mu} n\) is calculated.
    !> Otherwise, the full density response
    !> \[ \delta^{\bf q}_{I \mu} n({\bf r}) 
    !>    = \breve{\delta}\phantom{}^{\bf q}_{I \mu} n({\bf r})
    !>    - \sum\limits_{\kappa,\alpha} p^{I \mu}_{\kappa\alpha}({\bf q})\, 
    !>      \nabla_\alpha \left. n({\bf r}) \right|_{{\rm MT}\kappa} \]
    !> is returned.
    subroutine ph_rhopot_gen_drho_mt( pat, drho_mat, drho_mt, &
        soft )
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use mod_muffin_tin, only: nrmt
      !> displacement patterns \(p^{I \mu}_{\kappa\alpha}({\bf q})\)
      complex(dp), intent(in) :: pat(3, natmtot)
      !> muffin-tin density response matrix
      complex(dp), intent(in) :: drho_mat(:,:,:)
      !> (soft) muffin-tin density response as complex spherical harmonics expansion
      complex(dp), intent(out) :: drho_mt(:,:,:)
      !> compute only the soft density response (default: `.true.`)
      logical, optional, intent(in) :: soft

      integer :: is, ia, ias, ip, ir
      complex(dp) :: zflm(dfpt_lmmaxvr)
      logical :: dosoft

      dosoft = .true.
      if( present( soft ) ) dosoft = soft

      call dfpt_rhopot_gen_drho_mt( drho_mat, drho_mt )
      
      ! for full density response substract muffin tin gradient
      if( .not. dosoft ) then
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            do ip = 1, 3
              do ir = 1, nrmt(is)
                call rtozflm( dfpt_lmaxvr, grho_mt(:, ir, ip, ias), zflm )
                drho_mt(:, ir, ias) = drho_mt(:, ir, ias) - pat(ip, ias) * zflm
              end do
            end do
          end do
        end do
      end if
    end subroutine ph_rhopot_gen_drho_mt

    !> This subroutine symmetrizes a phonon-like perturbed unit cell function
    !> in the basis of an irrep. See also [[phonons_symmetry(module)]].
    !>
    !> The symmetrized function is given by
    !> \[ \delta^{\bf q}_{I \mu} f({\bf r})
    !>    = \frac{1}{N_S} \sum\limits_\mathcal{S}^{\mathcal{G}_{\bf q}} \sum_{\nu=1}^{d_I}
    !>      {\rm e}^{-{\rm i}{\bf q}\cdot{\bf \tau}_S}\, \texttt{S}^I_{\nu\mu} \,
    !>      \delta^{\bf q}_{I \nu} \bar{f}(\mathcal{S}^{-1}{\bf r}) \;, \]
    !> where \(\delta^{\bf q}_{I \nu} \bar{f}\) is the unsymmetrized function.
    subroutine ph_rhopot_symmetrize( zfun_mt, zfun_ir, vql, dirrep, nsym, isym, ivsym, symmat )
      use constants, only: zzero, twopi
      use mod_symmetry, only: lsplsymc, symlat, symlatc, ieqatom, vtlsymc, symapp_zfig
      use mod_atoms, only: natmtot, natmmax, nspecies, natoms, idxas
      use mod_muffin_tin, only: nrmtmax, nrmt
      use modinput
      !> complex muffin-tin function for each irrep member as complex spherical harmonics expansion
      complex(dp), intent(inout) :: zfun_mt(dfpt_lmmaxvr, nrmtmax, natmtot, *)
      !> complex interstitial function for each irrep member on real space FFT grid
      complex(dp), intent(inout) :: zfun_ir(dfpt_Gset%ngrtot, *)
      !> Bloch wavevector \({\bf q}\) in lattice coordinates
      real(dp), intent(in) :: vql(3)
      !> dimension \(d_I\) of the irrep
      integer, intent(in) :: dirrep
      !> number of symmetry operations in small group of \({\bf q}\)
      integer, intent(in) :: nsym
      !> indices of symmetries in global arrays
      integer, intent(in) :: isym(:)
      !> lattice vectors that map \({\bf \rm S}{\bf q}\) back to 1st BZ
      integer, intent(in) :: ivsym(3, *)
      !> matrix representation of symmetries in the basis of the irrep \({\bf \texttt{S}}^I\)
      complex(dp), intent(in) :: symmat(:,:,:)

      integer :: is, ia, ias, ja, jas, d, dd, s, lspl
      real(dp) :: sl(3,3), sc(3,3), a(3), b(3), phi

      complex(dp), allocatable :: zfun_mt1(:,:,:,:), zfun_mt2(:,:)
      complex(dp), allocatable :: zfun_ig1(:,:), zfun_ig2(:)

      ! **** muffin-tin part
      allocate( zfun_mt1(dfpt_lmmaxvr, nrmtmax, natmmax, dirrep) )
      allocate( zfun_mt2(dfpt_lmmaxvr, nrmtmax) )
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! make a copy of the input function
          do d = 1, dirrep
            zfun_mt1(:, :, ia, d) = zfun_mt(:, : ,ias, d)
            zfun_mt(:, :, ias, d) = zzero
          end do
        end do
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! loop over symmetries
          do s = 1, nsym
            lspl = lsplsymc(isym(s))                ! index of spatial rotation     
            sl = dble( symlat(:, :, lspl) )         ! lattice spatial rotation
            sc = symlatc(:, :, lspl)                ! Cartesian spatial rotation
            ja = ieqatom(ia, is, isym(s))           ! equivalent atom
            jas = idxas(ja, is)

            a = input%structure%speciesarray(is)%species%atomarray(ja)%atom%coord
            call r3mv( sl, a, b )
            b = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord - b
            phi = twopi * dot_product( vql, b )
            do d = 1, dirrep
              ! apply rotation to muffin-tin functions
              call rotzflm( sc, dfpt_lmaxvr, nrmt(is), dfpt_lmmaxvr, zfun_mt1(:, :, ja, d), zfun_mt2 )
              ! add contribution to input array
              do dd = 1, dirrep
                zfun_mt(:, :, ias, dd) = zfun_mt(:, :, ias, dd) + conjg( symmat(dd, d, s) ) * cmplx( cos( phi ), sin( phi ), dp ) * zfun_mt2
              end do
            end do
          end do
        end do
      end do
      ! normalize
      do d = 1, dirrep
        zfun_mt(:, :, :, d) = zfun_mt(:, :, :, d) / nsym
      end do
      deallocate( zfun_mt1, zfun_mt2 )

      ! **** interstitial part
      ! transform function to reciprocal space
      allocate( zfun_ig1(dfpt_Gset%ngrtot, dirrep), source=zzero )
      allocate( zfun_ig2(dfpt_Gset%ngrtot) )
      do d = 1, dirrep
        call zfftifc( 3, dfpt_Gset%ngrid, -1, zfun_ir(:, d) )
        call dfpt_Gset%igfft2ig( zfun_ir(:, d), zfun_ig2 )
        call dfpt_Gset%change_set( ph_Gqset, zfun_ig2, zfun_ig1(:, d), 'pull', ng=ph_Gqset%ngvec )
        zfun_ir(:, d) = zzero
      end do
      do s = 1, nsym
        lspl = lsplsymc(isym(s))
        phi = twopi * dot_product( vql + dble( ivsym(:, s) ), vtlsymc(:, isym(s)) ) ! phase
        do d = 1, dirrep
          ! apply symmetry to interstitial function
          zfun_ig2 = zzero
          call symapp_zfig( symlat(:, :, lspl), vtlsymc(:, isym(s)), vql, &
                 zfun_ig1(:, d), ph_Gqset%ngvec, ph_Gqset%ivg, ph_Gqset%igfft, .false., &
                 zfun_ig2, dfpt_Gset%intgv, dfpt_Gset%ivgig, dfpt_Gset%igfft, .true. )
          ! add contribution to input array
          do dd = 1, dirrep
            zfun_ir(:, dd) = zfun_ir(:, dd) + conjg( symmat(dd, d, s) ) * cmplx( cos( phi ), sin( phi ), dp) * zfun_ig2
          end do
        end do
      end do
      do d = 1, dirrep
        ! transform function to real space
        call zfftifc( 3, dfpt_Gset%ngrid, 1, zfun_ir(:, d) )
        ! normalize
        zfun_ir(:, d) = zfun_ir(:, d) / nsym
      end do
      deallocate( zfun_ig1, zfun_ig2 )
    end subroutine ph_rhopot_symmetrize

    !> This subroutine ccomputes the gradient of the Coulomb potential as well as
    !> the multipole moments of the density gradient.
    !>
    !> In the case of a canonical phonon-like perturbation \(\delta^{\bf q}_{\kappa\alpha}\)
    !> the response of the Hartree potential reads
    !> \[ \delta^{\bf q}_{\kappa\alpha} V_{\rm H}({\bf r})
    !>    = \int \frac{\breve{\delta}\phantom{}^{\bf q}_{\kappa\alpha}n({\bf r}')}
    !>      {|{\bf r} - {\bf r'}|} {\rm d}{\bf r}'
    !>    - \int\limits_{{\rm MT}\kappa} \frac{\nabla_\alpha n({\bf r}')}{|{\bf r} - {\bf r}'|} {\rm d}{\bf r}'
    !>    + \sum_{\bf R} {\rm e}^{{\rm i}{\bf q}\cdot{\bf R}} \oint\limits_{\partial{\rm MT}\kappa{\bf R}}
    !>      \frac{[n({\bf r}')]_{\rm SF}}{|{\bf r} - {\bf r}'|} \hat{e}_\alpha\, {\rm d}S \;, \]
    !> where \(\breve{\delta}\phantom{}^{\bf q}_{\kappa\alpha} n\) is the soft density response
    !> and the last integral comes from the density discontinuity on the muffin-tin surface
    !> where \([\dots]_{\rm SF}\) describes the difference of the muffin-tin and interstitial
    !> representation of the integrand.
    !>
    !> The response of the external potential is given by its negative gradient
    !> \[ \delta^{\bf q}_{\kappa\alpha} V_{\rm ext}({\bf r})
    !>    = \sum_{\bf R} {\rm e}^{{\rm i}{\bf q}\cdot{\bf R}}
    !>      \nabla_\alpha \frac{Z_\kappa}{|{\bf r}-{\bf \tau}_\kappa-{\bf R}|} \;. \]
    !>
    !> Both the gradient and the surface terms give rise to additional multipole moments
    !> when calculating the Coulomb potential response using Weinert's method which are 
    !> calculated by this routine and stored in `grho_mt_mp`.
    !>
    !> Further, this routine computes the solution to Poisson's equation in the 
    !> muffin-tin spheres from the density gradient (`gpot_coul_sph`) and the 
    !> gradient of the Coulomb potential on the muffin-tin surfaces (`gpot_coul_surf`).
    subroutine gen_gpot_coul
      use weinert
      use constants, only: zzero, fourpi, y00, zi
      use gaunt, only: gaunt_coeff_yyy, non_zero_gaunt_real
      use mod_potential_and_density, only: rho_ir => rhoir, rho_mt => rhomt
      use mod_atoms, only: natmtot, nspecies, natoms, idxas, spr, spzn
      use mod_muffin_tin, only: nrmtmax, rmt, nrmt
      use modinput

      integer :: ip, ig, igf
      integer :: is, ia, ias, ir, i, l, m, lm1, lm2, lm3
      real(dp) :: t1, rl(0:dfpt_lmaxvr)
      complex(dp) :: z1
      type(non_zero_gaunt_real), pointer :: yyy

      real(dp), allocatable :: pot_ext(:)
      complex(dp), allocatable :: grho_i_mp(:,:,:), drho_surf(:,:)
      complex(dp), allocatable :: zrho_ig(:), zrho_mt(:,:)
      complex(dp), allocatable :: gzrho_ig(:,:), gzrho_mt(:,:,:)
      
      allocate( grho_mt_mp(dfpt_lmmaxvr, natmtot, 3) )
      allocate( gpot_coul_sph(dfpt_lmmaxvr, nrmtmax, natmtot, 3) )
      allocate( gpot_coul_surf(dfpt_lmmaxvr, natmtot, 3) )
      allocate( drho_surf_int(natmtot, 3) )

      ! **** initialize auxilliary variables
      allocate( grho_i_mp(dfpt_lmmaxvr, natmtot, 3) )
      allocate( drho_surf(dfpt_lmmaxvr, natmtot) )

      ! **** find interstitial multipole moments
      allocate( zrho_ig(dfpt_Gset%ngrtot) )
      allocate( gzrho_ig(dfpt_Gset%ngrtot, 3) )
      ! Fourier transform density to reciprocal space
      zrho_ig = cmplx( rho_ir, 0.0_dp, dp )      
      call zfftifc( 3, dfpt_Gset%ngrid, -1, zrho_ig )
      gzrho_ig = zzero
      do ip = 1, 3
        ! buid interstital gradient of density
        do ig = 1, dfpt_Gset%ngvec
          igf = dfpt_Gset%igfft(ig)
          gzrho_ig(igf, ip) = dfpt_Gset%vgc(ip, ig) * cmplx( -aimag( zrho_ig(igf)) , dble( zrho_ig(igf) ), dp )
        end do
        ! find interstitial multipole moments of gradient
        call multipoles_ir( dfpt_lmaxvr, ph_Gqset%ngvec, ph_Gqset%gc, ph_Gqset%ivg, &
               ph_jlgqr, ph_ylmgq, ph_sfacgq, &
               dfpt_Gset%intgv, dfpt_Gset%ivgig, dfpt_Gset%igfft, gzrho_ig(:, ip), &
               grho_i_mp(:, :, ip) )
      end do
      ! find interstitial density on MT boundaries
      call surface_ir( dfpt_lmaxvr, ph_Gqset%ngvec, ph_Gqset%gc, ph_Gqset%ivg, &
             ph_jlgqr, ph_ylmgq, ph_sfacgq, &
             dfpt_Gset%intgv, dfpt_Gset%ivgig, dfpt_Gset%igfft, zrho_ig, &
             drho_surf )

      ! **** find the muffin tin multipole moments
      allocate( zrho_mt(dfpt_lmmaxvr, nrmtmax) )
      allocate( gzrho_mt(dfpt_lmmaxvr, nrmtmax, 3) )
      do is = 1, nspecies
        t1 = sqrt( fourpi / 3.0_dp ) * rmt(is)**2
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! convert real to complex spherical harmonic expansion
          do ir = 1, nrmt(is)
            call rtozflm( dfpt_lmaxvr, rho_mt(:, ir, ias), zrho_mt(:, ir) )
          end do
          ! store MT surface density discontinuity
          drho_surf(:, ias) = zrho_mt(:, nrmt(is)) - drho_surf(:, ias)
          ! store the surface integral of the density discontinuity
          drho_surf_int(ias, 1) = t1 * (dble( drho_surf(2, ias) ) - dble( drho_surf(4, ias) )) / sqrt( 2.0_dp )
          drho_surf_int(ias, 2) = t1 * (aimag( drho_surf(2, ias) ) + aimag( drho_surf(4, ias) )) / sqrt( 2.0_dp )
          drho_surf_int(ias, 3) = t1 *dble( drho_surf(3, ias) )
          ! build MT gradient of density
          call gradzfmt( dfpt_lmaxvr, nrmt(is), spr(1:nrmt(is), is), dfpt_lmmaxvr, nrmtmax, zrho_mt, gzrho_mt) 
          ! find MT multipole moments of gradient and solution to Poisson's equation in MT spheres
          do ip = 1, 3
            call poisson_and_multipoles_mt( dfpt_lmaxvr, nrmt(is), spr(1:nrmt(is), is), gzrho_mt(:, :, ip), &
                   gpot_coul_sph(:, :, ias, ip), grho_mt_mp(:, ias, ip) )
          end do
        end do
      end do

      ! **** add ionic potential gradient in MT spheres
      allocate( pot_ext(nrmtmax) )
      zrho_mt = zzero; gzrho_mt = zzero
      do is = 1, nspecies
        call potnucl( input%groundstate%ptnucl, nrmt(is), spr(1:nrmt(is), is), spzn(is), pot_ext )
        zrho_mt(1, :) = cmplx( pot_ext(:) / y00, 0, dp )
        call gradzfmt( 1, nrmt(is), spr(1:nrmt(is), is), dfpt_lmmaxvr, nrmtmax, zrho_mt, gzrho_mt )
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          do ip = 1, 3
            do m = -1, 1
              lm1 = 3 + m
              do ir = 1, nrmt(is)
                gpot_coul_sph(lm1, ir, ias, ip) = gpot_coul_sph(lm1, ir, ias, ip) + gzrho_mt(lm1, ir, ip) - gzrho_mt(lm1, nrmt(is), ip) * spr(ir, is) / rmt(is)
              end do
              ! if nuclei are not point-like, we take the numeric multipole moments
              if( .not. input%groundstate%ptnucl ) &
                grho_mt_mp(lm1, ias, ip) = grho_mt_mp(lm1, ias, ip) + 3.d0 * rmt(is)**2 / fourpi * gzrho_mt(lm1, nrmt(is), ip)
            end do
          end do
        end do
      end do
      deallocate( pot_ext )
      deallocate( zrho_mt, gzrho_mt )

      ! **** add multipoles of ionic density response and surface contribution
      yyy => gaunt_coeff_yyy
      t1 = sqrt( fourpi / 3.0_dp )
      do is = 1, nspecies
        rl(0) = 1.0_dp
        do l = 1, dfpt_lmaxvr
          rl(l) = rl(l-1) * rmt(is)
        end do
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! surface contribution
          do l = 0, dfpt_lmaxvr
            do m = -l, l
              lm1 = l * (l + 1) + m + 1

              lm3 = 3 !(l,m) = (1,0)
              do i = 1, yyy%num(lm1, lm3)
                lm2 = yyy%lm2(i, lm1, lm3)
                if( lm2 > dfpt_lmmaxvr ) exit
                z1 = drho_surf(lm2, ias) * rl(l) * t1 * yyy%val(i, lm1, lm3)
                grho_mt_mp(lm1, ias, 3) = grho_mt_mp(lm1, ias, 3) - z1
              end do
              lm3 = 2 !(l,m) = (1,-1)
              do i = 1, yyy%num(lm1, lm3)
                lm2 = yyy%lm2(i, lm1, lm3)
                if( lm2 > dfpt_lmmaxvr ) exit
                z1 = drho_surf(lm2, ias) * rl(l) * t1 * yyy%val(i, lm1, lm3) / sqrt( 2.0_dp )
                grho_mt_mp(lm1, ias, 1) = grho_mt_mp(lm1, ias, 1) - z1
                grho_mt_mp(lm1, ias, 2) = grho_mt_mp(lm1, ias, 2) - z1 * zi
              end do
              lm3 = 4 !(l,m) = (1,1)
              do i = 1, yyy%num(lm1, lm3)
                lm2 = yyy%lm2(i, lm1, lm3)
                if( lm2 > dfpt_lmmaxvr ) exit
                z1 = drho_surf(lm2, ias) * rl(l) * t1 * yyy%val(i, lm1, lm3) / sqrt( 2.0_dp )
                grho_mt_mp(lm1, ias, 1) = grho_mt_mp(lm1, ias, 1) + z1
                grho_mt_mp(lm1, ias, 2) = grho_mt_mp(lm1, ias, 2) - z1 * zi
              end do

            end do
          end do

          ! ionic contribution
          ! if nuclei are point-like, we take the exact multipole moments
          if( input%groundstate%ptnucl ) then
            z1 = -spzn(is) / t1
            lm1 = 3 ! (l,m) = (1,0)
            grho_mt_mp(lm1, ias, 3) = grho_mt_mp(lm1, ias, 3) + z1
            z1 = z1 / sqrt( 2.0_dp )
            lm1 = 2 ! (l,m) = (1,-1)
            grho_mt_mp(lm1, ias, 1) = grho_mt_mp(lm1, ias, 1) + z1
            grho_mt_mp(lm1, ias, 2) = grho_mt_mp(lm1, ias, 2) + z1 * zi
            lm1 = 4 ! (l,m) = (1,1)
            grho_mt_mp(lm1, ias, 1) = grho_mt_mp(lm1, ias, 1) - z1
            grho_mt_mp(lm1, ias, 2) = grho_mt_mp(lm1, ias, 2) + z1 * zi
          end if

        end do
      end do

      ! **** find interstitial gradient of Coulomb potential
      do ip = 1, 3
        ! find Fourier coefficient from multipole moments
        call poisson_ir( dfpt_lmaxvr, input%groundstate%npsden, ph_Gqset%ngvec, ph_Gqset%gc, ph_Gqset%ivg, &
               ph_jlgqr, ph_ylmgq, ph_sfacgq, dfpt_Gset%intgv, dfpt_Gset%ivgig, dfpt_Gset%igfft, &
               gzrho_ig(:, ip), grho_mt_mp(:, :, ip)-grho_i_mp(:, :, ip), zrho_ig )
        ! find potential gradient on MT sphere boundaries
        call surface_ir( dfpt_lmaxvr, ph_Gqset%ngvec, ph_Gqset%gc, ph_Gqset%ivg, &
               ph_jlgqr, ph_ylmgq, ph_sfacgq, &
               dfpt_Gset%intgv, dfpt_Gset%ivgig, dfpt_Gset%igfft, zrho_ig, &
               gpot_coul_surf(:, :, ip) )
      end do
      deallocate( zrho_ig, gzrho_ig )

      ! delete auxilliary variables
      deallocate( grho_i_mp, drho_surf )
    end subroutine gen_gpot_coul

    !> This subroutine computes the gradient of the exchange correlation potential
    !> in the muffin-tin spheres and stores them in `gpot_xc_mt`.
    subroutine gen_gpot_xc
      !use dfpt_density_potential, only: apply_xckernel_mt
      use mod_potential_and_density, only: pot_xc_mt => vxcmt
      use mod_atoms, only: natmtot, nspecies, natoms, idxas, spr
      use mod_muffin_tin, only: nrmtmax, nrmt

      integer :: is, ia, ias, nr

      ! compute gradient of xc-potential in muffin-tin spheres
      allocate( gpot_xc_mt(dfpt_lmmaxvr, nrmtmax, 3, natmtot) )
      do is = 1, nspecies
        nr = nrmt(is)
        do ia = 1, natoms(is)
          ias = idxas( ia, is)
          !TODO test with kernel
          !do i = 1, 3
          !  call apply_xckernel_mt( is, ia, grho_mt(:,:,i,ias), gpot_xc_mt(:,:,i,ias))
          !end do
          ! fallback
          call gradrfmt( dfpt_lmaxvr, nr, spr(1:nr, is), dfpt_lmmaxvr, nrmtmax, pot_xc_mt(:, :, ias), gpot_xc_mt(:, :, :, ias) )
        end do
      end do
    end subroutine gen_gpot_xc

    !> This subroutine rotates the canonical response of a cell-periodic function (e.g. density and potential) 
    !> at phonon wavevector \({\bf q}_0\) into the symmetry equivalent wavevector \({\bf q}\).
    !>
    !> The function is rotated using the symmetry operation \(\mathcal{S}=\lbrace \mathrm{\bf R} | {\bf \tau} \rbrace\)
    !> that connects \({\bf q}\) and \({\bf q}_0\) via \({\bf q} = \mathrm{\bf R}({\bf q}_0 + {\bf G})\).
    !> The rotated function is given by
    !> \[ \delta^{\bf q}_{\kappa\,\alpha} f({\bf r}) = 
    !>      \sum_{\kappa',\beta} \delta^{{\bf q}_0}_{\kappa'\,\beta} f(\mathcal{S}^{-1}{\bf r})\,
    !>      \Gamma^\ast_{\kappa\,\alpha,\kappa'\,\beta}(\mathcal{S};{\bf q}_0) \;, \]
    !> with the phonon symmetry matrix \({\bf \Gamma}(\mathcal{S};{\bf q}_0)\) (see [[ph_util_symmetry_G(subroutine)]]).
    subroutine ph_rhopot_rotate_q_canonical( vql0, vql, isym, dfun_mt, dfun_ir )
      use constants, only: zzero, zone, twopi
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use mod_symmetry, only: lsplsymc, symlat, symlatc, vtlsymc, ieqatom, symapp_zfig
      use mod_muffin_tin, only: nrmt
      use phonons_util, only: ph_util_symmetry_G
      use modinput
      !> phonon wavevector \({\bf q}_0\) in lattice coordinates
      real(dp), intent(in) :: vql0(3)
      !> phonon wavevector \({\bf q}\) in lattice coordinates
      real(dp), intent(in) :: vql(3)
      !> global symmetry index of \(\mathcal{S}\)
      integer, intent(in) :: isym
      !> on input: muffin-tin response at \({\bf q}_0\); 
      !> on output: muffin-tin response at \({\bf q}\)
      complex(dp), intent(inout) :: dfun_mt(:,:,:,:,:)
      !> on input: interstitial response at \({\bf q}_0\); 
      !> on output: interstitial response at \({\bf q}\)
      complex(dp), intent(inout) :: dfun_ir(:,:,:)

      integer :: is, ia, ias, js, ja, jas, ka, kas, ip, lspl, shift(3), lda
      real(dp) :: vql_rot(3), a(3), b(3), phi

      complex(dp), allocatable :: dfun_mt_rot(:,:,:,:,:), dfun_ir_rot(:,:,:), sym_G(:,:)

      ! allocate local arrays
      allocate( dfun_mt_rot, source=dfun_mt )
      allocate( dfun_ir_rot, source=dfun_ir )

      ! get shift vector G
      lspl = lsplsymc(isym)
      call r3mtv( dble( symlat(:, :, lspl) ), vql, vql_rot )
      shift = nint( vql_rot - vql0 )
      call assert( sum( abs( vql_rot - vql0 - shift ) ) < 1e-6_dp, &
        '`vql` and `vql0` are not related by given symmetry.' )

      ! get symmetry matrix Gamma
      allocate( sym_G(3*natmtot, 3*natmtot) )
      sym_G = ph_util_symmetry_G( isym, vql0 )

      ! apply symmetry to potential response
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          do ip = 1, 3

            ! muffin-tin part
            do js = 1, nspecies
              do ja = 1, natoms(js)
                jas = idxas(ja, js)
                ka = ieqatom(ja, js, isym)
                kas = idxas(ka, js)
                a = input%structure%speciesarray(js)%species%atomarray(ka)%atom%coord + vtlsymc(:, isym)
                call r3mv( dble( symlat(:, :, lspl) ), a, b )
                b = input%structure%speciesarray(js)%species%atomarray(ja)%atom%coord - b
                phi = twopi * dot_product( vql, b )
                call rotzflm( symlatc(:, :, lspl), dfpt_lmaxvr, nrmt(is), dfpt_lmmaxvr, &
                       dfun_mt(:, :, kas, ip, ias), dfun_mt_rot(:, :, jas, ip, ias) )
                dfun_mt_rot(:, :, jas, ip, ias) = dfun_mt_rot(:, :, jas, ip, ias) * cmplx( cos( phi ), sin( phi ), dp )
              end do
            end do
            ! interstitial region
            dfun_ir_rot(:, ip, ias) = zzero
            call zfftifc( 3, dfpt_Gset%ngrid, -1, dfun_ir(:, ip, ias) )
            call symapp_zfig( symlat(:, :, lspl), vtlsymc(:, isym), vql0, &
                   dfun_ir(:, ip, ias), dfpt_Gset%ngrtot, dfpt_Gset%ivg, dfpt_Gset%igfft, .true., &
                   dfun_ir_rot(:, ip, ias), dfpt_Gset%intgv, dfpt_Gset%ivgig, dfpt_Gset%igfft, .true. )
            call zfftifc( 3, dfpt_Gset%ngrid, 1, dfun_ir_rot(:, ip, ias) )

          end do
        end do
      end do

      lda = product( shape ( dfun_mt(:, :, :, 1, 1) ) )
      call zgemm( 'n', 'c', lda, 3*natmtot, 3*natmtot, zone, &
             dfun_mt_rot, lda, &
             sym_G, 3*natmtot, zzero, &
             dfun_mt, lda )
      lda = product( shape( dfun_ir(:, 1, 1) ) )
      call zgemm( 'n', 'c', lda, 3*natmtot, 3*natmtot, zone, &
             dfun_ir_rot, lda, &
             sym_G, 3*natmtot, zzero, &
             dfun_ir, lda )

      deallocate( dfun_mt_rot, dfun_ir_rot, sym_G )
    end subroutine ph_rhopot_rotate_q_canonical

end module phonons_density_potential
