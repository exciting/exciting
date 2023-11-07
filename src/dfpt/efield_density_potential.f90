!> This module contains procedures for the calculation of the response
!> of the electronic density and the effective potential upon a perturbing
!> electric field \(\delta_{\mathcal{E}_\alpha}\).
module efield_density_potential
  use dfpt_variables
  use dfpt_density_potential
  use efield_variables

  use precision, only: dp
  use asserts, only: assert

  implicit none
  private

  public :: ef_rhopot_gen_dpot
  public :: ef_rhopot_symmetrize

  contains

    !> This subroutine calculates the effective potential response from a given
    !> density response using Weinert's method. See also [[dfpt_rhopot_gen_dpot(subroutine)]].
    !>
    !> If `coulomb=.true.`, then the Coulomb potential response is included.
    !> If `xc=.true.`, then the exchange correlation potential response is included.
    subroutine ef_rhopot_gen_dpot( drho_mt, drho_ir, dpot_mt, dpot_ir, &
        coulomb, xc )
      !> muffin-tin density response as complex spherical harmonics expansion
      complex(dp), intent(in) :: drho_mt(:,:,:)
      !> interstitial density response on real space FFT grid
      complex(dp), intent(in) :: drho_ir(:)
      !> muffin-tin effective potential response as complex spherical harmonics expansion
      complex(dp), intent(out) :: dpot_mt(:,:,:)
      !> interstitial effective potential response on real space FFT grid
      complex(dp), intent(out) :: dpot_ir(:)
      !> include the Coulomb / exchange correlation potential response (default: `.true.`)
      logical, optional, intent(in) :: coulomb, xc

      logical :: vcoul, vxc

      vcoul = .true.
      if( present( coulomb ) ) vcoul = coulomb
      vxc = .true.
      if( present( xc ) ) vxc = xc

      call dfpt_rhopot_gen_dpot( drho_mt, drho_ir, dpot_mt, dpot_ir, &
             dfpt_Gset, ef_jlgr, ef_ylmg, ef_sfacg, &
             coulomb=vcoul, xc=vxc )

    end subroutine ef_rhopot_gen_dpot

    !> This subroutine symmetrizes an electric field perturbed unit cell function.
    !>
    !> The symmetrized function is given by
    !> \[ \delta_{\mathcal{E}_\alpha} f({\bf r})
    !>    = \frac{1}{N_S} \sum\limits_\mathcal{S}^{\mathcal{G}} \sum_{\beta=1}^3
    !>      S_{\alpha\beta} \,
    !>      \delta_{\mathcal{E}_\beta} \bar{f}(\mathcal{S}^{-1}{\bf r}) \;, \]
    !> where \(\delta_{\mathcal{E}_\alpha} \bar{f}\) is the unsymmetrized function.
    subroutine ef_rhopot_symmetrize( zfun_mt, zfun_ir, nsym, isym )
      use constants, only: zzero, twopi
      use mod_symmetry, only: lsplsymc, symlat, symlatc, ieqatom, vtlsymc, symapp_zfig
      use mod_atoms, only: natmtot, natmmax, nspecies, natoms, idxas
      use mod_muffin_tin, only: nrmtmax, nrmt
      use modinput
      !> complex muffin-tin function for each irrep member as complex spherical harmonics expansion
      complex(dp), intent(inout) :: zfun_mt(dfpt_lmmaxvr, nrmtmax, natmtot, *)
      !> complex interstitial function for each irrep member on real space FFT grid
      complex(dp), intent(inout) :: zfun_ir(dfpt_Gset%ngrtot, *)
      !> number of symmetry operations
      integer, intent(in) :: nsym
      !> indices of symmetries in global arrays
      integer, intent(in) :: isym(:)

      integer :: is, ia, ias, ja, ip, ipp, s, lspl
      real(dp) :: sc(3,3)

      complex(dp), allocatable :: zfun_mt1(:,:,:,:), zfun_mt2(:,:)
      complex(dp), allocatable :: zfun_ig1(:,:), zfun_ig2(:)

      ! **** muffin-tin part
      allocate( zfun_mt1(dfpt_lmmaxvr, nrmtmax, natmmax, 3) )
      allocate( zfun_mt2(dfpt_lmmaxvr, nrmtmax) )
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! make a copy of the input function
          do ip = 1, 3
            zfun_mt1(:, :, ia, ip) = zfun_mt(:, : ,ias, ip)
            zfun_mt(:, :, ias, ip) = zzero
          end do
        end do
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          ! loop over symmetries
          do s = 1, nsym
            lspl = lsplsymc(s)                ! index of spatial rotation     
            sc = symlatc(:, :, lspl)          ! Cartesian spatial rotation
            ja = ieqatom(ia, is, s)           ! equivalent atom

            do ip = 1, 3
              ! apply rotation to muffin-tin functions
              call rotzflm( sc, dfpt_lmaxvr, nrmt(is), dfpt_lmmaxvr, zfun_mt1(:, :, ja, ip), zfun_mt2 )
              ! add contribution to input array
              do ipp = 1, 3
                zfun_mt(:, :, ias, ipp) = zfun_mt(:, :, ias, ipp) + sc(ipp, ip) * zfun_mt2
              end do
            end do
          end do
        end do
      end do
      ! normalize
      do ip = 1, 3
        zfun_mt(:, :, :, ip) = zfun_mt(:, :, :, ip) / nsym
      end do
      deallocate( zfun_mt1, zfun_mt2 )

      ! **** interstitial part
      ! transform function to reciprocal space
      allocate( zfun_ig1(dfpt_Gset%ngrtot, 3), source=zzero )
      allocate( zfun_ig2(dfpt_Gset%ngrtot) )
      do ip = 1, 3
        call zfftifc( 3, dfpt_Gset%ngrid, -1, zfun_ir(:, ip) )
        call dfpt_Gset%igfft2ig( zfun_ir(:, ip), zfun_ig1(:, ip) )
        zfun_ir(:, ip) = zzero
      end do
      do s = 1, nsym
        lspl = lsplsymc(s)
        sc = symlatc(:, :, lspl)
        do ip = 1, 3
          ! apply symmetry to interstitial function
          zfun_ig2 = zzero
          call symapp_zfig( symlat(:, :, lspl), vtlsymc(:, s), [0.0_dp, 0.0_dp, 0.0_dp], &
                 zfun_ig1(:, ip), dfpt_Gset%ngvec, dfpt_Gset%ivg, dfpt_Gset%igfft, .false., &
                 zfun_ig2, dfpt_Gset%intgv, dfpt_Gset%ivgig, dfpt_Gset%igfft, .true. )
          ! add contribution to input array
          do ipp = 1, 3
            zfun_ir(:, ipp) = zfun_ir(:, ipp) + sc(ipp, ip) * zfun_ig2
          end do
        end do
      end do
      do ip = 1, 3
        ! transform function to real space
        call zfftifc( 3, dfpt_Gset%ngrid, 1, zfun_ir(:, ip) )
        ! normalize
        zfun_ir(:, ip) = zfun_ir(:, ip) / nsym
      end do
      deallocate( zfun_ig1, zfun_ig2 )
    end subroutine ef_rhopot_symmetrize
end module efield_density_potential
