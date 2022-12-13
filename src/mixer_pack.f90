!> Routines for packing / unpacking functions for mixing
module mixer_pack
  use precision, only: dp

  implicit none
  private

  interface pack_fun
    module procedure :: pack_rfun, pack_zfun
  end interface

  interface unpack_fun
    module procedure :: unpack_rfun, unpack_zfun
  end interface

  public :: pack_fun, unpack_fun

  contains

    !> Pack a real function given in the muffin-tin spheres and the
    !> interstitial region to a real 1D array.
    subroutine pack_rfun( fmt, lmax, rstep, fir, n, f1d)
      use mod_atoms, only: nspecies, natoms, idxas
      use mod_muffin_tin, only: nrmt
      !> muffin-tin function
      real(dp), intent(in) :: fmt(:,:,:)
      !> maximum l for muffin-tin expansion
      integer, intent(in) :: lmax
      !> radial step for muffin-tin expansion
      integer, intent(in) :: rstep
      !> interstitial function
      real(dp), intent(in) :: fir(:)
      !> number of points in 1D array
      integer, intent(inout) :: n
      !> function packed to real 1D array
      real(dp), intent(out) :: f1d(*)

      integer :: is, ia, ias, ir, lm, lmmax

      lmmax = (lmax + 1)**2

      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do ir = 1, nrmt(is), rstep
            do lm = 1, lmmax
              n = n + 1
              f1d(n) = fmt(lm,ir,ias)
            end do
          end do
        end do
      end do
      do ir = 1, size(fir)
        n = n + 1
        f1d(n) = fir(ir)
      end do
    end subroutine pack_rfun

    !> Pack a complex function given in the muffin-tin spheres and the
    !> interstitial region to a real 1D array.
    subroutine pack_zfun( fmt, lmax, rstep, fir, n, f1d)
      use mod_atoms, only: nspecies, natoms, idxas
      use mod_muffin_tin, only: nrmt
      !> muffin-tin function
      complex(dp), intent(in) :: fmt(:,:,:)
      !> maximum l for muffin-tin expansion
      integer, intent(in) :: lmax
      !> radial step for muffin-tin expansion
      integer, intent(in) :: rstep
      !> interstitial function
      complex(dp), intent(in) :: fir(:)
      !> number of points in 1D array
      integer, intent(inout) :: n
      !> function packed to real 1D array
      real(dp), intent(out) :: f1d(*)

      integer :: is, ia, ias, ir, lm, lmmax

      lmmax = (lmax + 1)**2

      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do ir = 1, nrmt(is), rstep
            do lm = 1, lmmax
              n = n + 2
              f1d(n-1) = dble(fmt(lm,ir,ias))
              f1d(n)  = aimag(fmt(lm,ir,ias))
            end do
          end do
        end do
      end do
      do ir = 1, size(fir)
        n = n + 2
        f1d(n-1) = dble(fir(ir))
        f1d(n)  = aimag(fir(ir))
      end do
    end subroutine pack_zfun

    !> Unpack a real function given in the muffin-tin spheres and the
    !> interstitial region from a real 1D array.
    subroutine unpack_rfun( fmt, lmax, rstep, fir, n, f1d)
      use mod_atoms, only: nspecies, natoms, idxas
      use mod_muffin_tin, only: nrmt
      !> muffin-tin function
      real(dp), intent(out) :: fmt(:,:,:)
      !> maximum l for muffin-tin expansion
      integer, intent(in) :: lmax
      !> radial step for muffin-tin expansion
      integer, intent(in) :: rstep
      !> interstitial function
      real(dp), intent(out) :: fir(:)
      !> number of points in 1D array
      integer, intent(inout) :: n
      !> function packed to real 1D array
      real(dp), intent(in) :: f1d(*)

      integer :: is, ia, ias, ir, lm, lmmax

      lmmax = (lmax + 1)**2

      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do ir = 1, nrmt(is), rstep
            do lm = 1, lmmax
              n = n + 1
              fmt(lm,ir,ias) = f1d(n)
            end do
          end do
        end do
      end do
      do ir = 1, size(fir)
        n = n + 1
        fir(ir) = f1d(n)
      end do
    end subroutine unpack_rfun

    !> Unpack a complex function given in the muffin-tin spheres and the
    !> interstitial region from a real 1D array.
    subroutine unpack_zfun( fmt, lmax, rstep, fir, n, f1d)
      use mod_atoms, only: nspecies, natoms, idxas
      use mod_muffin_tin, only: nrmt
      !> muffin-tin function
      complex(dp), intent(out) :: fmt(:,:,:)
      !> maximum l for muffin-tin expansion
      integer, intent(in) :: lmax
      !> radial step for muffin-tin expansion
      integer, intent(in) :: rstep
      !> interstitial function
      complex(dp), intent(out) :: fir(:)
      !> number of points in 1D array
      integer, intent(inout) :: n
      !> function packed to real 1D array
      real(dp), intent(in) :: f1d(*)

      integer :: is, ia, ias, ir, lm, lmmax

      lmmax = (lmax + 1)**2

      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do ir = 1, nrmt(is), rstep
            do lm = 1, lmmax
              n = n + 2
              fmt(lm,ir,ias) = cmplx(f1d(n-1),f1d(n),dp)
            end do
          end do
        end do
      end do
      do ir = 1, size(fir)
        n = n + 2
        fir(ir) = cmplx(f1d(n-1),f1d(n),dp)
      end do
    end subroutine unpack_zfun

end module mixer_pack
