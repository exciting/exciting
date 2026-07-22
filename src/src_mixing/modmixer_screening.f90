!> Screened Coulomb helpers used by Kerker and density-Pulay mixing.
Module modmixer_screening
  Use constants, Only: fourpi
  Use m_zfftifc, Only: zfftifc
  Use modinput, Only: input
  Use modmain, Only: gc, igfft, idxas, lmmaxvr, natoms, natmtot, ngvec, ngrid, ngrtot, &
 &    nrmt, nrmtmax, nspecies, omega, rmt, sfacg, spnrmax, spr, y00, ylmg, zil
  Use modmixer_special_functions, Only: msbesseli, msbesselk
  Use precision, Only: dp
  implicit none

Contains

  !> Compute the screened correction for muffin-tin and interstitial residuals.
  Subroutine potcorr(residualmt, residualir, lambda, svclmt, svclir)
    implicit none
    !> Residual inside muffin tins in real spherical-harmonic representation.
    real(dp), intent(in) :: residualmt(lmmaxvr, nrmtmax, natmtot)
    !> Residual in the interstitial region on the FFT grid.
    real(dp), intent(in) :: residualir(ngrtot)
    !> Screening parameter \( \lambda \).
    real(dp), intent(in) :: lambda
    !> Screened muffin-tin correction in real spherical-harmonic representation.
    real(dp), intent(out) :: svclmt(lmmaxvr, nrmtmax, natmtot)
    !> Screened interstitial correction on the FFT grid.
    real(dp), intent(out) :: svclir(ngrtot)
    complex(dp), allocatable :: szvclir(:), zrhoir(:), szvclmt(:, :, :), zrhomt(:, :, :)
    complex(dp) :: szrho0
    integer :: ia, ias, ir, is

    allocate(zrhomt(lmmaxvr, nrmtmax, natmtot))
    allocate(zrhoir(ngrtot))
    allocate(szvclmt(lmmaxvr, nrmtmax, natmtot))
    allocate(szvclir(ngrtot))

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        do ir = 1, nrmt(is)
          call rtozflm(input%groundstate%lmaxvr, residualmt(:, ir, ias), zrhomt(:, ir, ias))
        end do
      end do
    end do
    zrhoir = residualir

    call szpotcorr(nrmt, nrmtmax, spnrmax, spr, 1, gc, lambda, ylmg, sfacg, zrhomt, zrhoir, &
 &    szvclmt, szvclir, szrho0)

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        do ir = 1, nrmt(is)
          call ztorflm(input%groundstate%lmaxvr, szvclmt(:, ir, ias), svclmt(:, ir, ias))
        end do
      end do
    end do
    svclir = real(szvclir, dp)

    deallocate(zrhomt, zrhoir, szvclmt, szvclir)
  End Subroutine

  !> Solve the screened Poisson problem inside a single muffin tin.
  Subroutine szpotclmt(lmax, nr, r, lambda, ld, zrhomt, szvclmt)
    implicit none
    !> Maximum angular momentum of the spherical-harmonic expansion.
    integer, intent(in) :: lmax
    !> Number of radial mesh points.
    integer, intent(in) :: nr
    !> Leading dimension of the spherical-harmonic coefficient array.
    integer, intent(in) :: ld
    !> Radial mesh points for the muffin tin.
    real(dp), intent(in) :: r(nr)
    !> Screening parameter \( \lambda \).
    real(dp), intent(in) :: lambda
    !> Complex muffin-tin source term in spherical-harmonic representation.
    complex(dp), intent(in) :: zrhomt(ld, nr)
    !> Complex screened potential inside the muffin tin.
    complex(dp), intent(out) :: szvclmt(ld, nr)
    integer :: ir, l, lm, m
    real(dp) :: cf(3, nr), fr1(nr), fr2(nr), fr3(nr), fr4(nr), gr1(nr), gr2(nr), gr3(nr), gr4(nr)
    real(dp) :: il(0:lmax, nr), kl(0:lmax, nr), t1, t2, t3, x, z1, z2, z3, z4, z5, z6

    do ir = 1, nr
      do l = 0, lmax
        x = r(ir) * lambda
        call msbesseli(l, x, il(:, ir))
        call msbesselk(l, x, kl(:, ir))
      end do
    end do

    lm = 0
    do l = 0, lmax
      t1 = fourpi * lambda
      do m = -l, l
        lm = lm + 1
        do ir = 1, nr
          t2 = il(l, ir) * r(ir)**2
          t3 = kl(l, ir) * r(ir)**2
          fr1(ir) = t2 * real(zrhomt(lm, ir), dp)
          fr2(ir) = t2 * aimag(zrhomt(lm, ir))
          fr3(ir) = t3 * real(zrhomt(lm, ir), dp)
          fr4(ir) = t3 * aimag(zrhomt(lm, ir))
        end do
        call fderiv(-1, nr, r, fr1, gr1, cf)
        call fderiv(-1, nr, r, fr2, gr2, cf)
        call fderiv(-1, nr, r, fr3, gr3, cf)
        call fderiv(-1, nr, r, fr4, gr4, cf)
        z1 = gr3(nr)
        z2 = gr4(nr)
        do ir = 1, nr
          z3 = kl(l, ir)
          z4 = il(l, ir)
          z5 = z3 * gr1(ir) + z4 * (z1 - gr3(ir))
          z6 = z3 * gr2(ir) + z4 * (z2 - gr4(ir))
          szvclmt(lm, ir) = t1 * cmplx(z5, z6, kind=dp)
        end do
      end do
    end do
  End Subroutine

  !> Assemble screened Coulomb corrections in reciprocal and muffin-tin space.
  Subroutine szpotcorr(nr, nrmax, ld, r, igp0, gpc, lambda, ylmgp, sfacgp, zrhomt, zrhoir, &
 &    szvclmt, szvclir, szrho0)
    implicit none
    !> Number of radial mesh points for each species.
    integer, intent(in) :: nr(nspecies)
    !> Maximum radial mesh size.
    integer, intent(in) :: nrmax
    !> Leading dimension of the radial mesh array.
    integer, intent(in) :: ld
    !> Index of the \( G=0 \) vector in the reciprocal grid list.
    integer, intent(in) :: igp0
    !> Species-dependent radial meshes.
    real(dp), intent(in) :: r(ld, nspecies)
    !> Reciprocal-grid magnitudes.
    real(dp), intent(in) :: gpc(ngvec)
    !> Screening parameter \( \lambda \).
    real(dp), intent(in) :: lambda
    !> Spherical harmonics on reciprocal grid directions.
    complex(dp), intent(in) :: ylmgp(lmmaxvr, ngvec)
    !> Structure factors on the reciprocal grid.
    complex(dp), intent(in) :: sfacgp(ngvec, natmtot)
    !> Complex charge residual in muffin-tin representation.
    complex(dp), intent(in) :: zrhomt(lmmaxvr, nrmax, natmtot)
    !> Complex charge residual in interstitial representation.
    complex(dp), intent(in) :: zrhoir(ngrtot)
    !> Complex screened potential in muffin-tin representation.
    complex(dp), intent(out) :: szvclmt(lmmaxvr, nrmax, natmtot)
    !> Complex screened potential in interstitial representation.
    complex(dp), intent(out) :: szvclir(ngrtot)
    !> \( G=0 \) component of the screened interstitial potential before it is removed.
    complex(dp), intent(out) :: szrho0
    integer :: ia, ifg, ig, ias, ir, is, l, lm, local_npsd, m
    real(dp), allocatable :: ilmt(:, :), jlgpr(:, :, :), klmt(:, :)
    real(dp) :: fpo, ill(0:input%groundstate%lmaxvr, nrmax, nspecies)
    real(dp) :: rl(0:input%groundstate%lmaxvr, nrmax), t1, t2, x1, x2, z1
    complex(dp) :: qi(lmmaxvr, natmtot), qmt(lmmaxvr, natmtot), vilm(lmmaxvr), zrp(lmmaxvr), zsum1, zsum2, zt1, zt2
    real(dp) :: factnm
    external factnm

    local_npsd = nint(0.25d0 * input%groundstate%gmaxvr * maxval(rmt(1:nspecies)))
    fpo = fourpi / omega
    allocate(ilmt(0:input%groundstate%lmaxvr + local_npsd + 1, nspecies))
    allocate(jlgpr(0:input%groundstate%lmaxvr + local_npsd + 1, ngvec, nspecies))
    allocate(klmt(0:input%groundstate%lmaxvr + local_npsd + 1, nspecies))

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        call szpotclmt(input%groundstate%lmaxvr, nr(is), r(:, is), lambda, lmmaxvr, zrhomt(:, :, ias), szvclmt(:, :, ias))
      end do
    end do

    do is = 1, nspecies
      do ig = 1, ngvec
        t1 = gpc(ig) * rmt(is)
        do l = 0, input%groundstate%lmaxvr + local_npsd + 1
          call sbessel(l, t1, jlgpr(:, ig, is))
        end do
      end do
    end do

    do is = 1, nspecies
      x1 = rmt(is) * lambda
      do l = 0, input%groundstate%lmaxvr + local_npsd + 1
        call msbesselk(l, x1, klmt(:, is))
        call msbesseli(l, x1, ilmt(:, is))
      end do
    end do

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        lm = 0
        do l = 0, input%groundstate%lmaxvr
          t1 = factnm(2 * l + 1, 2) / (lambda**l)
          t2 = 1d0 / (klmt(l, is) * fourpi * lambda)
          do m = -l, l
            lm = lm + 1
            qmt(lm, ias) = t1 * t2 * szvclmt(lm, nr(is), ias)
          end do
        end do
      end do
    end do

    szvclir = zrhoir
    call zfftifc(3, ngrid, -1, szvclir)
    qi = 0.0d0
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        do ig = 1, ngvec
          ifg = igfft(ig)
          if (gpc(ig) .gt. input%structure%epslat) then
            zt1 = szvclir(ifg) * sfacgp(ig, ias) / (gpc(ig)**2 + lambda**2)
            lm = 0
            do l = 0, input%groundstate%lmaxvr
              if (l .eq. 0) then
                zt2 = zt1 * fourpi * rmt(is)**2 * ((lambda * jlgpr(0, ig, is) * &
 &                (ilmt(1, is) + ilmt(0, is) / (lambda * rmt(is)))) - &
 &                (gpc(ig) * ((jlgpr(0, ig, is) / (gpc(ig) * rmt(is))) - jlgpr(1, ig, is)) * ilmt(0, is)))
              else
                zt2 = zt1 * fourpi * zil(l) * rmt(is)**2 * factnm(2 * l + 1, 2) / (lambda**l) * &
 &                ((lambda * jlgpr(l, ig, is) * ilmt(l - 1, is)) - (gpc(ig) * jlgpr(l - 1, ig, is) * ilmt(l, is)))
              end if
              do m = -l, l
                lm = lm + 1
                qi(lm, ias) = qi(lm, ias) + zt2 * conjg(ylmgp(lm, ig))
              end do
            end do
          else
            t1 = fourpi * y00 * rmt(is)**2 * ilmt(1, is) / lambda
            qi(1, ias) = qi(1, ias) + t1 * szvclir(ifg)
          end if
        end do
      end do
    end do

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        lm = 0
        do l = 0, input%groundstate%lmaxvr
          t1 = 1d0 / factnm(2 * l + 1, 2)
          do m = -l, l
            lm = lm + 1
            zrp(lm) = (qmt(lm, ias) - qi(lm, ias)) * t1
          end do
        end do
        do ig = 1, ngvec
          ifg = igfft(ig)
          if (gpc(ig) .gt. input%structure%epslat) then
            z1 = gpc(ig)
            zt1 = fpo * conjg(sfacgp(ig, ias)) / (z1**(local_npsd + 1))
            lm = 0
            do l = 0, input%groundstate%lmaxvr
              zsum1 = jlgpr(local_npsd + l + 1, ig, is) * (lambda**(l + local_npsd + 1)) / &
 &              ilmt(local_npsd + l + 1, is) * conjg(zil(l))
              lm = lm + 1
              zsum2 = zrp(lm) * ylmgp(lm, ig)
              do m = -l + 1, l
                lm = lm + 1
                zsum2 = zsum2 + zrp(lm) * ylmgp(lm, ig)
              end do
              szvclir(ifg) = szvclir(ifg) + zt1 * zsum1 * zsum2
            end do
          else
            z1 = (fpo * y00 * rmt(is)**(local_npsd + 1) * lambda**(local_npsd + 1)) / &
 &            ((factnm(2 * local_npsd + 3, 2)**2) * ilmt(local_npsd + 1, is))
            szvclir(ifg) = szvclir(ifg) + z1 * zrp(1)
          end if
        end do
      end do
    end do

    ifg = igfft(igp0)
    szrho0 = szvclir(ifg)
    szvclir(ifg) = 0.d0
    do ig = 1, ngvec
      ifg = igfft(ig)
      szvclir(ifg) = fourpi * szvclir(ifg) / (gpc(ig)**2 + lambda**2)
    end do

    do is = 1, nspecies
      do ir = 1, nr(is)
        x2 = r(ir, is) * lambda
        do l = 0, input%groundstate%lmaxvr
          call msbesseli(l, x2, ill(:, ir, is))
          rl(l, ir) = ill(l, ir, is) / ilmt(l, is)
        end do
      end do
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        vilm = 0.d0
        do ig = 1, ngvec
          ifg = igfft(ig)
          zt1 = fourpi * szvclir(ifg) * sfacgp(ig, ias)
          lm = 0
          do l = 0, input%groundstate%lmaxvr
            zt2 = jlgpr(l, ig, is) * zt1 * zil(l)
            do m = -l, l
              lm = lm + 1
              vilm(lm) = vilm(lm) + zt2 * conjg(ylmgp(lm, ig))
            end do
          end do
        end do
        lm = 0
        do l = 0, input%groundstate%lmaxvr
          do m = -l, l
            lm = lm + 1
            zt1 = vilm(lm) - szvclmt(lm, nr(is), ias)
            do ir = 1, nr(is)
              szvclmt(lm, ir, ias) = szvclmt(lm, ir, ias) + zt1 * rl(l, ir)
            end do
          end do
        end do
      end do
    end do

    call zfftifc(3, ngrid, 1, szvclir)
    deallocate(ilmt, jlgpr, klmt)
  End Subroutine

End Module
