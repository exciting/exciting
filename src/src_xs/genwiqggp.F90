! Copyright(C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genwiqggp
! !INTERFACE:
!
subroutine genwiqggp(flag, iq, igq1, igq2, clwt)
! !USES:
  use modinput
  use modmpi, only: terminate
  use mod_constants, only: pi, twopi, fourpi
  use mod_lattice, only: omega, bvec, binv
  use mod_qpoint, only: nqpt, ngridq
  use modxs, only: vgqc, sptclg, gqc
  use m_genfilname
  use m_getunit
! !DESCRIPTION:
!   Effective integrals of Coulomb interaction. See routine {\tt genwiq2}.
!
! !REVISION HISTORY:
!   Created February 2008 (Sagmeister)
!   Added comments, formatting and truncated potentials (Aurich)
!EOP
!BOC

  implicit none

  ! Arguments
  integer, intent(in) :: flag, iq, igq1, igq2
  real(8), intent(out) :: clwt

  ! Local variables
  integer, parameter :: ns0 = 10, nss = 20
  integer :: np, ns, i1, i2, i3, i, ip, nrbox
  real(8) :: d(3), dv, sum2, t1, t11, t22, t33, t2
  real(8) :: blim(2), blen, vllim(3), ran(3), ranl(3), omegabox
  real(8) :: qsz
  real(8), allocatable :: xa(:), ya(:), c(:)
  real(8) :: vsc(3), v01(3), v02(3), v11(3), v21(3), v31(3), v32(3)
  real(8) :: smod1, smod2, gpq1, gpq2, kxy1, kz1, kxy2, kz2, rccut
  real(8) :: avecz(3)

  ! External functions
  real(8), external :: polynom

  ! Determine G+q-vectors
  v01(:) = vgqc(:, igq1, iq)
  v02(:) = vgqc(:, igq2, iq)

  ! Integrate out singularities and/or improve accuracy in summations
  select case(flag)

    case(0)
      ! Modified coulomb potential in reciprocal space
      clwt = sptclg(igq1, iq) * sptclg(igq2, iq)

    case(1)
      ! Radius of sphere with same volume than subcell
      qsz = (6*pi**2/(omega*nqpt)) ** (1.d0/3.d0)
      if((igq1 .ne. 1) .or. (igq2 .ne. 1)) then
        ! Integrate out 1/q singularity by spherical volume
        clwt = (qsz**2*omega*nqpt/pi) / gqc(igq2, iq)
      else if((igq1 .eq. 1) .and. (igq2 .eq. 1)) then
        ! Integrate out 1/q^2 singularity by spherical volume
        clwt = 2 * qsz * omega * nqpt / pi
      else
        write(*,*)
        write(*, '("Error(genwiqggp): Analytic method chosen for regular case")')
        write(*,*)
        call terminate
      end if

    case(2)
      np = 2
      ! Higher order extrapolation for 1/q^2 term
      if((igq1 .eq. 1) .and. (igq2 .eq. 1)) np = 3
      allocate(xa(np), ya(np), c(np))
      ! Loop over different subdivisions
      ns = ns0
      do ip = 1, np
        ! Subdivision vectors in lattice coordinates
        do i = 1, 3
          d(i) = 1.d0 / (dble(ngridq(i)*2*ns))
        end do
        ! Smallest volume element (we drop the (2pi)^3/omega factor!)
        dv = d(1) * d(2) * d(3)
        ! Compute the integral of 1/(|p+q+g||p+q+gp|) over the small
        ! fraction of the Brillouin zone centered at the gamma-point
        sum2 = 0.d0
        ! The p-point grid is started here
        do i1 = - ns, ns - 1
          t11 = dble(i1) * d(1)
          do i2 = - ns, ns - 1
            t22 = dble(i2) * d(2)
            do i3 = - ns, ns - 1
              t33 = dble(i3) * d(3)
              ! P-vector
              vsc(:) = t11 * bvec(:, 1) + t22 * bvec(:, 2) + t33 * bvec(:, 3)
              ! P+q+g vector
              v31(:) = v01(:) + vsc(:)
              ! P+q+gp vector
              v32(:) = v02(:) + vsc(:)
              t2 = sqrt(sum(v31**2)*sum(v32**2))
              ! Check if integrand would approach singularity
              if(t2 .gt. 1.d-14) then
                sum2 = sum2 + 1.d0 / t2
              end if
            end do
          end do
        end do
        sum2 = sum2 * dv
        xa(ip) = dv ** (1.d0/3.d0)
        ya(ip) = sum2
        ! Increment number of subdivisions
        ns = ns + nss
      end do
      ! Extrapolate the volume element to zero with a polynomial
      clwt = polynom(0, np, xa, ya, c, 0.d0) * fourpi * nqpt
      deallocate(xa, ya, c)

    case(3)
      ! Rim method: Documentation of self code by Andrea Marini
      ! Find maximum extension of small Brillouin zone
      blim(:) = 0.d0
      do i1 = - 1, 1, 2
        v11(:) = dble(i1) * bvec(:, 1) / dble(ngridq(1)*2)
        do i2 = - 1, 1, 2
          v21(:) = v11(:) + dble(i2) * bvec(:, 2) / dble(ngridq(2)*2)
          do i3 = - 1, 1, 2
            v31(:) = v21(:) + dble(i3) * bvec(:, 3) / dble(ngridq(3)*2)
            do i = 1, 3
              ! Lower limit for box
              if(v31(i) .lt. blim(1)) blim(1) = v31(i)
              ! Upper limit for box
              if(v31(i) .gt. blim(2)) blim(2) = v31(i)
            end do
          end do
        end do
      end do
      ! Box length
      blen = blim(2) - blim(1)
      ! Limits of sbz in lattice
      vllim(:) = 1.d0 / dble(2*ngridq(:))
      ! Needs high value(above 1e6) to converge - inferior to method nr. 2
      nrbox = 1000000
      t1 = 0.d0
      omegabox = blen ** 3
      do i = 1, 10000
        call random_number(t2)
      end do
      do i = 1, nrbox
        call random_number(ran)
        ! Map random number to box
        ran(:) = ran(:) * blen
        ran(:) = ran(:) + blim(1)
        ! Check if random vector is in sbz
        ranl = matmul(binv, ran)
        if(all(ranl .gt.-vllim) .and. (all(ranl .lt. vllim))&
          & .and. (sum(abs(ran)) .gt. 1.d-14)) then
          t2 = sum((v01+ran)**2) * sum((v02+ran)**2)
          t1 = t1 + 1.d0 / sqrt(t2)
        end if
      end do
      clwt = t1 * (omegabox/nrbox) * fourpi * nqpt * omega / (twopi) ** 3

    ! 0d cutoff for Coulomb potential
    case(4)

      ! Get a_z
      avecz(:) = input%structure%crystal%basevect(:,3)
      ! Get cut parameter 
      rccut = 0.5d0*dsqrt(dot_product(avecz,avecz))
      ! Get length of G+q and G'+q
      gpq1 = gqc(igq1, iq)
      gpq2 = gqc(igq2, iq)
      ! Make modifier for square root of Coulomb potential
      smod1 = dsqrt(1.d0-dcos(dsqrt(gpq1)*rccut))
      smod2 = dsqrt(1.d0-dcos(dsqrt(gpq2)*rccut))
      ! GG' symmetrized and truncated Coulomb potential 
      clwt = sptclg(igq1, iq)*smod1 * sptclg(igq2, iq)*smod2

    ! 2d cutoff for Coulomb potential
    case(5)

      ! Get a_z
      avecz(:) = input%structure%crystal%basevect(:,3)
      ! Get cut parameter 
      rccut = 0.5d0*dsqrt(dot_product(avecz,avecz))
      ! Modified coulomb potential in reciprocal space
      kxy1 = dsqrt(v01(1)*v01(1)+v01(2)*v01(2))
      kz1 = dabs(v01(3))
      kxy2 = dsqrt(v02(1)*v02(1)+v02(2)*v02(2))
      kz2 = dabs(v02(3))
      smod1 = dsqrt(1.d0-dexp(-kxy1*rccut)*dcos(kz1*rccut))
      smod2 = dsqrt(1.d0-dexp(-kxy2*rccut)*dcos(kz2*rccut))
      clwt = sptclg(igq1, iq)*smod1 * sptclg(igq2, iq)*smod2

    case default
      write(*,*) "Error(genwiqggq): Invalid flag passed."
      call terminate

  end select
end subroutine genwiqggp
!EOC
