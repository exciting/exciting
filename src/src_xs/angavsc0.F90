! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: angavsc0
! !INTERFACE:
subroutine angavsc0(n, nmax, scrnh, scrnw, scrn, scieff)
! !USES:
  use modmpi
  use mod_constants, only: zzero, pi, twopi, fourpi
  use mod_qpoint, only: ngridq
  use mod_lattice, only: omega, binv
  use modinput, only: input
  use modxs, only: dielten0, dielten, symt2, tleblaik,&
                 & lmmaxdielt, sptclg
  use invert
! !INPUT/OUTPUT PARAMETERS:
!   IN:
!   n, integer               : Number of G+q vectors for current q-point
!   nmax, integer            : Maximum number of G+q vectors over all q-points
!   scrnh(3,3), complex(8)   : Head of dielectric matrix
!   scrnw(n,2,3), complex(8) : Wings of dielectric matrix
!   scrn(n,n), complex(8)    : Body of dielectric matrix
!   OUT:
!   scieff(nmax,nmax), complex(8) : Screened coulomb potential 
!   
! !DESCRIPTION:
! This routine deals with the divergences of 
! $\frac{\tilde{\varepsilon}_{G,G'}(q,\omega=0)}{|G+q||G'+q|}$
! for $G=G'=q=0$ by averaging around the Gamma point.
!
! !REVISION HISTORY:
! Added to documentation scheme. (Aurich)
!
!EOP
!BOC

  implicit none

  ! Arguments
  integer, intent(in) :: n, nmax
  complex(8), intent(in) :: scrn(n, n), scrnw(n, 2, 3), scrnh(3, 3)
  complex(8), intent(out) :: scieff(nmax, nmax)

  ! Local variables
  integer, parameter :: nsphcov = 5810, iq0 = 1
  integer :: iop, jop, j1, j2, itp, lm, ntpsph
  real(8) :: vomega, t00, r, qsz
  complex(8) :: dtns(3, 3), w1, w2,scol(3),zt
  real(8), allocatable :: plat(:, :), p(:), tp(:, :), spc(:,:), w(:)
  complex(8), allocatable :: m00lm(:), mx0lm(:), mxxlm(:)
  complex(8), allocatable :: ei00 (:), eix0 (:), ei0x(:), eixx(:)
  complex(8), allocatable :: eitmp(:), eilmtmp(:)
  complex(8), allocatable :: ei00lm(:), eix0lm(:), ei0xlm(:), eixxlm(:), spc2(:)
  complex(8), allocatable :: ylm(:), zylm(:, :)
  complex(8), allocatable :: b(:, :), bi(:, :), u(:, :), v(:, :), s(:, :), t(:, :)
  complex(8), allocatable :: e3 (:, :), ie3(:, :)
  integer :: i1, i2

  scieff=zzero

  ! Crystal volume
  vomega = omega * product(ngridq)

  ! Wigner-Seitz radius and spherical approximation to 1/q^2 average
  qsz = (6*pi**2/vomega) ** (1.d0/3.d0)

  ! Weight for sqrt(4pi)/q based on Wigner-Seitz radius
  w1 = qsz ** 2 * vomega / (4.d0*pi**2) * sqrt(fourpi)

  ! Weight for 4pi/q^2 based on Wigner-Seitz radius
  w2 = 2 * qsz * vomega / pi

  ! Calculate RPA dielectric tensor including local field effects
  dielten0(:, :) = scrnh(:, :)

  q0: if(n .gt. 1) then

    allocate(b(n-1, n-1), bi(n-1, n-1), u(n-1, 3), v(3, n-1), s(n-1, 3), t(3, n-1))

    ! Body of dielectric matrix
    b(:, :) = scrn(2:, 2:)

    ! Column wing
    u(:, :) = scrnw(2:, 2, :)

    ! Row wing
    v(:, :) = transpose(scrnw(2:, 1, :))

    select case(input%xs%bse%scrherm)
      case(0) ! Default
        ! Use full matrix (both wings and full body)
      case(1)
        ! Hermitian average matrix (average both wings and body)
        b = 0.5d0 * (b+conjg(transpose(b)))
        u = 0.5d0 * (u+conjg(transpose(v)))
        v = conjg(transpose(u))
      case(2)
        ! Use upper triangle (take row wing, assign column wing)
        u = conjg(transpose(v))
      case(3)
        ! Use lower triangle (take column wing, assign row wing)
        v = conjg(transpose(u))
      case default
        write(*,*)
        write(*, '("Error(angavsc0): Not a valid flag:", i6)') input%xs%bse%scrherm
        write(*,*)
        call terminate
    end select

    ! Invert body (optionally including hermitian average)
    call zinvert_hermitian(input%xs%bse%scrherm, b, bi)
    s = matmul(bi, u)
    t = matmul(v, bi)
    dielten = dielten0 - matmul(v, s)

  else q0

     dielten = dielten0

  end if q0

  ! Symmetrize the dielectric tensor
  dtns(:, :) = dielten(:, :)
  do iop = 1, 3
    do jop = 1, 3
      call symt2app(iop, jop, 1, symt2, dtns, dielten(iop, jop))
    end do
  end do

  ! Calculate averaged screened coulomb interaction in Fourier space at gamma point
  select case(trim(input%xs%bse%sciavtype))
    ! Default
    case('spherical')
      ! Scaling factor
      t00 = (omega/(twopi)**3) * product(ngridq)
      ! Number of points on sphere (tleblaik is true by default)
      if(tleblaik) then
        ntpsph = input%xs%bse%nleblaik
      else
        ntpsph = nsphcov
      end if

      ! Sanity check
      if(lmmaxdielt .gt. ntpsph) then
        write(*,*)
        write(*, '("Error(angavdm0): lmmaxdielt.gt.ntpsph: ", 2i6)') lmmaxdielt, ntpsph
        write(*,*)
        stop
      end if

      ! Allocate local arrays
      allocate(plat(3, ntpsph), p(ntpsph))
      allocate(m00lm(lmmaxdielt), mx0lm(lmmaxdielt), mxxlm(lmmaxdielt))
      allocate(ei00(ntpsph), eix0(ntpsph), ei0x(ntpsph), eixx(ntpsph))
      allocate(ei00lm(lmmaxdielt), eix0lm(lmmaxdielt))
      allocate(ei0xlm(lmmaxdielt), eixxlm(lmmaxdielt))
      allocate(ylm(lmmaxdielt), zylm(ntpsph, lmmaxdielt))
      allocate(tp(2, ntpsph), spc(3, ntpsph))
      allocate(spc2(ntpsph))
      allocate(w(ntpsph))

      ! Default
      if(tleblaik) then
        ! Generate lebedev laikov grid
        call leblaik(ntpsph, spc, w)
        ! Generate tetha and phi angles
        do itp = 1, ntpsph
          call sphcrd(spc(:, itp), r, tp(:, itp))
        end do
      else
        ! Distribution is assumed to be uniform
        w(:) = 1.d0 / ntpsph

        ! Generate spherical covering set(angles and coordinates)
        call sphcover(ntpsph, tp)
        spc(1, :) = sin(tp(1, :)) * cos(tp(2, :))
        spc(2, :) = sin(tp(1, :)) * sin(tp(2, :))
        spc(3, :) = cos(tp(1, :))
      end if

      ! Generate spherical harmonics on covering set
      do itp = 1, ntpsph
        call genylm(input%xs%bse%lmaxdielt, tp(:, itp), ylm)
        zylm(itp, :) = ylm(:)
      end do

      ! Unit vectors of spherical covering set in lattice coordinates
      plat = matmul(binv, spc)

      ! Distances to subcell cell boundaries in reciprocal space
      do itp = 1, ntpsph
        p(itp) = 1.d0 / (2.d0*maxval(abs(ngridq(:)*plat(:, itp)), 1))
      end do

      ! Calculate function on covering set
      do itp = 1, ntpsph
        ! Head, 1/(p*l*p)
        ei00(itp) = 1.d0 / dot_product(spc(:, itp), matmul(dielten, spc(:, itp)))
      end do
 
      ! Calculate lm-expansion coefficients
      do lm = 1, lmmaxdielt
        ei00lm(lm) = fourpi * dot_product(zylm(:, lm), ei00*w)
        m00lm(lm) = fourpi * dot_product(zylm(:, lm), p*w)
        mx0lm(lm) = fourpi * dot_product(zylm(:, lm), p**2/2.d0*w)
        mxxlm(lm) = fourpi * dot_product(zylm(:, lm), p**3/3.d0*w)
      end do

      ! Subcell average(head)
      scieff(1, 1) = fourpi * t00 * dot_product(m00lm, ei00lm)

      fullinv: if(input%xs%bse%scrherm .eq. 0) then

        ! Loop over(g,gp) indices
        do j1 = 2, n
          scol(1:3)=s(j1-1,1:3)

          do itp = 1, ntpsph
            spc2(itp) = (spc(1, itp)* scol(1)+spc(2, itp)&
              &* scol(2)+spc(3, itp)* scol(3)) * ei00(itp)
          end do

          do itp = 1, ntpsph
            ! Wing, -p*s/(p*l*p)
            eix0 (itp) = - spc2(itp) *w(itp) !* ei00(itp)
            ! Wing, -p*t/(p*l*p)
            ei0x(itp) = - dot_product(spc(:, itp), t(:, j1-1)) * ei00(itp)*w(itp)
          end do

          call zgemv('c',ntpsph,lmmaxdielt,fourpi,zylm,ntpsph,eix0,1,zzero,eix0lm,1)
          call zgemv('c',ntpsph,lmmaxdielt,fourpi,zylm,ntpsph,ei0x,1,zzero,ei0xlm,1)
          ! Subcell average (wings)
          scieff(j1, 1) = sqrt(fourpi)*sptclg(j1, iq0)*t00*dot_product(mx0lm, eix0lm)
          scieff(1, j1) = sqrt(fourpi)*sptclg(j1, iq0)*t00*dot_product(mx0lm, ei0xlm)
        end do

        ! Average body of screened coulomb interaction 
        avbody: if(input%xs%bse%sciavbd) then

          allocate(eitmp(ntpsph),eilmtmp(lmmaxdielt))

          do i1=1,3
            do i2=1,3
              eitmp(:)=w(:)*spc(i1,:)*spc(i2,:)*ei00(:)
              call zgemv('c',ntpsph,lmmaxdielt,fourpi,zylm,&
                & ntpsph,eitmp,1,zzero,eilmtmp,1)
              zt=t00 * dot_product(mxxlm, eilmtmp)
              do j2=2,n
                do j1=2,n
                  scieff(j1, j2) = scieff(j1, j2)+sptclg(j1, iq0)&
                   &* sptclg(j2, iq0) * zt *s(j1-1,i1)*t(i2,j2-1)
                end do
              end do
            end do
          end do 

          eitmp(:)=cmplx(w(:),0d0,8)
          call zgemv('c',ntpsph,lmmaxdielt,fourpi,zylm,ntpsph,eitmp,1,zzero,eilmtmp,1)
          zt=t00 * dot_product(mxxlm, eilmtmp)
          do j2=2,n
            do j1=2,n
              scieff(j1, j2) = scieff(j1, j2)+ sptclg(j1, iq0)&
                &* sptclg(j2, iq0) * zt *bi(j1-1,j2-1)
            end do
          end do

          deallocate(eitmp,eilmtmp)

        else avbody

          ! No subcell average(body)
          scieff(j1, 2:n) = bi(j1-1, :)

        end if avbody

      else fullinv


        do j1 = 2, n

          scol(1:3)=s(j1-1,1:3)

          do itp = 1, ntpsph
            spc2(itp) = spc(1, itp)* scol(1)+spc(2, itp)* scol(2)+spc(3, itp)* scol(3)
          end do

          do itp = 1, ntpsph
            ! Wing, -p*s/(p*l*p)
            eix0 (itp) = - spc2(itp) * ei00 (itp)
          end do

          call zgemv('c',ntpsph,lmmaxdielt,fourpi,zylm,ntpsph,ei0x,1,zzero,ei0xlm,1)

          ! Subcell average(wings)
          scieff(j1, 1) = sqrt(fourpi)*sptclg(j1, iq0)*t00 * dot_product(mx0lm, eix0lm)
          scieff(1, j1) = conjg(scieff(j1, 1))

          ! Average body of coulomb interaction
          if(input%xs%bse%sciavbd) then
            do j2 = j1, n
              do itp = 1, ntpsph
                ! Body, b^-1 + p*s p*t/(p*l*p)
                eixx(itp) = w(itp)*( bi(j1-1, j2-1)&
                  &+ (spc(1, itp)* scol(1)+spc(2, itp)* scol(2)+spc(3, itp)* scol(3))&
                  &* (spc(1, itp)* t(1,j2-1)+spc(2, itp)* t(2,j2-1)&
                  &+ spc(3, itp)* t(3,j2-1))* ei00(itp) )
              end do
              call zgemv('c',ntpsph,lmmaxdielt,fourpi,zylm,ntpsph,eixx,1,zzero,eixxlm,1)
              ! Subcell average(body)
              scieff(j1, j2) = sptclg(j1, iq0) * sptclg(j2, iq0)&
               &* t00 * dot_product(mxxlm, eixxlm)
              scieff(j2, j1) = conjg(scieff(j1, j2))
            end do
          else
            ! No subcell average(body)
            scieff(j1, 2:n) = bi(j1-1, :)
          end if

        end do

      end if fullinv

      deallocate(ei00, eix0, ei0x, eixx, ei00lm, eix0lm, ei0xlm, m00lm, mx0lm, mxxlm)
      deallocate(ylm, zylm, tp, spc, w, plat, p)
      deallocate(spc2)

    case('screendiag', 'invscreendiag')

      if(input%xs%bse%sciavbd) then
        write(*,*)
        write(*, '("Error(angavsc0): (inv)screendiag-method does not allow&
          & for averaging the body of w")')
        write(*,*)
        stop
      end if
      allocate(e3(n+2, n+2), ie3(n+2, n+2))
      ! Invert dielectric matrix including 3 times g=0 according to the limits
      ! Q->0_x, q->0_y, q->0_z
      ! G=0, g'=0 elements
      e3(1:3, 1:3) = dielten0(:, :)
      ! G!=0, g'=0 components and vice versa
      if(n .gt. 1) then
        do i1 = 1, 3
          do j2 = 2, n
            e3(i1, j2+2) = scrnw(j2, 1, i1)
          end do
        end do
        do j1 = 2, n
          do i2 = 1, 3
            e3(j1+2, i2) = scrnw(j1, 2, i2)
          end do
        end do
        do j1 = 2, n
          do j2 = 2, n
            e3(j1+2, j2+2) = scrn(j1, j2)
          end do
        end do
      end if

      call zinvert_hermitian(input%xs%bse%scrherm, e3, ie3)

      ! Select again
      select case(trim(input%xs%bse%sciavtype))
        case('screendiag')
          ! Head
          scieff(1, 1) = w2 * 1.d0 / ((e3(1, 1)+ie3(2, 2)+ie3(3, 3))/3.d0)
          if(n .gt. 1) then
            ! Wings, set to zero in this approximation
            scieff(1, 2:n) = zzero
            scieff(2:n, 1) = zzero
            ! Body, only diagonal is assigned
            scieff(2:n, 2:n) = zzero
            forall(j1=2:n)
              scieff(j1, j1) = sptclg(j1, iq0) ** 2 / e3 (j1+2, j1+2)
            end forall
          end if

        case('invscreendiag')
          ! Head
          scieff(1, 1) = w2 * (ie3(1, 1)+ie3(2, 2)+ie3(3, 3)) / 3.d0
          ! Wings
          if(n .gt. 1) then
            forall(j1=2:n)
              scieff(j1, 1) = w1 * sptclg(j1, iq0)&
                &* (ie3(j1+2,1)+ie3(j1+2, 2)+ie3(j1+3, 3)) / 3.d0
              scieff(1, j1) = conjg(scieff(j1, 1))
            end forall
            ! Body
            forall(j1=2:n, j2=2:n)
              scieff(j1, j2) = sptclg(j1, iq0) * sptclg(j2, iq0) * ie3 (j1+2, j2+2)
            end forall
          end if
      end select

      deallocate(e3, ie3)

    case('none')
      iop = 1 !!!only x-component here!!!
      ! Longitudinal treatment, three components of vanishing q(direction)
      allocate(e3(n, n), ie3(n, n))
      e3(1, 1) = scrnh(iop, iop)
      e3(1, 2:n) = scrnw(2:, 1, iop)
      e3(2:n, 1) = scrnw(2:, 2, iop)
      e3(2:n, 2:n) = scrn(2:, 2:)
      call zinvert_hermitian(input%xs%bse%scrherm, e3, ie3)
      write(*,*) 'eps^{-1}_{00}=', ie3 (1, 1)
      ! Head
      scieff(1, 1) = w2 * ie3 (1, 1)
      ! Wings
      if(n .gt. 1) then
        forall(j1=2:n)
          scieff(j1, 1) = w1 * sptclg(j1, iq0) * ie3 (j1, 1)
          scieff(1, j1) = conjg(scieff(j1, 1))
        end forall
        ! Body
        forall(j1=2:n, j2=2:n)
          scieff(j1, j2) = sptclg(j1, iq0) * sptclg(j2, iq0) * ie3 (j1+2, j2+2)
        end forall
      end if
      deallocate(e3, ie3)

    case default
      write(*,*)
      write(*, '("Error(angavsc0): invalid averaging method")')
      write(*,*)
      stop

  end select

  if(n .gt. 1) deallocate(b, bi, u, v, s, t)

  call writedielt('DIELTENS', 1, [0.d0], dielten, 1)
  call writedielt('DIELTENS_NOSYM', 1, [0.d0], dtns, 1)

end subroutine angavsc0
!EOC
