
! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine angavsc0 (n, nmax, scrnh, scrnw, scrn, scieff)
      Use modmain
      Use modinput
      Use modxs
      Use invert
      Implicit None
  ! arguments
      Integer, Intent (In) :: n, nmax
      Complex (8), Intent (In) :: scrn (n, n), scrnw (n, 2, 3), scrnh &
     & (3, 3)
      Complex (8), Intent (Out) :: scieff (nmax, nmax)
  ! local variables
      Integer, Parameter :: nsphcov = 5810, iq0 = 1
      Integer :: iop, jop, j1, j2, ji, itp, lm, ntpsph
      Real (8) :: vomega, t00, r, qsz
      Complex (8) :: dtns (3, 3), w1, w2
      Real (8), Allocatable :: plat (:, :), p (:), tp (:, :), spc (:, &
     & :), w (:)
      Complex (8), Allocatable :: m00lm (:), mx0lm (:), mxxlm (:)
      Complex (8), Allocatable :: ei00 (:), eix0 (:), ei0x (:), eixx &
     & (:)
      Complex (8), Allocatable :: ei00lm (:), eix0lm (:), ei0xlm (:), &
     & eixxlm (:)
      Complex (8), Allocatable :: ylm (:), zylm (:, :)
      Complex (8), Allocatable :: b (:, :), bi (:, :), u (:, :), v (:, &
     & :), s (:, :), t (:, :), e3 (:, :), ie3 (:, :)
      Integer :: i1, i2
  ! crystal volume
      vomega = omega * product (ngridq)
  ! Wigner-Seitz radius and spherical approximation to 1/q^2 average
      qsz = (6*pi**2/vomega) ** (1.d0/3.d0)
  ! weight for sqrt(4pi)/q based on Wigner-Seitz radius
      w1 = qsz ** 2 * vomega / (4.d0*pi**2) * Sqrt (fourpi)
  ! weight for 4pi/q^2 based on Wigner-Seitz radius
      w2 = 2 * qsz * vomega / pi
  ! calculate RPA dielectric tensor including local field effects
      dielten0 (:, :) = scrnh (:, :)
      If (n .Gt. 1) Then
         Allocate (b(n-1, n-1), bi(n-1, n-1), u(n-1, 3), v(3, n-1), &
        & s(n-1, 3), t(3, n-1))
     ! body of dielectric matrix
         b (:, :) = scrn (2:, 2:)
     ! column wing
         u (:, :) = scrnw (2:, 2, :)
     ! row wing
         v (:, :) = transpose (scrnw(2:, 1, :))
         Select Case (input%xs%BSE%scrherm)
         Case (0)
        ! use full matrix (both wings and full body)
         Case (1)
        ! Hermitian average matrix (average both wings and body)
            b = 0.5d0 * (b+conjg(transpose(b)))
            u = 0.5d0 * (u+conjg(transpose(v)))
            v = conjg (transpose(u))
         Case (2)
        ! use upper triangle (take row wing, assign column wing)
            u = conjg (transpose(v))
         Case (3)
        ! use lower triangle (take column wing, assign row wing)
            v = conjg (transpose(u))
         Case Default
            Write (*,*)
            Write (*, '("Error(angavsc0): not a valid flag:", i6)') &
           & input%xs%BSE%scrherm
            Write (*,*)
            Call terminate
         End Select
     ! invert body (optionally including Hermitian average)
         Call zinvert_hermitian (input%xs%BSE%scrherm, b, bi)
         s = matmul (bi, u)
         t = matmul (v, bi)
         dielten = dielten0 - matmul (v, s)
      Else
         dielten = dielten0
      End If
  ! symmetrize the dielectric tensor
      dtns (:, :) = dielten (:, :)
      Do iop = 1, 3
         Do jop = 1, 3
            Call symt2app (iop, jop, 1, symt2, dtns, dielten(iop, jop))
         End Do
      End Do
  ! calculate averaged screened Coulomb interaction in Fourier space at Gamma point
      Select Case (trim(input%xs%BSE%sciavtype))
      Case ('spherical')
     ! scaling factor
         t00 = (omega/(twopi)**3) * product (ngridq)
     ! number of points on sphere
         If (tleblaik) Then
            ntpsph = input%xs%BSE%nleblaik
         Else
            ntpsph = nsphcov
         End If
         If (lmmaxdielt .Gt. ntpsph) Then
            Write (*,*)
            Write (*, '("Error(angavdm0): lmmaxdielt.gt.ntpsph: ", 2i6)&
           &') lmmaxdielt, ntpsph
            Write (*,*)
            Stop
         End If
         Allocate (plat(3, ntpsph), p(ntpsph))
         Allocate (m00lm(lmmaxdielt), mx0lm(lmmaxdielt), &
        & mxxlm(lmmaxdielt))
         Allocate (ei00(ntpsph), eix0(ntpsph), ei0x(ntpsph), &
        & eixx(ntpsph))
         Allocate (ei00lm(lmmaxdielt), eix0lm(lmmaxdielt), &
        & ei0xlm(lmmaxdielt), eixxlm(lmmaxdielt))
         Allocate (ylm(lmmaxdielt), zylm(ntpsph, lmmaxdielt))
         Allocate (tp(2, ntpsph), spc(3, ntpsph))
         Allocate (w(ntpsph))
         If (tleblaik) Then
        ! generate Lebedev Laikov grid
            Call leblaik (ntpsph, spc, w)
        ! generate tetha and phi angles
            Do itp = 1, ntpsph
               Call sphcrd (spc(:, itp), r, tp(:, itp))
            End Do
         Else
        ! distribution is assumed to be uniform
            w (:) = 1.d0 / ntpsph
        ! generate spherical covering set (angles and coordinates)
            Call sphcover (ntpsph, tp)
            spc (1, :) = Sin (tp(1, :)) * Cos (tp(2, :))
            spc (2, :) = Sin (tp(1, :)) * Sin (tp(2, :))
            spc (3, :) = Cos (tp(1, :))
         End If
     ! generate spherical harmonics on covering set
         Do itp = 1, ntpsph
            Call genylm (input%xs%BSE%lmaxdielt, tp(:, itp), ylm)
            zylm (itp, :) = ylm (:)
         End Do
     ! unit vectors of spherical covering set in lattice coordinates
         plat = matmul (binv, spc)
     ! distances to subcell cell boundaries in reciprocal space
         Do itp = 1, ntpsph
            p (itp) = 1.d0 / (2.d0*maxval(Abs(ngridq(:)*plat(:, itp)), &
           & 1))
         End Do
     ! calculate function on covering set
         Do itp = 1, ntpsph
        ! head, 1/(p*L*p)
            ei00 (itp) = 1.d0 / dot_product (spc(:, itp), &
           & matmul(dielten, spc(:, itp)))
         End Do
     ! calculate lm-expansion coefficients
         Do lm = 1, lmmaxdielt
            ei00lm (lm) = fourpi * dot_product (zylm(:, lm), ei00*w)
            m00lm (lm) = fourpi * dot_product (zylm(:, lm), p*w)
            mx0lm (lm) = fourpi * dot_product (zylm(:, lm), &
           & p**2/2.d0*w)
            mxxlm (lm) = fourpi * dot_product (zylm(:, lm), &
           & p**3/3.d0*w)
         End Do
     ! subcell average (head)
         scieff (1, 1) = fourpi * t00 * dot_product (m00lm, ei00lm)
     ! loop over (G,Gp) indices
         Do j1 = 2, n
            Do itp = 1, ntpsph
           ! wing, -p*S/(p*L*p)
               eix0 (itp) = - dot_product (spc(:, itp), s(j1-1, :)) * &
              & ei00 (itp)
           ! wing, -p*T/(p*L*p)
               If (input%xs%BSE%scrherm .Eq. 0) ei0x (itp) = - &
              & dot_product (spc(:, itp), t(:, j1-1)) * ei00 (itp)
            End Do
            Do lm = 1, lmmaxdielt
               eix0lm (lm) = fourpi * dot_product (zylm(:, lm), eix0*w)
               If (input%xs%BSE%scrherm .Eq. 0) ei0xlm (lm) = fourpi * &
              & dot_product (zylm(:, lm), ei0x*w)
            End Do
        ! subcell average (wings)
            scieff (j1, 1) = Sqrt (fourpi) * sptclg (j1, iq0) * t00 * &
           & dot_product (mx0lm, eix0lm)
            If (input%xs%BSE%scrherm .Eq. 0) Then
               scieff (1, j1) = Sqrt (fourpi) * sptclg (j1, iq0) * t00 &
              & * dot_product (mx0lm, ei0xlm)
            Else
               scieff (1, j1) = conjg (scieff(j1, 1))
            End If
            If (input%xs%BSE%sciavbd) Then
               ji = j1
               If (input%xs%BSE%scrherm .Eq. 0) ji = 2
               Do j2 = ji, n
                  Do itp = 1, ntpsph
                 ! body, B^-1 + p*S p*T/(p*L*p)
                     eixx (itp) = bi (j1-1, j2-1) + dot_product (spc(:, &
                    & itp), s(j1-1, :)) * dot_product (spc(:, itp), &
                    & t(:, j2-1)) * ei00 (itp)
                  End Do
                  Do lm = 1, lmmaxdielt
                     eixxlm (lm) = fourpi * dot_product (zylm(:, lm), &
                    & eixx*w)
                  End Do
              ! subcell average (body)
                  scieff (j1, j2) = sptclg (j1, iq0) * sptclg (j2, iq0) &
                 & * t00 * dot_product (mxxlm, eixxlm)
                  If (input%xs%BSE%scrherm .Ne. 0) scieff (j2, j1) = &
                 & conjg (scieff(j1, j2))
               End Do
            Else
           ! no subcell average (body)
               scieff (j1, 2:n) = bi (j1-1, :)
            End If
         End Do
         Deallocate (ei00, eix0, ei0x, eixx, ei00lm, eix0lm, ei0xlm, &
        & m00lm, mx0lm, mxxlm)
         Deallocate (ylm, zylm, tp, spc, w, plat, p)
      Case ('screendiag', 'invscreendiag')
         If (input%xs%BSE%sciavbd) Then
            Write (*,*)
            Write (*, '("Error(angavsc0): (inv)screendiag-method does n&
           &ot allow for averaging the body of W")')
            Write (*,*)
            Stop
         End If
         Allocate (e3(n+2, n+2), ie3(n+2, n+2))
     ! invert dielectric matrix including 3 times G=0 according to the limits
     ! q->0_x, q->0_y, q->0_z
     ! G=0, G'=0 elements
         e3 (1:3, 1:3) = dielten0 (:, :)
     ! G!=0, G'=0 components and vice versa
         If (n .Gt. 1) Then
            Do i1 = 1, 3
               Do j2 = 2, n
                  e3 (i1, j2+2) = scrnw (j2, 1, i1)
               End Do
            End Do
            Do j1 = 2, n
               Do i2 = 1, 3
                  e3 (j1+2, i2) = scrnw (j1, 2, i2)
               End Do
            End Do
            Do j1 = 2, n
               Do j2 = 2, n
                  e3 (j1+2, j2+2) = scrn (j1, j2)
               End Do
            End Do
         End If
         Call zinvert_hermitian (input%xs%BSE%scrherm, e3, ie3)
     ! select again
         Select Case (trim(input%xs%BSE%sciavtype))
         Case ('screendiag')
        ! head
            scieff (1, 1) = w2 * 1.d0 / ((e3(1, 1)+ie3(2, 2)+ie3(3, &
           & 3))/3.d0)
            If (n .Gt. 1) Then
           ! wings, set to zero in this approximation
               scieff (1, 2:n) = zzero
               scieff (2:n, 1) = zzero
           ! body, only diagonal is assigned
               scieff (2:n, 2:n) = zzero
               Forall (j1=2:n)
                  scieff (j1, j1) = sptclg (j1, iq0) ** 2 / e3 (j1+2, &
                 & j1+2)
               End Forall
            End If
         Case ('invscreendiag')
        ! head
            scieff (1, 1) = w2 * (ie3(1, 1)+ie3(2, 2)+ie3(3, 3)) / 3.d0
        ! wings
            If (n .Gt. 1) Then
               Forall (j1=2:n)
                  scieff (j1, 1) = w1 * sptclg (j1, iq0) * (ie3(j1+2, &
                 & 1)+ie3(j1+2, 2)+ie3(j1+3, 3)) / 3.d0
                  scieff (1, j1) = conjg (scieff(j1, 1))
               End Forall
           ! body
               Forall (j1=2:n, j2=2:n)
                  scieff (j1, j2) = sptclg (j1, iq0) * sptclg (j2, iq0) &
                 & * ie3 (j1+2, j2+2)
               End Forall
            End If
         End Select
         Deallocate (e3, ie3)
      Case ('none')
         iop = 1 !!!only x-component here!!!
     ! longitudinal treatment, three components of vanishing q (direction)
         Allocate (e3(n, n), ie3(n, n))
         e3 (1, 1) = scrnh (iop, iop)
         e3 (1, 2:n) = scrnw (2:, 1, iop)
         e3 (2:n, 1) = scrnw (2:, 2, iop)
         e3 (2:n, 2:n) = scrn (2:, 2:)
         Call zinvert_hermitian (input%xs%BSE%scrherm, e3, ie3)
         Write (*,*) 'eps^{-1}_{00}=', ie3 (1, 1)
     ! head
         scieff (1, 1) = w2 * ie3 (1, 1)
     ! wings
         If (n .Gt. 1) Then
            Forall (j1=2:n)
               scieff (j1, 1) = w1 * sptclg (j1, iq0) * ie3 (j1, 1)
               scieff (1, j1) = conjg (scieff(j1, 1))
            End Forall
        ! body
            Forall (j1=2:n, j2=2:n)
               scieff (j1, j2) = sptclg (j1, iq0) * sptclg (j2, iq0) * &
              & ie3 (j1+2, j2+2)
            End Forall
         End If
         Deallocate (e3, ie3)
      Case Default
         Write (*,*)
         Write (*, '("Error(angavsc0): invalid averaging method")')
         Write (*,*)
         Stop
      End Select

      If (n .Gt. 1) deallocate (b, bi, u, v, s, t)

      Call writedielt ('DIELTENS', 1, 0.d0, dielten, 1)
      Call writedielt ('DIELTENS_NOSYM', 1, 0.d0, dtns, 1)

End Subroutine angavsc0
