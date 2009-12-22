!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_ftfun
      Implicit None
Contains
!
!
      Subroutine ftfun (ng, lmax, tir, tmt, gir, gmt, ftg)
         Use modmain
         Use modxs
         Implicit None
    ! arguments
         Integer, Intent (In) :: ng
         Integer, Intent (In) :: lmax
         Logical, Intent (In) :: tir
         Logical, Intent (In) :: tmt
         Complex (8), Intent (In) :: gir (:)
         Complex (8), Intent (In) :: gmt (:, :, :)
         Complex (8), Intent (Out) :: ftg (:)
    ! local variales
         Complex (8), Allocatable :: zfft (:)
         Complex (8) :: zt1, zt2
         Real (8), Allocatable :: jbesslh (:), jbessl (:, :)
         Real (8), Allocatable :: r1 (:), r2 (:), fr (:), fr2 (:), gr &
        & (:), cf (:, :)
         Integer :: ig, ifg, is, ia, ias, ir, nr, l, m, lm
!
    ! interstitial part
         ftg (:) = zzero
         If (tir) Then
            Allocate (zfft(ngrtot))
       ! multiply effective potential with smooth characteristic function
            zfft (:) = gir (:) * cfunir (:)
       ! Fourier transform to G-space
            Call zfftifc (3, ngrid,-1, zfft)
            Do ig = 1, ng
               ifg = igfft (ig)
               ftg (ig) = ftg (ig) + zfft (ifg)
            End Do
            Deallocate (zfft)
         End If
!
    ! muffin-tin part
         If (tmt) Then
       ! allocate array for Bessel functions
            Allocate (jbessl(0:lmax, nrmtmax))
            Allocate (jbesslh(0:lmax))
            Allocate (r1(nrmtmax), r2(nrmtmax), fr(nrmtmax), &
           & fr2(nrmtmax), gr(nrmtmax), cf(3, nrmtmax))
       ! loop over G vectors
            Do ig = 1, ng
          ! loop over species
               Do is = 1, nspecies
                  nr = nrmt (is)
                  Do ir = 1, nr
                     r1 (ir) = spr (ir, is)
                     r2 (ir) = r1 (ir) ** 2
                  End Do
             ! calculate bessel functions j_l(|G||r|)
                  Do ir = 1, nr
                     Call sbessel (lmax, gc(ig)*r1(ir), jbesslh)
                     jbessl (:, ir) = jbesslh (:)
                  End Do
             ! loop over atoms
                  Do ia = 1, natoms (is)
                     ias = idxas (ia, is)
                     zt1 = zzero
                     zt2 = zzero
                     Do l = 0, lmax
                        Do m = - l, l
                           lm = idxlm (l, m)
                           Do ir = 1, nr
                         ! integrand
                              fr (ir) = r2 (ir) * jbessl (l, ir) * dble &
                             & (gmt(lm, ir, ias))
                              fr2 (ir) = r2 (ir) * jbessl (l, ir) * &
                             & aimag (gmt(lm, ir, ias))
                           End Do
                           If (any(fr .Ne. 0.d0)) Then
                         ! integration
                              Call fderiv (-1, nr, spr(1, is), fr, gr, &
                             & cf)
                              zt1 = zt1 + gr (nr) * conjg (zil(l)) * &
                             & ylmg (lm, ig)
                           End If
                           If (any(fr2 .Ne. 0.d0)) Then
                         ! integration
                              Call fderiv (-1, nr, spr(1, is), fr2, gr, &
                             & cf)
                              zt2 = zt2 + gr (nr) * conjg (zil(l)) * &
                             & ylmg (lm, ig)
                           End If
                        End Do ! m
                     End Do ! l
                ! form factor summation
                     ftg (ig) = ftg (ig) + (fourpi/omega) * conjg &
                    & (sfacg(ig, ias)) * (zt1+zi*zt2)
                  End Do ! ia
               End Do ! is
            End Do ! ig
            Deallocate (jbessl, jbesslh, r1, r2, fr, fr2, gr, cf)
         End If
!
      End Subroutine ftfun
!
End Module m_ftfun
