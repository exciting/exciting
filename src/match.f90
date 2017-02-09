!
!
!
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: match
! !INTERFACE:
!
!
Subroutine match (ngp, gpc, tpgpc, sfacgp, apwalm)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer)
!   gpc    : length of G+p-vectors (in,real(ngkmax))
!   tpgpc  : (theta, phi) coordinates of G+p-vectors (in,real(2,ngkmax))
!   sfacgp : structure factors of G+p-vectors (in,complex(ngkmax,natmtot))
!   apwalm : APW matching coefficients
!            (out,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
! !DESCRIPTION:
!   Computes the $({\bf G+p})$-dependent matching coefficients for the APW basis
!   functions. Inside muffin-tin $\alpha$, the APW functions are given by
!   $$ \phi^{\alpha}_{\bf G+p}({\bf r})=\sum_{l=0}^{l_{\rm max}}
!    \sum_{m=-l}^{l}\sum_{j=1}^{M^{\alpha}_l}A^{\alpha}_{jlm}({\bf G+p})
!    u^{\alpha}_{jl}(r)Y_{lm}(\hat{{\bf r}}), $$
!   where $A^{\alpha}_{jlm}({\bf G+p})$ is the matching coefficient,
!   $M^{\alpha}_l$ is the order of the APW and $u^{\alpha}_{jl}$ is the radial
!   function. In the interstitial region, an APW function is a plane wave,
!   $\exp(i({\bf G+p})\cdot{\bf r})/\sqrt{\Omega}$, where $\Omega$ is the unit
!   cell volume. Ensuring continuity up to the $(M^{\alpha}_l-1)$th derivative
!   across the muffin-tin boundary therefore requires that the matching
!   coefficients satisfy
!   $$ \sum_{j=1}^{M^{\alpha}_l}D_{ij}A^{\alpha}_{jlm}({\bf G+p})=b_i\;, $$
!   where
!   $$ D_{ij}=\left.\frac{d^{i-1}u^{\alpha}_{jl}(r)}{dr^{i-1}}
!    \right|_{r=R_{\alpha}} $$
!   and
!   $$ b_i=\frac{4\pi i^l}{\sqrt{\Omega}}|{\bf G+p}|^{i-1}j^{(i-1)}_l
!    (|{\bf G+p}|R_{\alpha})\exp(i({\bf G+p})\cdot{\bf r}_{\alpha})Y^*_{lm}
!    (\widehat{{\bf G+p}}), $$
!   with ${\bf r}_{\alpha}$ the atomic position and $R_{\alpha}$ the muffin-tin
!   radius. See routine {\tt wavefmt}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Fixed documentation, June 2006 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: ngp
      Real (8), Intent (In) :: gpc (ngkmax)
      Real (8), Intent (In) :: tpgpc (2, ngkmax)
      Complex (8), Intent (In) :: sfacgp (ngkmax, natmtot)
      Complex (8), Intent (Out) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
! local variables
      Integer :: np, is, ia, ias, omax
      Integer :: l, m, lm, io1, io2
      Integer :: i, ir, igp, info
      Real (8) :: fpso, t1
      Complex (8) zt1, zt2
! allocatable arrays
      Integer, Allocatable :: ipiv (:)
      Real (8), Allocatable :: c (:)
      Real (8), Allocatable :: djl (:, :, :)
      Complex (8), Allocatable :: ylmgp (:, :)
      Complex (8), Allocatable :: zd (:, :)
      Complex (8), Allocatable :: zb (:, :)
! external functions
      Real (8) :: polynom
      External polynom
      fpso = fourpi / Sqrt (omega)
! polynomial order
      np = Max (apwordmax+1, 4)
      Allocate (ipiv(np))
      Allocate (c(np))
      Allocate (djl(0:input%groundstate%lmaxapw, apwordmax, ngp))
      Allocate (ylmgp(lmmaxapw, ngp))
      Allocate (zd(apwordmax, apwordmax))
      Allocate (zb(apwordmax, ngp*(2*input%groundstate%lmaxapw+1)))
! compute the spherical harmonics of the G+p-vectors
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(igp)
!$OMP DO
#endif
      Do igp = 1, ngp
         Call genylm (input%groundstate%lmaxapw, tpgpc(:, igp), &
        & ylmgp(:, igp))
      End Do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
! begin loops over atoms and species
      Do is = 1, nspecies
! evaluate the spherical Bessel function derivatives for all G+p-vectors
         omax = 0
         Do l = 0, input%groundstate%lmaxapw
            omax = Max (omax, apword(l, is))
         End Do
!
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(igp,t1,io1)
!$OMP DO
#endif
         Do igp = 1, ngp
            t1 = gpc (igp) * rmt (is)
            Do io1 = 1, omax
               Call sbesseldm (io1-1, input%groundstate%lmaxapw, t1, &
              & djl(:, io1, igp))
            End Do
            t1 = 1.d0
            Do io1 = 2, omax
               t1 = t1 * gpc (igp)
               djl (:, io1, igp) = t1 * djl (:, io1, igp)
            End Do
         End Do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
!
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! begin loop over l
            Do l = 0, input%groundstate%lmaxapw
               zt1 = fpso * zil (l)
! set up matrix of derivatives
!
               Do io2 = 1, apword (l, is)
                  ir = nrmt (is) - np + 1
                  Do io1 = 1, apword (l, is)
!
                     zd (io1, io2) = polynom (io1-1, np, spr(ir, is), &
                    & apwfr(ir, 1, io2, l, ias), c, rmt(is))
                  End Do
               End Do
!
               i = 0
               Do igp = 1, ngp
                  zt2 = zt1 * sfacgp (igp, ias)
                  Do m = - l, l
                     lm = idxlm (l, m)
                     i = i + 1
                     Do io1 = 1, apword (l, is)
                        zb (io1, i) = djl (l, io1, igp) * zt2 * conjg (ylmgp(lm, igp))
                     End Do
                  End Do
               End Do
! solve the general complex linear systems
!
               Call zgesv (apword(l, is), i, zd, apwordmax, ipiv, zb, &
              & apwordmax, info)
!
               If (info .Ne. 0) Then
                  Write (*,*)
                  Write (*, '("Error(match): could not find APW matchin&
                 &g coefficients")')
                  Write (*, '(" for species ", I4)') is
                  Write (*, '(" and atom ", I4)') ia
                  Write (*, '(" ZGESV returned INFO = ", I8)') info
                  Write (*,*)
       ! stop
               End If
               i = 0
               Do igp = 1, ngp
                  Do m = - l, l
                     lm = idxlm (l, m)
                     i = i + 1
                     Do io1 = 1, apword (l, is)
                        apwalm (igp, io1, lm, ias) = zb (io1, i)
                     End Do
                  End Do
               End Do
! end loop over l
            End Do
! end loops over atoms and species
         End Do
!
      End Do
      Deallocate (ipiv, c, djl, ylmgp, zd, zb)
      Return
End Subroutine
!EOC
