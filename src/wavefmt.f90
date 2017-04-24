! NOTE "FANCY" CODE: 
!   Not compilable with debug options because a complex array is passed to a real array dummy argument.
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: wavefmt
! !INTERFACE:
!
!
Subroutine wavefmt (lrstp, lmax, is, ia, ngp, apwalm, evecfv, ld, wfmt)
! !USES:
      Use modinput
      use mod_atoms, only: idxas, natmtot
      use mod_muffin_tin, only: idxlm, nrmt, lmmaxapw
      use mod_eigensystem, only: idxlo, nmatmax_ptr
      use mod_APW_LO, only: apword, apwordmax, apwfr,&
                          & lofr, nlorb, lorbl
      use mod_Gkvector, only: ngkmax_ptr
! !INPUT/OUTPUT PARAMETERS:
!   lrstp  : radial step length (in,integer)
!   lmax   : maximum angular momentum required (in,integer)
!   is     : species number (in,integer)
!   ia     : atom number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evecfv : first-variational eigenvector (in,complex(nmatmax))
!   ld     : leading dimension (in,integer)
!   wfmt   : muffin-tin wavefunction (out,complex(ld,*))
! !DESCRIPTION:
!   Calculates the first-variational wavefunction in the muffin-tin in terms of
!   a spherical harmonic expansion. For atom $\alpha$ and a particular $k$-point
!   ${\bf p}$, the $r$-dependent $(l,m)$-coefficients of the wavefunction for
!   the $i$th state are given by
!   $$ \Phi^{i{\bf p}}_{\alpha lm}(r)=\sum_{\bf G}b^{i{\bf p}}_{\bf G}
!    \sum_{j=1}^{M^{\alpha}_l}A^{\alpha}_{jlm}({\bf G+p})u^{\alpha}_{jl}(r)
!    +\sum_{j=1}^{N^{\alpha}}b^{i{\bf p}}_{(\alpha,j,m)}v^{\alpha}_j(r)
!    \delta_{l,l_j}, $$
!   where $b^{i{\bf p}}$ is the $i$th eigenvector returned from routine
!   {\tt seceqn}; $A^{\alpha}_{jlm}({\bf G+p})$ is the matching coefficient;
!   $M^{\alpha}_l$ is the order of the APW; $u^{\alpha}_{jl}$ is the APW radial
!   function; $N^{\alpha}$ is the number of local-orbitals; $v^{\alpha}_j$ is
!   the $j$th local-orbital radial function; and $(\alpha,j,m)$ is a compound
!   index for the location of the local-orbital in the eigenvector. See routines
!   {\tt genapwfr}, {\tt genlofr}, {\tt match} and {\tt seceqn}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Fixed description, October 2004 (C. Brouder)
!   Removed argument ist, November 2006 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lrstp
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax_ptr, apwordmax, lmmaxapw, &
     & natmtot)
      Complex (8), Intent (In) :: evecfv (nmatmax_ptr)
      Integer, Intent (In) :: ld
      Complex (8), Intent (Out) :: wfmt (ld,*)
! local variables
      Integer :: ias, l, m, lm, i
      Integer :: ir, nr, io, ilo
      Real (8) :: a, b
      Complex (8) zt1
! external functions
      Complex (8) zdotu
      External zdotu
      If (lmax .Gt. input%groundstate%lmaxapw) Then
         Write (*,*)
         Write (*, '("Error(wavefmt): lmax > lmaxapw : ", I8)') lmax
         Write (*,*)
         Stop
      End If
      ias = idxas (ia, is)
! zero the wavefunction
      nr = 0
      Do ir = 1, nrmt (is), lrstp
         nr = nr + 1
         wfmt (:, nr) = 0.d0
      End Do
! APW functions
      Do l = 0, lmax
         Do m = - l, l
            lm = idxlm (l, m)
            Do io = 1, apword (l, is)
               zt1 = zdotu (ngp, evecfv, 1, apwalm(:, io, lm, ias), 1)
               a = dble (zt1)
               b = aimag (zt1)
               Call wavefmt_add (nr, ld, wfmt(lm, 1), a, b, lrstp, &
              & apwfr(:, :, io, l, ias))
            End Do
         End Do
      End Do
! local-orbital functions
      Do ilo = 1, nlorb (is)
         l = lorbl (ilo, is)
         If (l .Le. lmax) Then
            Do m = - l, l
               lm = idxlm (l, m)
               i = ngp + idxlo (lm, ilo, ias)
               a = dble (evecfv(i))
               b = aimag (evecfv(i))
               Call wavefmt_add (nr, ld, wfmt(lm, 1), a, b, lrstp, &
              & lofr(:, :, ilo, ias))
            End Do
         End If
      End Do
      Return
End Subroutine
!EOC
!
!BOP
! !ROUTINE: wavefmt_add
! !INTERFACE:
!
!
Subroutine wavefmt_add (nr, ld, wfmt, a, b, lrstp, fr)
! !INPUT/OUTPUT PARAMETERS:
!   nr     : number of radial mesh points (in,integer)
!   ld     : leading dimension (in,integer)
!   wfmt   : complex muffin-tin wavefunction passed in as a real array
!            (inout,real(2*ld,*))
!   a      : real part of complex constant (in,real)
!   b      : imaginary part of complex constant (in,real)
!   lrstp  : radial step length (in,integer)
!   fr     : real radial function (in,real(lrstp,*))
! !DESCRIPTION:
!   Adds a real function times a complex constant to a complex muffin-tin
!   wavefunction as efficiently as possible. See routine {\tt wavefmt}.
!
! !REVISION HISTORY:
!   Created December 2006 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: nr
      Integer, Intent (In) :: ld
      Real (8), Intent (Inout) :: wfmt (2*ld,*)
      Real (8), Intent (In) :: a
      Real (8), Intent (In) :: b
      Integer, Intent (In) :: lrstp
      Real (8), Intent (In) :: fr (lrstp,*)
! local variables
      Integer :: ir
! values smaller than eps are taken to be zero
      Real (8), Parameter :: eps = 1.d-14
      If (Abs(b) .Lt. eps) Then
! zero constant
         If (Abs(a) .Lt. eps) Return
! pure real constant
         Do ir = 1, nr
            wfmt (1, ir) = wfmt (1, ir) + a * fr (1, ir)
         End Do
      Else If (Abs(a) .Lt. eps) Then
! pure imaginary constant
         Do ir = 1, nr
            wfmt (2, ir) = wfmt (2, ir) + b * fr (1, ir)
         End Do
      Else
! general complex constant
         Do ir = 1, nr
            wfmt (1, ir) = wfmt (1, ir) + a * fr (1, ir)
            wfmt (2, ir) = wfmt (2, ir) + b * fr (1, ir)
         End Do
      End If
      Return
End Subroutine
!EOC
