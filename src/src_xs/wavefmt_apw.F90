!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: wavefmt_apw
! !INTERFACE:
!
!
Subroutine wavefmt_apw (lrstp, lmax, is, ia, ngp, apwalm, evecfv, ld, &
& wfmt)
! !USES:
      Use modinput
      Use modmain
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
!   Muffin-tin wavefunction built up by APW contribution only.
!   Based upon the routine {\tt wavefmt}.
!
! !REVISION HISTORY:
!   Created May 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lrstp
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Complex (8), Intent (In) :: evecfv (nmatmax)
      Integer, Intent (In) :: ld
      Complex (8), Intent (Out) :: wfmt (ld,*)
! local variables
      Integer :: ias, l, m, lm
      Integer :: ir, nr, io
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
               zt1 = zdotu (ngp, evecfv, 1, apwalm(1, io, lm, ias), 1)
               a = dble (zt1)
               b = aimag (zt1)
               Call wavefmt_add (nr, ld, wfmt(lm, 1), a, b, lrstp, &
              & apwfr(1, 1, io, l, ias))
            End Do
         End Do
      End Do
      Return
End Subroutine
!EOC
