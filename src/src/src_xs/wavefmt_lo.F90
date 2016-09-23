!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: wavefmt_lo
! !INTERFACE:
!
!
Subroutine wavefmt_lo (lrstp, lmax, is, ia, ngp, evecfv, ld, wfmt)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp  : radial step length (in,integer)
!   lmax   : maximum angular momentum required (in,integer)
!   is     : species number (in,integer)
!   ia     : atom number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   evecfv : first-variational eigenvector (in,complex(nmatmax))
!   ld     : leading dimension (in,integer)
!   wfmt   : muffin-tin wavefunction (out,complex(ld,*))
! !DESCRIPTION:
!   Muffin-tin wavefunction built up by local orbital contribution only.
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
      Complex (8), Intent (In) :: evecfv (nmatmax)
      Integer, Intent (In) :: ld
      Complex (8), Intent (Out) :: wfmt (ld,*)
! local variables
      Integer :: ias, l, m, lm, i
      Integer :: ir, nr, ilo
      Real (8) :: a, b
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
              & lofr(1, 1, ilo, ias))
            End Do
         End If
      End Do
      Return
End Subroutine
!EOC
