!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: olprad
! !INTERFACE:
!
!
Subroutine olprad
! !USES:
      Use modmain
! !DESCRIPTION:
!   Calculates the radial overlap integrals of the APW and local-orbital basis
!   functions. In other words, for spin $\sigma$ and atom $j$ of species $i$, it
!   computes integrals of the form
!   $$ o^{\sigma;ij}_{qp}=\int_0^{R_i}u^{\sigma;ij}_{q;l_p}(r)v^{\sigma;ij}_p(r)
!    r^2dr $$
!   and
!   $$ o^{\sigma;ij}_{pp'}=\int_0^{R_i}v^{\sigma;ij}_p(r)v^{\sigma;ij}_{p'}(r)
!    r^2dr,\quad l_p=l_{p'} $$
!   where $u^{\sigma;ij}_{q;l}$ is the $q$th APW radial function for angular
!   momentum $l$; and $v^{\sigma;ij}_p$ is the $p$th local-orbital radial
!   function and has angular momentum $l_p$.
!
! !REVISION HISTORY:
!   Created November 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, ir, nr
      Integer :: l, ilo, ilo1, ilo2, io
! automatic arrays
      Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax)
      Do is = 1, nspecies
         nr = nrmt (is)
         Do ir = 1, nr
            r2 (ir) = spr (ir, is) ** 2
         End Do
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
!--------------------------------------!
!     APW-local-orbital integtrals     !
!--------------------------------------!
            Do ilo = 1, nlorb (is)
               l = lorbl (ilo, is)
               Do io = 1, apword (l, is)
                  Do ir = 1, nr
                     fr (ir) = apwfr (ir, 1, io, l, ias) * lofr (ir, 1, &
                    & ilo, ias) * r2 (ir)
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  oalo (io, ilo, ias) = gr (nr)
               End Do
            End Do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
            Do ilo1 = 1, nlorb (is)
               l = lorbl (ilo1, is)
               Do ilo2 = 1, nlorb (is)
                  If (lorbl(ilo2, is) .Eq. l) Then
                     Do ir = 1, nr
                        fr (ir) = lofr (ir, 1, ilo1, ias) * lofr (ir, &
                       & 1, ilo2, ias) * r2 (ir)
                     End Do
                     Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                     ololo (ilo1, ilo2, ias) = gr (nr)
                  End If
               End Do
            End Do
         End Do
      End Do
      Return
End Subroutine
!EOC
