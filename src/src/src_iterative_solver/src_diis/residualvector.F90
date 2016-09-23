
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
!BOP
! !ROUTINE: seceqn
!
!
Subroutine residualvectors (n, iunconverged, h, s, evalfv, r, rnorms)
      Use modmain, Only: nmatmax, zone, zzero
!
! !INPUT/OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  the residual is
!  \newcommand{\bra}[1]{\langle #1|}
!\newcommand{\ket}[1]{|#1\rangle}
!\newcommand{\braket}[2]{\langle #1|#2\rangle}
!$$
!\ket{\mathbf{R}\left(\ket{\mathbf{A}^{ap}},E^{ap}\right)}=(\mathbf{H}-E^{ap}\mathbf{S})\ket{ \mathbf{A}^{ap}}
!$$
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
      Use modmain, Only: nstfv
      Use diisinterfaces
      Implicit None
      Integer, Intent (In) :: n, iunconverged
!packed ut
      Complex (8), Intent (In) :: h (n, nstfv), s (n, nstfv)
      Complex (8), Intent (Out) :: r (n, nstfv)
      Real (8), Intent (In) :: evalfv (nstfv)
      Real (8), Intent (Out) :: rnorms (nstfv)
      Integer :: i
      Complex (8) :: z
!
      Do i = 1, iunconverged
         Call zcopy (n, h(1, i), 1, r(1, i), 1)
         z = cmplx (-evalfv(i), 0)
         Call zaxpy (n, z, s(1, i), 1, r(1, i), 1)
         rnorms (i) = Sqrt (dble(zdotc(n, r(1, i), 1, r(1, i), 1)))
      End Do
!
End Subroutine
