!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: rschrodapp
! !INTERFACE:
!
!
Subroutine rschrodapp (l, nr, r, vr, p0, q0, q1, hp0)
use modinput
! !INPUT/OUTPUT PARAMETERS:
!   l   : angular momentum quantum number (in,integer)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   p0  : m th energy derivative of P (in,real(nr))
!   q0  : m th energy derivative of Q (in,real(nr))
!   q1  : radial derivative of q0 (in,real(nr))
!   hp0 : H applied to P (out,real(nr))
! !DESCRIPTION:
!   Applies the scalar relativistic radial Hamiltonian, $H$, to a radial
!   wavefunction, $P_l$. This is an approximation since we assume $P_l$ is a
!   scalar wavefunction, normalisable to unity. A Hamiltonian which satisfies
!   $H P_l=E P_l$ is given implicitly by
!   $$ H P_l=\left[\frac{l(l+1)}{2Mr^2}+V\right]P_l-\frac{1}{r}Q_l
!    -\frac{d}{dr}Q_l, $$
!   where $V$ is the external potential, $M=1-V/2c^2$ and $Q_l$ is obtained from
!   integrating the coupled scalar relativistic equations. See the routine
!   {\tt rschrodint} for further details.
!
! !REVISION HISTORY:
!   Created October 2003 (JKD)
!EOP
!BOC
      Implicit None
      Integer, Intent (In) :: l
      Integer, Intent (In) :: nr
      Real (8), Intent (In) :: r (nr)
      Real (8), Intent (In) :: vr (nr)
      Real (8), Intent (In) :: p0 (nr)
      Real (8), Intent (In) :: q0 (nr)
      Real (8), Intent (In) :: q1 (nr)
      Real (8), Intent (Out) :: hp0 (nr)
! local variables
      Integer :: ir
! fine structure constant
      Real (8), Parameter :: alpha = 1.d0 / 137.03599911d0
      Real (8) :: rm, t1, rmfactor
      if (input%groundstate%ValenceRelativity.eq."zora") then
         rmfactor=1d0
      else
         rmfactor=0d0
      endif
      Do ir = 1, nr
         rm = 1.d0 - 0.5d0 * (alpha**2) * vr (ir)*rmfactor
         t1 = dble (l*(l+1)) / (2.d0*rm*r(ir)**2)
         hp0 (ir) = (t1+vr(ir)) * p0 (ir) - q0 (ir) / r (ir) - q1 (ir)
      End Do
      Return
End Subroutine
!EOC
