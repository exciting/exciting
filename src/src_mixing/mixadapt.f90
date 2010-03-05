
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: mixadapt
! !INTERFACE:
!
!
Subroutine mixadapt (iscl, beta0, betainc, betadec, n, nu, mu, beta, f, &
& d)
! !INPUT/OUTPUT PARAMETERS:
!   iscl    : self-consistent loop number (in,integer)
!   beta0   : initial value for mixing parameter (in,real)
!   betainc : mixing parameter increase (in,real)
!   betadec : mixing parameter decrease (in,real)
!   n       : vector length (in,integer)
!   nu      : current output vector as well as the next input vector of the
!             self-consistent loop (inout,real(n))
!   mu      : used internally (inout,real(n))
!   beta    : used internally (inout,real(n))
!   f       : used internally (inout,real(n))
!   d       : RMS difference between old and new output vectors (out,real)
! !DESCRIPTION:
!   Given the input vector $\mu^i$ and output vector $\nu^i$ of the $i$th
!   self-consistent loop, this routine generates the next input vector to the
!   loop using an adaptive mixing scheme. The $j$th component of the output
!   vector is mixed with a fraction of the same component of the input vector:
!   $$ \mu^{i+1}_j=\beta^i_j\nu^i_j+(1-\beta^i_j)\mu^i_j, $$
!   where $\beta^i_j$ is set to $\beta_0$ at initialisation and increased by
!   scaling with $\beta_{\rm inc}$ ($>1$) if $f^i_j\equiv\nu^i_j-\mu^i_j$ does
!   not change sign between loops. If $f^i_j$ does change sign, then $\beta^i_j$
!   is scaled by $\beta_{\rm dec}$ ($<1$). Note that the array {\tt nu} serves
!   for both input and output, and the arrays {\tt mu}, {\tt beta} and {\tt f}
!   are used internally and should not be changed between calls. The routine is
!   initialised at the first iteration and is thread-safe so long as each thread
!   has its own independent work array. Complex arrays may be passed as real
!   arrays with $n$ doubled.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!   Modified, September 2008 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: iscl
      Real (8), Intent (In) :: beta0
      Real (8), Intent (In) :: betainc
      Real (8), Intent (In) :: betadec
      Integer, Intent (In) :: n
      Real (8), Intent (Inout) :: nu (n)
      Real (8), Intent (Inout) :: mu (n)
      Real (8), Intent (Inout) :: beta (n)
      Real (8), Intent (Inout) :: f (n)
      Real (8), Intent (Out) :: d
! local variables
      Integer :: i
      Real (8) :: t1
!
      d = 0.d0
      Do i = 1, n
         d = d + (nu(i)-mu(i)) ** 2
      End Do
      d = Sqrt (d/dble(n))
!
!
      If (iscl .Le. 1) Then
         mu (:) = nu (:)
         f (:) = 0.d0
         beta (:) = beta0
         d = 1.d0
         Return
      End If
      Do i = 1, n
         t1 = nu (i) - mu (i)
         If (t1*f(i) .Gt. 0.d0) Then
            beta (i) = beta (i) * betainc
            If (beta(i) .Gt. 1.d0) beta (i) = 1.d0
         Else
            beta (i) = beta (i) * betadec
         End If
         f (i) = t1
      End Do
      nu (:) = beta (:) * nu (:) + (1.d0-beta(:)) * mu (:)
!
      mu (:) = nu (:)
      Return
End Subroutine
!EOC
