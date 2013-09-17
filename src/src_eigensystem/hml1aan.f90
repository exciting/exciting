!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: hmlaa
! !INTERFACE:
!
!
Subroutine hml1aan (hamilton, is, ia, ngp, apwalm)
! !USES:
      Use modinput
      Use modmain
      Use modfvsystem
! !INPUT/OUTPUT PARAMETERS:
!   tapp   : .true. if the Hamiltonian is to be applied to the input vector,
!            .false. if the full matrix is to be calculated (in,logical)
!   is     : species number (in,integer)
!   ia     : atom number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   v      : input vector to which H is applied if tapp is .true., otherwise
!            not referenced (in,complex(nmatmax))
!   h      : H applied to v if tapp is .true., otherwise it is the Hamiltonian
!            matrix in packed form (inout,complex(npmatmax))
! !DESCRIPTION:
!   Calculates the APW-APW contribution to the Hamiltonian matrix.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
      Implicit None
!
! arguments
      Type (HermitianMatrix), Intent (Inout) :: hamilton
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Complex (8) :: x (ngp), y (ngp)
!
! local variables
      Integer :: ias, io1, io2
      Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
      Real (8) :: t1
      Complex (8) zt1, zsum, zt2
      Real(8) :: alpha,a2,energyref
      Parameter (alpha=1d0 / 137.03599911d0)

! automatic arrays
      Complex (8) zv (ngp)
! external functions
      Real (8) :: polynom
      Complex (8) zdotc
      External polynom, zdotc
!      write(*,*) gntyry (1, 1, 1)
!      stop0.5d0*alpha**2
      energyref=input%groundstate%energyref
      if (input%groundstate%ValenceRelativity.eq."scalar") then
         a2=0.5d0*alpha**2
      else
         a2=0d0
      endif
      t1 = 0.25d0 * rmt (is) ** 2
      ias = idxas (ia, is)
      Do l1 = 0, input%groundstate%lmaxmat
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do io1 = 1, apword (l1, is)
               zv (:) = 0.d0
               Do io2 = 1, apword (l1, is)
                  zsum = 0.d0
                  zt1 = 0.5d0*gntyry (lm1, 1, lm1) * h1aa (io1, io2, l1, ias)
!If (lm1 .Eq. lm3) zsum = zsum * 0.5d0
!0.5d0 factor is now above
                  zsum = zsum + zt1
                  Call zaxpy (ngp, zsum, apwalm(1, io2, lm1, ias), 1, zv, 1)
               End Do
!               Call zaxpy (ngp, zsum, apwalm(1, io2, lm1, ias), 1, zv, 1)
               x = conjg (apwalm(1:ngp, io1, lm1, ias))
               y = conjg (zv)
               Call Hermitianmatrix_rank2update (hamilton, ngp, zone, x, y)
            End Do
         End Do
      End Do
      Return
End Subroutine
!EOC
