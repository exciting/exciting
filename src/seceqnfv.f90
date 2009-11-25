!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: seceqnfv
! !INTERFACE:
!
!
Subroutine seceqnfv (nmatp, ngp, igpig, vgpc, apwalm, evalfv, evecfv)
  ! !USES:
      Use modinput
      Use modmain
      Use modfvsystem
!
  ! !INPUT/OUTPUT PARAMETERS:
  !   nmatp  : order of overlap and Hamiltonian matrices (in,integer)
  !   ngp    : number of G+k-vectors for augmented plane waves (in,integer)
  !   igpig  : index from G+k-vectors to G-vectors (in,integer(ngkmax))
  !   vgpc   : G+k-vectors in Cartesian coordinates (in,real(3,ngkmax))
  !   apwalm : APW matching coefficients
  !            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
  !   evalfv : first-variational eigenvalues (out,real(nstfv))
  !   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
  ! !DESCRIPTION:
  !   Solves the secular equation,
  !   $$ (H-\epsilon O)b=0, $$
  !   for the all the first-variational states of the input $k$-point.
  !
  ! !REVISION HISTORY:
  !   Created March 2004 (JKD)
  !EOP
  !BOC
      Implicit None
  ! arguments
      Integer, Intent (In) :: nmatp
      Integer, Intent (In) :: ngp
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Real (8), Intent (Out) :: evalfv (nstfv)
      Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv)
  ! local variables
      Type (evsystem) :: system
      Logical :: packed
      Integer :: is, ia, i, m, np, info
      Real (8) :: vl, vu
      Real (8) :: ts0, ts1
  ! allocatable arrays
      Integer, Allocatable :: iwork (:)
      Integer, Allocatable :: ifail (:)
      Real (8), Allocatable :: w (:)
      Real (8), Allocatable :: rwork (:)
      Complex (8), Allocatable :: v (:)
      Complex (8), Allocatable :: work (:)
      np = (nmatp*(nmatp+1)) / 2
!
  !----------------------------------------!
  !     Hamiltonian and overlap set up     !
  !----------------------------------------!
!
      packed = .True.
      Call newsystem (system, packed, nmatp)
      Call hamiltonandoverlapsetup (system, ngp, apwalm, igpig, vgpc)
!
  !------------------------------------!
  !     solve the secular equation     !
  !------------------------------------!
!
      Call timesec (ts0)
      vl = 0.d0
      vu = 0.d0
  ! LAPACK 3.0 call
!
      Allocate (iwork(5*nmatp))
      Allocate (ifail(nmatp))
      Allocate (w(nmatp))
      Allocate (rwork(7*nmatp))
      Allocate (v(1))
      Allocate (work(2*nmatp))
      Call zhpgvx (1, 'V', 'I', 'U', nmatp, system%hamilton%zap, &
     & system%overlap%zap, vl, vu, 1, nstfv, &
     & input%groundstate%solver%evaltol, m, w, evecfv, nmatmax, work, &
     & rwork, iwork, ifail, info)
      evalfv (1:nstfv) = w (1:nstfv)
!
!
!
      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(seceqnfv): diagonalisation failed")')
         Write (*, '(" ZHPGVX returned INFO = ", I8)') info
         If (info .Gt. nmatp) Then
            i = info - nmatp
            Write (*, '(" The leading minor of the overlap matrix of or&
           &der ", I8)') i
            Write (*, '("  is not positive definite")')
            Write (*, '(" Order of overlap matrix : ", I8)') nmatp
            Write (*,*)
         End If
         Stop
      End If
      Call timesec (ts1)
  !$OMP CRITICAL
      timefv = timefv + ts1 - ts0
  !$OMP END CRITICAL
      Call deleteystem (system)
      Deallocate (iwork, ifail, w, rwork, v, work)
      Return
End Subroutine seceqnfv
!EOC
