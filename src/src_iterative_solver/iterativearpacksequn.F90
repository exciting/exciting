
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine iterativearpacksecequn (ik, ispn, apwalm, vgpc, evalfv, &
& evecfv)
!
      Use modinput
  !USES:
      Use modfvsystem
      Use modmpi
      Use mod_eigensystem
      Use mod_timing
      Use mod_Gkvector
      Use mod_potential_and_density
      Use mod_muffin_tin
      Use mod_atoms, Only: natmtot
!
      Use mod_spin, Only: nspnfv
      Use mod_APW_LO, Only: apwordmax
      Use mod_eigenvalue_occupancy, Only: nstfv
      Use sclcontroll
!
  ! !INPUT/OUTPUT PARAMETERS:
  !   ik     : k-point number (in,integer)
  !   ispn   : first-variational spin index (in,integer)
  !   apwalm : APW matching coefficients
  !            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
  !   vgpc   : G+k-vectors in Cartesian coordinates
  !   evalfv : first-variational eigenvalues (out,real(nstfv))
  !   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
  ! !DESCRIPTION:
  ! This routine will perform several ARPACK iterations
!
  !BOC
      Implicit None
#ifdef DEBUG
  !include declarations for timing output of ARPACK
#include "./debugf90.h"
#endif
  ! arguments
      Integer, Intent (In) :: ik
      Integer, Intent (In) :: ispn
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Real (8), Intent (Inout) :: evalfv (nstfv, nspnfv)
      Complex (8), Intent (Inout) :: evecfv (nmatmax, nstfv, nspnfv)
!
  ! local variables
      Logical :: packed
      Type (evsystem) :: system
!
      Integer :: n
      Real :: cpu0, cpu1, cpu2
      Complex (8) :: zero, one
      Parameter (zero=(0.0D+0, 0.0D+0), one=(1.0D+0, 0.0D+0))
  !IO vars
      Integer :: koffset, recl
      Character (256) :: outfilenamestring, filetag
      External outfilenamestring
!
  !ARPACK Interface vars
      Integer :: ido, nev, ncv, lworkl, info, infoznaupd, info2, j, i
      Integer :: nevmax, ncvmax, nmax
      Integer :: nconv, maxitr, ishfts, mode, ldv
      Integer :: iparam (11), ipntr (14)
      Complex (8), Allocatable :: resid (:), v (:, :), workev (:), &
     & workl (:), d (:)
      Complex (8), Pointer :: workd (:)
      Real (8), Allocatable :: rwork (:), rd (:)
      Integer, Allocatable :: idx (:)
      Complex (8) :: sigma
      Character :: bmat * 1, which * 2
      Real (8) :: tol
      Logical :: rvec
      Logical :: select (nmat(ispn, ik))
      Complex (8), Pointer :: vin (:), vout (:)
!
#ifdef DEBUG
      ndigit = - 3
      logfil = 6
      mngets = 1
      mnaitr = 1
      mnapps = 1
      mnaupd = 1
      mnaup2 = 1
      mneigh = 1
      mneupd = 1
      Open (logfil, File="ARPACK.OUT", Action="WRITE")
#endif
!
!
  !##################
  !ARPACK parameters
  !##################
      nev = nstfv
      ncv = 2 * nev
      ncv = Min (2*nev, maxncv)
      ncv = Max (ncv, nev+2)
      nevmax = nev
      ncvmax = ncv
      nmax = nmatmax
      n = nmat (ispn, ik)
      ldv = n
      lworkl = 3 * ncvmax * ncvmax + 5 * ncvmax
      Allocate (workd(3*nmax))
      Allocate (resid(nmax))
      Allocate (v(ldv, ncvmax))
      Allocate (workev(2*ncvmax))
      Allocate (workl(lworkl))
      Allocate (d(ncvmax))
      Allocate (rwork(ncvmax))
      Allocate (rd(ncvmax), idx(ncvmax))
      bmat = 'G'
      which = 'LM'
      If (lowesteval .Eq.-1.d0) Then
         Call minenergy (sigma)
      Else
         sigma = dcmplx (lowesteval, 0)
      End If
!
      resid (:) = 0.0
      tol = input%groundstate%solver%epsarpack
      ido = 0
      info = 0
      ishfts = 1
      maxitr = 40 * nstfv
      mode = 3
      iparam (1) = ishfts
      iparam (3) = maxitr
      iparam (7) = mode
  !################################
  !open file with previous residual
  !################################
      Inquire (IoLength=Recl) resid
      koffset = ik - firstk (procofk(ik)) + 1
      infoznaupd = 0
!
  !##################
  !setup hamiltonian#
  !##################
!
!
      If (associated(input%groundstate%solver)) Then
         packed = input%groundstate%solver%packedmatrixstorage
      Else
         packed = .True.
      End If
!
      Call newsystem (system, packed, n)
      Call hamiltonandoverlapsetup (system, ngk(ispn, ik), apwalm, &
     & igkig(1, ispn, ik), vgpc)
!
!
      Call cpu_time (cpu0)
  !#######################################################################
  !calculate LU decomposition to be used in the reverse communication loop
  !#######################################################################
!
      Call HermiteanMatrixAXPY (-sigma, system%overlap, &
     & system%hamilton)
      Call HermiteanMatrixLU (system%hamilton)
      Call cpu_time (cpu1)
  !################################################
  !# reverse comunication loop of arpack library: #
  !################################################
!
      Do i = 1, maxitr
         Call znaupd (ido, bmat, n, which, nev, tol, resid, ncv, v, &
        & ldv, iparam, ipntr, workd, workl, lworkl, rwork, infoznaupd)
         vin => workd (ipntr(1) :ipntr(1)+n-1)
         vout => workd (ipntr(2) :ipntr(2)+n-1)
!
         If (ido .Eq.-1 .Or. ido .Eq. 1) Then
!
            Call Hermiteanmatrixvector (system%overlap, one, vin, zero, &
           & vout)
            Call Hermiteanmatrixlinsolve (system%hamilton, vout)
         Else If (ido .Eq. 1) Then
            Call zcopy (n, workd(ipntr(3)), 1, vout, 1)
            Call Hermiteanmatrixlinsolve (system%hamilton, vout)
         Else If (ido .Eq. 2) Then
            Call Hermiteanmatrixvector (system%overlap, one, vin, zero, &
           & vout)
!
         Else
            Exit
         End If
      End Do
  !###############
  ! errorhandling
  !###############
      If (infoznaupd .Ne. 0) Then
         Print *, ' '
         Print *, ' Error with znaupd, info = ', infoznaupd
         Print *, ' Check the documentation of znaupd'
         Print *, ' '
         Stop
      Else
!
         If (i .Gt. maxitr) Then
            Write (*,*) "Error reached maximum iteration count in arpac&
           &k."
            Stop
         End If
     !########################
     !post processing of evec
     !########################
         rvec = .True.
         select = .True.
         Call zneupd (rvec, 'A', select, d, v, n, sigma, workev, bmat, &
        & n, which, nev, tol, resid, ncv, v, n, iparam, ipntr, workd, &
        & workl, lworkl, rwork, info2)
         If (info2 .Ne. 0) Then
            Print *, ' '
            Print *, ' Error with zneupd, info = ', info2
            Print *, ' Check the documentation of zneupd'
            Print *, ' '
            Write (*,*) "eval", d (1:nev)
            Write (*,*) "iter", i
            Stop
         End If
!
      End If
      Call cpu_time (cpu2)
      timefv = timefv + cpu2 - cpu0
#ifdef DEBUG
      Close (logfil)
#endif
      If (rank .Eq. 0) write (60,*) "k=", ik, "ARPACK iterations", i
      If (rank .Eq. 0) write (60,*) "matrixsize", n, "time LU", cpu1 - &
     & cpu0, "iterations", cpu2 - cpu1
      If (rank .Eq. 0) write (60,*) "minenergy (inversioncenter)", dble &
     & (sigma)
!
  !##########################
  !sort and copy eigenvectors
  !##########################
      rd = real (d)
      Call sortidx (nstfv, rd(:), idx(:))
      Do j = 1, nstfv
         evecfv (1:n, j, ispn) = v (1:n, idx(j))
         evalfv (j, ispn) = rd (idx(j))
      End Do
      Call deleteystem (system)
      Deallocate (workd, resid, v, workev, workl, d)
      Deallocate (rwork, rd, idx)
      Return
End Subroutine iterativearpacksecequn
!EOC
