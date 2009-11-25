!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: writepmatxs
! !INTERFACE:
!
!
Subroutine writepmatxs
! !USES:
      Use modinput
      Use modmain, Only: nkpt, ngkmax, apwordmax, lmmaxapw, natmtot, &
     & nmatmax, nstfv, nstsv, nlotot, nlomax, lolmax, task, vkl, vgkl, &
     & ngk, gkc, tpgkc, sfacgk, igkig, vgkc
      Use modmpi
      Use modxs
      Use m_putpmat
      Use m_genfilname
! !DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMAT\_XS.OUT}. Derived from
!   the routine {\tt writepmat}.
!
! !REVISION HISTORY:
!   Created 2006 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! local variables
      Character (*), Parameter :: thisnam = 'writepmatxs'
      Integer :: ik
      Character (32) :: fnam
      Complex (8), Allocatable :: apwalmt (:, :, :, :)
      Complex (8), Allocatable :: evecfvt (:, :)
      Complex (8), Allocatable :: evecsvt (:, :)
      Complex (8), Allocatable :: pmat (:, :, :)
      If (tscreen) Then
         fnam = 'PMAT'
         Call genfilname (basename=trim(fnam), appfilext=.True., &
        & filnam=fnpmat)
         Call genfilname (basename=trim(fnam), procs=procs, rank=rank, &
        & appfilext=.True., filnam=fnpmat_t)
      Else
         fnam = 'PMAT_XS'
         Call genfilname (basename=trim(fnam), filnam=fnpmat)
         Call genfilname (basename=trim(fnam), procs=procs, rank=rank, &
        & filnam=fnpmat_t)
      End If
  ! initialise universal variables
      Call init0
      Call init1
      Call init2
  ! generate index ranges for parallel execution
      Call genparidxran ('k', nkpt)
  ! k-point interval for process
      kpari = firstofset (rank, nkpt)
      kparf = lastofset (rank, nkpt)
      Allocate (apwalmt(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfvt(nmatmax, nstfv))
      Allocate (evecsvt(nstsv, nstsv))
  ! allocate the momentum matrix elements array
      Allocate (pmat(3, nstsv, nstsv))
  ! get eigenvectors for q=0
      If ( .Not. tscreen) Call genfilname (iqmt=0, setfilext=.True.)
  ! generate band combinations
      Call ematbdcmbs (1)
      If (input%xs%fastpmat) Then
         If (allocated(apwcmt)) deallocate (apwcmt)
         Allocate (apwcmt(nstsv, apwordmax, lmmaxapw, natmtot))
         If (allocated(ripaa)) deallocate (ripaa)
         Allocate (ripaa(apwordmax, lmmaxapw, apwordmax, lmmaxapw, &
        & natmtot, 3))
         If (nlotot .Gt. 0) Then
            If (allocated(locmt)) deallocate (locmt)
            Allocate (locmt(nstsv, nlomax,-lolmax:lolmax, natmtot))
            If (allocated(ripalo)) deallocate (ripalo)
            Allocate (ripalo(apwordmax, lmmaxapw, &
           & nlomax,-lolmax:lolmax, natmtot, 3))
            If (allocated(riploa)) deallocate (riploa)
            Allocate (riploa(nlomax,-lolmax:lolmax, apwordmax, &
           & lmmaxapw, natmtot, 3))
            If (allocated(riplolo)) deallocate (riplolo)
            Allocate (riplolo(nlomax,-lolmax:lolmax, &
           & nlomax,-lolmax:lolmax, natmtot, 3))
         End If
     ! calculate gradient of radial functions times spherical harmonics
         Call pmatrad
      End If
      Do ik = kpari, kparf
         Call chkpt (2, (/ task, ik /), 'ematqk: task, k - point index;&
        & momentum matrix elements')
     ! get the eigenvectors and values from file
         Call getevecfv (vkl(1, ik), vgkl(1, 1, 1, ik), evecfvt)
         Call getevecsv (vkl(1, ik), evecsvt)
     ! find the matching coefficients
         Call match (ngk(1, ik), gkc(1, 1, ik), tpgkc(1, 1, 1, ik), &
        & sfacgk(1, 1, 1, ik), apwalmt)
         If (input%xs%fastpmat) Then
        ! generate APW expansion coefficients for muffin-tin
            Call genapwcmt (input%groundstate%lmaxapw, ngk(1, ik), 1, &
           & nstfv, apwalmt, evecfvt, apwcmt)
        ! generate local orbital expansion coefficients for muffin-tin
            If (nlotot .Gt. 0) Call genlocmt (ngk(1, ik), 1, nstfv, &
           & evecfvt, locmt)
        ! calculate the momentum matrix elements
            Call genpmat2 (ngk(1, ik), igkig(1, 1, ik), vgkc(1, 1, 1, &
           & ik), evecfvt, evecsvt, pmat)
         Else
        ! calculate the momentum matrix elements
            Call genpmat (ngk(1, ik), igkig(1, 1, ik), vgkc(1, 1, 1, &
           & ik), apwalmt, evecfvt, evecsvt, pmat)
         End If
     ! parallel write
         Call putpmat (ik, .True., trim(fnpmat), pmat)
      End Do
      Call barrier
      Deallocate (apwalmt, evecfvt, evecsvt, pmat)
      If (input%xs%fastpmat) Then
         Deallocate (apwcmt)
         If (nlotot .Gt. 0) Then
            Deallocate (locmt)
            Deallocate (ripaa, ripalo, riploa, riplolo)
         End If
      End If
      Call barrier
      Write (unitout, '(a)') "Info(" // trim (thisnam) // "): momentum &
     &matrix elements finished"
  ! reset global file extension to default
      Call genfilname (setfilext=.True.)
End Subroutine writepmatxs
!EOC
