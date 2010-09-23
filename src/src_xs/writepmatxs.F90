
! Copyright (C) 2005-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writepmatxs
! !INTERFACE:
Subroutine writepmatxs
! !USES:
      Use modinput
      Use modmain, Only: nkpt, ngkmax, apwordmax, lmmaxapw, natmtot, &
     & nmatmax, nstfv, nstsv, nlotot, nlomax, lolmax, task, vkl, vgkl, &
     & ngk, gkc, tpgkc, sfacgk, igkig, vgkc, filext
      Use modmpi
      Use modxs
      Use m_putpmat
      Use m_genfilname
! !DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMAT.OUT}, {\tt PMAT\_XS.OUT} or
!   {\tt PMAT\_SCR.OUT} depending on the context of the execution.
!
! !REVISION HISTORY:
!   Created 2006 (S. Sagmeister)
!   Modifications, August 2010 (S. Sagmeister)
!EOP
!BOC
      Implicit None
  ! local variables
      Integer :: ik,recl
      Character (32) :: fnam
      logical :: fast
      Complex (8), Allocatable :: apwalmt (:, :, :, :)
      Complex (8), Allocatable :: evecfvt (:, :)
      Complex (8), Allocatable :: evecsvt (:, :)
      Complex (8), Allocatable :: pmat (:, :, :)
      Logical, External :: tqgamma
      fast=(task.ne.120).or.((task.eq.120).and.input%properties%momentummatrix%fastpmat)
write(500,*) 'sag: 1'; call flushifc(500)
      if (task .ne. 120) then
        if (.not.tqgamma(1)) return
      end if
write(500,*) 'sag: 2'; call flushifc(500)
      tscreen=(task .Ge. 400) .And. (task .Le. 499)
      If ((task .eq. 120).or. tscreen) Then
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
write(500,*) 'sag: 3'; call flushifc(500)
  ! initialise universal variables
      Call init0
      Call init1
write(500,*) 'sag: 4'; call flushifc(500)
      if (task .ne. 120) Call init2
  ! generate index ranges for parallel execution
      Call genparidxran ('k', nkpt)
write(500,*) 'sag: 5'; call flushifc(500)
  ! k-point interval for process
      kpari = firstofset (rank, nkpt)
      kparf = lastofset (rank, nkpt)
      Allocate (apwalmt(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfvt(nmatmax, nstfv))
      Allocate (evecsvt(nstsv, nstsv))
  ! allocate the momentum matrix elements array
      Allocate (pmat(3, nstsv, nstsv))
write(500,*) 'sag: 6'; call flushifc(500)
      if (task .eq. 120) then
  ! read in the density and potentials from file
        call readstate
  ! find the new linearisation energies
        call linengy
  ! generate the APW radial functions
        call genapwfr
  ! generate the local-orbital radial functions
        call genlofr
  ! find the record length
        inquire(iolength=recl) pmat
        open(50,file='PMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
         status='REPLACE',recl=recl)
      end if
write(500,*) 'sag: 7'; call flushifc(500)
  ! get eigenvectors for q=0
      If ((.Not. tscreen) .and. (task .ne. 120)) Call genfilname (iqmt=0, setfilext=.True.)
  ! generate band combinations
write(500,*) 'sag: 8'; call flushifc(500)
      if (task .eq. 120) then
        call ematbdcmbs(0)
      else
        Call ematbdcmbs(1)
      end if
write(500,*) 'sag: 9'; call flushifc(500)
      if (fast) then
         If (allocated(apwcmt)) deallocate (apwcmt)
         Allocate (apwcmt(nstsv, apwordmax, lmmaxapw, natmtot))
write(500,*) 'sag: 9.1'; call flushifc(500)
         If (allocated(ripaa)) deallocate (ripaa)
         Allocate (ripaa(apwordmax, lmmaxapw, apwordmax, lmmaxapw,natmtot, 3))
write(500,*) 'sag: 9.2'; call flushifc(500)
         If (nlotot .Gt. 0) Then
            If (allocated(locmt)) deallocate (locmt)
            Allocate (locmt(nstsv, nlomax,-lolmax:lolmax, natmtot))
write(500,*) 'sag: 9.3'; call flushifc(500)
            If (allocated(ripalo)) deallocate (ripalo)
            Allocate (ripalo(apwordmax, lmmaxapw, nlomax,-lolmax:lolmax, natmtot, 3))
write(500,*) 'sag: 9.4'; call flushifc(500)
            If (allocated(riploa)) deallocate (riploa)
            Allocate (riploa(nlomax,-lolmax:lolmax, apwordmax, lmmaxapw, natmtot, 3))
write(500,*) 'sag: 9.5'; call flushifc(500)
            If (allocated(riplolo)) deallocate (riplolo)
            Allocate (riplolo(nlomax,-lolmax:lolmax, nlomax,-lolmax:lolmax, natmtot, 3))
write(500,*) 'sag: 9.6'; call flushifc(500)
         End If
     ! calculate gradient of radial functions times spherical harmonics
         Call pmatrad
write(500,*) 'sag: 9.7'; call flushifc(500)
      End If
write(500,*) 'sag: 10'; call flushifc(500)
      Do ik = kpari, kparf
write(500,*) 'sag: 10.1'; call flushifc(500)
         if (task .ne. 120) Call chkpt (2, (/ task, ik /), 'ematqk: task, k - point index;&
        & momentum matrix elements')
     ! get the eigenvectors and values from file
         Call getevecfv (vkl(1, ik), vgkl(1, 1, 1, ik), evecfvt)
write(500,*) 'sag: 10.2'; call flushifc(500)
         Call getevecsv (vkl(1, ik), evecsvt)
write(500,*) 'sag: 10.3'; call flushifc(500)
     ! find the matching coefficients
         Call match (ngk(1, ik), gkc(1, 1, ik), tpgkc(1, 1, 1, ik), &
        & sfacgk(1, 1, 1, ik), apwalmt)
         If (fast) Then
        ! generate APW expansion coefficients for muffin-tin
            Call genapwcmt (input%groundstate%lmaxapw, ngk(1, ik), 1, &
           & nstfv, apwalmt, evecfvt, apwcmt)
        ! generate local orbital expansion coefficients for muffin-tin
            If (nlotot .Gt. 0) Call genlocmt (ngk(1, ik), 1, nstfv, &
           & evecfvt, locmt)
        ! calculate the momentum matrix elements
            Call genpmatxs (ngk(1, ik), igkig(1, 1, ik), vgkc(1, 1, 1, &
           & ik), evecfvt, evecsvt, pmat)
         Else
        ! calculate the momentum matrix elements
            Call genpmat (ngk(1, ik), igkig(1, 1, ik), vgkc(1, 1, 1, &
           & ik), apwalmt, evecfvt, evecsvt, pmat)
         End If
         if (task .eq. 120) then
     ! write the matrix elements to direct-access file
           write(50,rec=ik) pmat
         else
     ! parallel write
           Call putpmat (ik, .True., trim(fnpmat), pmat)
         end if
      End Do
write(500,*) 'sag: 11'; call flushifc(500)
      Call barrier
      Deallocate (apwalmt, evecfvt, evecsvt, pmat)
      If (fast) Then
         Deallocate (apwcmt)
         Deallocate (ripaa)
         If (nlotot .Gt. 0) Then
            Deallocate (locmt)
            Deallocate (ripalo, riploa, riplolo)
         End If
      End If
      Call barrier
write(500,*) 'sag: 12'; call flushifc(500)
      if (task .eq. 120) then
        close(50)
        Write (*,*)
        Write (*, '("Info(writepmatxs):")')
        Write (*, '(" momentum matrix elements written to file PMAT.OUT")')
        Write (*,*)
      else
        Write (unitout, '(a)') "Info(writepmatxs): momentum matrix elements finished"
      end if
  ! reset global file extension to default
      Call genfilname (setfilext=.True.)
write(500,*) 'sag: 13'; call flushifc(500)
End Subroutine writepmatxs
!EOC
