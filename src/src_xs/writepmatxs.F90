
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
      Complex (8), Allocatable :: apwalmt (:, :, :, :)
      Complex (8), Allocatable :: evecfvt (:, :)
      Complex (8), Allocatable :: evecsvt (:, :)
      Complex (8), Allocatable :: pmat (:, :, :)
      Character(256) :: string
      logical :: fast
      Logical, External :: tqgamma
      ! check if fast (default) version of matrix elements is used
      fast=.false.
      if (associated(input%properties)) then
        if (associated(input%properties%momentummatrix)) then
	  if (input%properties%momentummatrix%fastpmat) fast=.true.
	end if
      end if	
      fast=(task.ne.120).or.((task.eq.120).and.fast)
      ! check if Q-point is Gamma point
      if (task .ne. 120) then
        if (.not.tqgamma(1)) return
      end if
      ! check if this routine is called for the screening
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
  ! initialise universal variables
      Call init0
      Call init1
      if (task .ne. 120) Call init2
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
      if (task .eq. 120) then
! read density and potentials from file
        If (hybridhf) Then
! in case of HF hybrids use PBE potential
            string=filext
            filext='_PBE.OUT'
            Call readstate
            filext=string
        Else
           Call readstate
        End If
  ! find the new linearisation energies
        call linengy
  ! generate the APW radial functions
        call genapwfr
  ! generate the local-orbital radial functions
        call genlofr
! update potential in case if HF Hybrids
        If (hybridhf) Call readstate
  ! find the record length
        inquire(iolength=recl) pmat
        open(50,file='PMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
         status='REPLACE',recl=recl)
      end if
  ! get eigenvectors for q=0
      If ((.Not. tscreen) .and. (task .ne. 120)) Call genfilname (iqmt=0, setfilext=.True.)
  ! generate band combinations
      if (task .eq. 120) then
        call ematbdcmbs(0)
      else
        Call ematbdcmbs(1)
      end if
      if (fast) then
         If (allocated(apwcmt)) deallocate (apwcmt)
         Allocate (apwcmt(nstfv, apwordmax, lmmaxapw, natmtot))
         If (allocated(ripaa)) deallocate (ripaa)
         Allocate (ripaa(apwordmax, lmmaxapw, apwordmax, lmmaxapw,natmtot, 3))
         If (nlotot .Gt. 0) Then
            If (allocated(locmt)) deallocate (locmt)
            Allocate (locmt(nstfv, nlomax,-lolmax:lolmax, natmtot))
            If (allocated(ripalo)) deallocate (ripalo)
            Allocate (ripalo(apwordmax, lmmaxapw, nlomax,-lolmax:lolmax, natmtot, 3))
            If (allocated(riploa)) deallocate (riploa)
            Allocate (riploa(nlomax,-lolmax:lolmax, apwordmax, lmmaxapw, natmtot, 3))
            If (allocated(riplolo)) deallocate (riplolo)
            Allocate (riplolo(nlomax,-lolmax:lolmax, nlomax,-lolmax:lolmax, natmtot, 3))
         End If
     ! calculate gradient of radial functions times spherical harmonics
         Call pmatrad
      End If
      Do ik = kpari, kparf
         if (task .ne. 120) Call chkpt (2, (/ task, ik /), 'ematqk: task, k - point index;&
        & momentum matrix elements')
     ! get the eigenvectors and values from file
         Call getevecfv (vkl(1, ik), vgkl(1, 1, 1, ik), evecfvt)
         Call getevecsv (vkl(1, ik), evecsvt)
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
        !    write(*,*)"putpmat ik ",ik,"done"
         end if
      End Do
      Call barrier

      Inquire (IoLength=Recl) vkl (:, ik), nstsv, pmat
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
      !  write(*,*)"cp pmat"
      if (.not. input%sharedfs) call  cpFileToNodes( trim(fnpmat))

      if (task .eq. 120) then
        close(50)
        if (rank==0) then
          Write (*,*)
          Write (*, '("Info(writepmatxs):")')
          Write (*, '(" momentum matrix elements written to file PMAT.OUT")')
          Write (*,*)
        end if
      else
        Write (unitout, '(a)') "Info(writepmatxs): momentum matrix elements finished"
      end if
  ! reset global file extension to default
      Call genfilname (setfilext=.True.)
End Subroutine writepmatxs
!EOC
