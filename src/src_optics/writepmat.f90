
!BOP
!!ROUTINE: writepmat
!!INTERFACE:
!
Subroutine writepmat
!
!!USES:
      Use modinput
      Use modmain, Only: nkpt, ngkmax, apwordmax, lmmaxapw, natmtot, &
      &  nmatmax, nstfv, nstsv, nlotot, nlomax, lolmax, task, vkl, vgkl, &
      &  ngk, gkc, tpgkc, sfacgk, igkig, vgkc, filext
      Use modmpi
      Use modxs

!!DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMAT.OUT}, {\tt PMAT\_XS.OUT} or
!   {\tt PMAT\_SCR.OUT} depending on the context of the execution.

!!REVISION HISTORY:
!   Created 2006 (S. Sagmeister)
!   Modifications, August 2010 (S. Sagmeister)
!   Readjusted and re-parallelized, April 2013 (DIN)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ik, recl
      Complex(8), Allocatable :: apwalmt (:, :, :, :)
      Complex(8), Allocatable :: evecfvt (:, :)
      Complex(8), Allocatable :: evecsvt (:, :)
      Complex(8), Allocatable :: pmat (:, :, :, :)
      logical :: fast
      
      integer :: n, iproc, ikloc, nkloc
      integer, allocatable :: ikmap(:,:)

      if (.not.associated(input%properties%momentummatrix)) &
      &  input%properties%momentummatrix => getstructmomentummatrix(emptynode)
      fast=input%properties%momentummatrix%fastpmat

! initialise universal variables
      Call init0
      Call init1

! generating mapping arrays
      allocate(ikmap(nkpt,2))
      ikmap(:,:) = 0
      do iproc = 0, procs-1
        ! k-point interval for all processes
        kpari = firstofset(iproc, nkpt)
        kparf = lastofset(iproc, nkpt)
        ikloc = 0
        do ik = kpari, kparf
          ikloc = ikloc+1
          ikmap(ik,1) = ikloc
          ikmap(ik,2) = iproc
        end do ! ik
      end do ! iproc

! k-point interval for a given process
      kpari = firstofset(rank, nkpt)
      kparf = lastofset(rank, nkpt)
      nkloc = kparf-kpari+1

! Debugging info      
!      write(*,'(">>>> process ", i4, " takes ", i4," k-points in range [", 2i8, "]")') rank, nkloc, kpari, kparf
!      if (rank==0) then
!        do ik = 1, nkpt
!          write(*,*) 'k-point ', ik, ' mapk=', mapk(ik,:)
!        end do
!      end if
      
      Allocate (apwalmt(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (evecfvt(nmatmax, nstfv))
      Allocate (evecsvt(nstsv, nstsv))
! read in the density and potentials from file
      call readstate
! find the new linearisation energies
      call linengy
! generate the APW radial functions
      call genapwfr
! generate the local-orbital radial functions
      call genlofr

      if (fast) then
         If (allocated(apwcmt)) deallocate(apwcmt)
         Allocate(apwcmt(nstfv,apwordmax,lmmaxapw,natmtot))
         If (allocated(ripaa)) deallocate(ripaa)
         Allocate(ripaa(apwordmax,lmmaxapw,apwordmax,lmmaxapw,natmtot,3))
         If (nlotot .Gt. 0) Then
            If (allocated(locmt)) deallocate(locmt)
            Allocate(locmt(nstfv,nlomax,-lolmax:lolmax,natmtot))
            If (allocated(ripalo)) deallocate(ripalo)
            Allocate(ripalo(apwordmax,lmmaxapw,nlomax,-lolmax:lolmax,natmtot,3))
            If (allocated(riploa)) deallocate (riploa)
            Allocate(riploa(nlomax,-lolmax:lolmax,apwordmax,lmmaxapw,natmtot,3))
            If (allocated(riplolo)) deallocate(riplolo)
            Allocate(riplolo(nlomax,-lolmax:lolmax,nlomax,-lolmax:lolmax,natmtot,3))
         End If
! calculate gradient of radial functions times spherical harmonics
         Call pmatrad
      End If

!---------------------------------
! Main loop over k-points
!---------------------------------

      allocate(pmat(3,nstsv,nstsv,nkloc))
      ikloc = 0
      Do ik = kpari, kparf
         ikloc = ikloc+1
! get the eigenvectors and values from file
         Call getevecfv (vkl(1, ik), vgkl(1, 1, 1, ik), evecfvt)
         Call getevecsv (vkl(1, ik), evecsvt)
! find the matching coefficients
         Call match (ngk(1, ik), gkc(1, 1, ik), tpgkc(1, 1, 1, ik), &
        &  sfacgk(1, 1, 1, ik), apwalmt)
         If (fast) Then
            ! generate APW expansion coefficients for muffin-tin
            Call genapwcmt(input%groundstate%lmaxapw, ngk(1, ik), 1, &
           &  nstfv, apwalmt, evecfvt, apwcmt)
            ! generate local orbital expansion coefficients for muffin-tin
            If (nlotot.Gt.0) Call genlocmt(ngk(1,ik), 1, nstfv, evecfvt, locmt)
            ! calculate the momentum matrix elements
            Call genpmatxs(ngk(1,ik), igkig(1,1,ik), vgkc(1,1,1,ik), evecfvt, evecsvt, pmat(:,:,:,ikloc))
         Else
        ! calculate the momentum matrix elements
            Call genpmat(ngk(1,ik), igkig(1,1,ik), vgkc(1,1,1,ik), &
           &  apwalmt, evecfvt, evecsvt, pmat(:,:,:,ikloc))
         End If
      End Do ! ik
      call barrier

! find the record length
      inquire(iolength=recl) pmat(:,:,:,1)
      open(50,file='PMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
     &  status='REPLACE',recl=recl)
      do ik = 1, nkpt
        ikloc = ikmap(ik,1)
        iproc = ikmap(ik,2)
        if (rank==iproc) then
          write(50,rec=ik) pmat(:,:,:,ikloc)
        end if
        call barrier
      end do
      close(50)
!
      deallocate(apwalmt, evecfvt, evecsvt, pmat, ikmap)
      if (fast) then
         deallocate(apwcmt)
         deallocate(ripaa)
         if (nlotot > 0) then
            deallocate(locmt)
            deallocate(ripalo,riploa,riplolo)
         end if
      End If

      if (rank==0) then
          Write (*,*)
          Write (*, '("Info(writepmat):")')
          Write (*, '("  Momentum matrix elements written to file PMAT.OUT")')
          Write (*,*)
      end if

End Subroutine writepmat
!EOC
