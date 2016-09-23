
!BOP
!!ROUTINE: writepmat
!!INTERFACE:
!
Subroutine writepmat
!
!!USES:
      Use modinput
      Use modmain
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
      
      integer :: i, j, n, iproc, ikloc, nkloc, isym, lspl, ilspl, iv(3)
      integer, allocatable :: ikpmap(:,:)
      real(8) :: t1, v(3), v1(3), v2(3), pm(9), sl(3,3), sc(3,3)
      complex(8) :: p(3), o(6)
! external function      
      real(8) :: r3taxi

      if (.not.associated(input%properties%momentummatrix)) &
      &  input%properties%momentummatrix => getstructmomentummatrix(emptynode)
      fast=input%properties%momentummatrix%fastpmat

! initialise universal variables
      Call init0
      Call init1

! generating mapping arrays
      allocate(ikpmap(nkpt,2))
      ikpmap(:,:) = 0
      do iproc = 0, procs-1
        ! k-point interval for all processes
        kpari = firstofset(iproc, nkpt)
        kparf = lastofset(iproc, nkpt)
        ikloc = 0
        do ik = kpari, kparf
          ikloc = ikloc+1
          ikpmap(ik,1) = ikloc
          ikpmap(ik,2) = iproc
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
         Call getevecfv(vkl(1,ik),vgkl(1,1,1,ik),evecfvt)
         !call getsymevecfv(vkl(1,ik),vgkl(1,1,1,ik),evecfvt)
         Call getevecsv(vkl(1,ik),evecsvt)
! find the matching coefficients
         Call match (ngk(1, ik), gkc(1, 1, ik), tpgkc(1, 1, 1, ik), &
         &           sfacgk(1, 1, 1, ik), apwalmt)
         If (fast) Then
            ! generate APW expansion coefficients for muffin-tin
            Call genapwcmt(input%groundstate%lmaxapw, ngk(1, ik), 1, &
            &              nstfv, apwalmt, evecfvt, apwcmt)
            ! generate local orbital expansion coefficients for muffin-tin
            If (nlotot.Gt.0) Call genlocmt(ngk(1,ik), 1, nstfv, evecfvt, locmt)
            ! calculate the momentum matrix elements
            Call genpmatxs(ngk(1,ik), igkig(1,1,ik), vgkc(1,1,1,ik), evecfvt, evecsvt, pmat(:,:,:,ikloc))
         Else
        ! calculate the momentum matrix elements
            Call genpmat(ngk(1,ik), igkig(1,1,ik), vgkc(1,1,1,ik), &
           &  apwalmt, evecfvt, evecsvt, pmat(:,:,:,ikloc))
         End If

         ! symmetrize |pmat|^2
         do i = 1, nstsv
         do j = 1, nstsv
           o(:) = zzero
           do isym = 1, nsymcrys
             lspl = lsplsymc(isym)
             sc(:,:) = symlatc(:,:,lspl)
             call r3mv(sc,dble(pmat(:,i,j,ikloc)),v1)
             call r3mv(sc,aimag(pmat(:,i,j,ikloc)),v2)
             p(:) = cmplx(v1(:),v2(:),8)
             o(1) = o(1)+p(1)*conjg(p(1))
             o(2) = o(2)+p(2)*conjg(p(2))
             o(3) = o(3)+p(3)*conjg(p(3))
             if (symlatd(lspl)<0) then
               o(4) = o(4)+conjg(p(2))*p(1)
               o(5) = o(5)+conjg(p(3))*p(1)
               o(6) = o(6)+conjg(p(3))*p(2)
             else
               o(4) = o(4)+p(2)*conjg(p(1))
               o(5) = o(5)+p(3)*conjg(p(1))
               o(6) = o(6)+p(3)*conjg(p(2))
             end if
           end do ! isym
           pm(1) = dble(o(1))/dble(nsymcrys)
           pm(2) = dble(o(2))/dble(nsymcrys)
           pm(3) = dble(o(3))/dble(nsymcrys)
           pm(4) = dble(o(4))/dble(nsymcrys)
           pm(5) = dble(o(5))/dble(nsymcrys)
           pm(6) = dble(o(6))/dble(nsymcrys)
           pm(7) = aimag(o(4))/dble(nsymcrys)
           pm(8) = aimag(o(5))/dble(nsymcrys)
           pm(9) = aimag(o(6))/dble(nsymcrys)
           !if (rank==0) call write_pmatij(i,j,pm(1:3))
         end do
         end do
         
      End Do ! ik
      call barrier

! find the record length
      inquire(iolength=recl) pmat(:,:,:,1)
      open(50,file='PMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
     &  status='REPLACE',recl=recl)
      do ik = 1, nkpt
        ikloc = ikpmap(ik,1)
        iproc = ikpmap(ik,2)
        if (rank==iproc) then
          write(50,rec=ik) pmat(:,:,:,ikloc)
        end if
        call barrier
      end do
      close(50)
!
      deallocate(apwalmt, evecfvt, evecsvt, pmat, ikpmap)
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

contains

    subroutine write_pmatij(i,j,p)
      integer, intent(in) :: i, j
      real(8), intent(in) :: p(3)
      write(77,*) "band ", i, " - ", j
      write(77,'(a,3f18.6)') "|p_aa|^2 =", p(1), p(2), p(3)
    end subroutine      
      
End Subroutine writepmat
!EOC
