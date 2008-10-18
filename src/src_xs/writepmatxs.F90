
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writepmatxs
! !INTERFACE:
subroutine writepmatxs
! !USES:
  use modmain
  use modmpi
  use modxs
  use m_putpmat
  use m_genfilname
! !DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMAT\_XS.OUT}. Derived from
!   the routine {\tt writepmat}.
!
! !REVISION HISTORY:
!   Created 2006 (Sagmeister)
!EOP
!BOC
  implicit none
  ! local variables
  character(*), parameter :: thisnam='writepmatxs'
  integer :: ik
  character(32) :: fnam
  complex(8), allocatable :: apwalmt(:,:,:,:)
  complex(8), allocatable :: evecfvt(:,:)
  complex(8), allocatable :: evecsvt(:,:)
  complex(8), allocatable :: pmat(:,:,:)
  if (tscreen) then
     fnam='PMAT'
     call genfilname(basename=trim(fnam),appfilext=.true.,filnam=fnpmat)
     call genfilname(basename=trim(fnam),procs=procs,rank=rank, &
          appfilext=.true.,filnam=fnpmat_t)
  else
     fnam='PMAT_XS'
     call genfilname(basename=trim(fnam),filnam=fnpmat)
     call genfilname(basename=trim(fnam),procs=procs,rank=rank,filnam=fnpmat_t)
  end if
  ! initialise universal variables
  call init0
  call init1
  call init2
  ! generate index ranges for parallel execution
  call genparidxran('k',nkpt)
  ! k-point interval for process
  kpari=firstofset(rank,nkpt)
  kparf=lastofset(rank,nkpt)
  allocate(apwalmt(ngkmax,apwordmax,lmmaxapw,natmtot))
  allocate(evecfvt(nmatmax,nstfv))
  allocate(evecsvt(nstsv,nstsv))
  ! allocate the momentum matrix elements array
  allocate(pmat(3,nstsv,nstsv))
  ! get eigenvectors for q=0
  if (.not.tscreen) call genfilname(iqmt=0,setfilext=.true.)
  ! generate band combinations
  call ematbdcmbs(1)
  if (fastpmat) then
     if (allocated(apwcmt)) deallocate(apwcmt)
     allocate(apwcmt(nstsv,apwordmax,lmmaxapw,natmtot))
     if (allocated(ripaa)) deallocate(ripaa)
     allocate(ripaa(apwordmax,lmmaxapw,apwordmax,lmmaxapw,natmtot,3))
     if (nlotot.gt.0) then
        if (allocated(locmt)) deallocate(locmt)
        allocate(locmt(nstsv,nlomax,-lolmax:lolmax,natmtot))
        if (allocated(ripalo)) deallocate(ripalo)
        allocate(ripalo(apwordmax,lmmaxapw,nlomax,-lolmax:lolmax,natmtot,3))
        if (allocated(riploa)) deallocate(riploa)
        allocate(riploa(nlomax,-lolmax:lolmax,apwordmax,lmmaxapw,natmtot,3))
        if (allocated(riplolo)) deallocate(riplolo)
        allocate(riplolo(nlomax,-lolmax:lolmax,nlomax,-lolmax:lolmax,natmtot,3))
     end if
     ! calculate gradient of radial functions times spherical harmonics
     call pmatrad
  end if
  do ik=kpari,kparf
     if ((modulo(ik-kpari+1,max((kparf-kpari+1)/10,1)).eq.0).or.(ik.eq.kparf)) &
          write(*,'("Info(",a,"): ",I6," of ",I6,I6," k-points")') thisnam,ik, &
          kpari,kparf
     ! get the eigenvectors and values from file
     call getevecfv(vkl(1,ik),vgkl(1,1,1,ik),evecfvt)
     call getevecsv(vkl(1,ik),evecsvt)
     ! find the matching coefficients
     call match(ngk(1,ik),gkc(1,1,ik),tpgkc(1,1,1,ik),sfacgk(1,1,1,ik), &
          apwalmt)
     if (fastpmat) then
        ! generate APW expansion coefficients for muffin-tin
        call genapwcmt(lmaxapw,ngk(1,ik),1,nstfv,apwalmt,evecfvt,apwcmt)
        ! generate local orbital expansion coefficients for muffin-tin
        if (nlotot.gt.0) call genlocmt(ngk(1,ik),1,nstfv,evecfvt,locmt)
        ! calculate the momentum matrix elements
        call genpmat2(ngk(1,ik),igkig(1,1,ik),vgkc(1,1,1,ik),evecfvt,evecsvt, &
	     pmat)
     else
        ! calculate the momentum matrix elements
        call genpmat(ngk(1,ik),igkig(1,1,ik),vgkc(1,1,1,ik),apwalmt,evecfvt, &
             evecsvt,pmat)
     end if
     ! parallel write
     call putpmat(ik,.true.,trim(fnpmat),pmat)
  end do
  call barrier
  deallocate(apwalmt,evecfvt,evecsvt,pmat)
  if (fastpmat) then
     deallocate(apwcmt)
     if (nlotot.gt.0) then
        deallocate(locmt)
        deallocate(ripaa,ripalo,riploa,riplolo)
     end if
  end if
  call barrier
  write(unitout,'(a)') "Info("//trim(thisnam)//"): momentum matrix elements &
       &finished"
  ! reset global file extension to default
  call genfilname(setfilext=.true.)
end subroutine writepmatxs
!EOC
