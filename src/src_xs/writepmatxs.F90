
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writepmatxs
! !INTERFACE:
subroutine writepmatxs(lgather)
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
  ! arguments
  logical :: lgather
  ! local variables
  character(*), parameter :: thisnam='writepmatxs'
  integer :: ik
  character(32) :: fnam
  complex(8), allocatable :: apwalmt(:,:,:,:)
  complex(8), allocatable :: apwcmt(:,:,:,:)
  complex(8), allocatable :: locmt(:,:,:,:)
  complex(8), allocatable :: evecfvt(:,:)
  complex(8), allocatable :: evecsvt(:,:)
  complex(8), allocatable :: pmat(:,:,:)
  real(8), allocatable :: ripaa(:,:,:,:,:,:)
  real(8), allocatable :: ripalo(:,:,:,:,:,:)
  real(8), allocatable :: riploa(:,:,:,:,:,:)
  real(8), allocatable :: riplolo(:,:,:,:,:,:)
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
  ! analytic evaluation of interstitial contribution
  if (pmatira) then
     write(unitout,'(a)') 'Info('//thisnam//'): using an analytic method for &
          &calculation of momentum'
     write(unitout,'(a)') ' matrix elements in the interstitial'
  end if
  ! initialise universal variables
  if (calledxs.eq.1) call init0
  call init1
  call init2xs
  ! generate index ranges for parallel execution
  call genparidxran('k')
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
  if (lgather) goto 10
  if (pmatstrat.ne.0) then
     allocate(ripaa(apwordmax,lmmaxapw,apwordmax,lmmaxapw,natmtot,3))
     allocate(apwcmt(nstsv,apwordmax,lmmaxapw,natmtot))
     if (nlotot.gt.0) then
        allocate(ripalo(apwordmax,lmmaxapw,nlomax,-lolmax:lolmax,natmtot,3))
        allocate(riploa(nlomax,-lolmax:lolmax,apwordmax,lmmaxapw,natmtot,3))
        allocate(riplolo(nlomax,-lolmax:lolmax,nlomax,-lolmax:lolmax,natmtot,3))
        allocate(locmt(nstsv,nlomax,-lolmax:lolmax,natmtot))
     end if
     ! calculate gradient of radial functions times spherical harmonics
     call pmatrad(ripaa,ripalo,riploa,riplolo)
  end if
  do ik=kpari,kparf
     if ((modulo(ik-kpari+1,max((kparf-kpari+1)/10,1)).eq.0).or.(ik.eq.kparf)) &
          write(*,'("Info(",a,"): ",I6," of ",I6,I6," k-points")') thisnam,ik, &
          kpari,kparf
     ! get the eigenvectors and values from file
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfvt)
     call getevecsv(vkl(1,ik),evecsvt)
     ! calculate the momentum matrix elements
     if (pmatstrat.eq.0) then
        ! find the matching coefficients
        call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1), &
             apwalmt)
        call genpmat(ngk(ik,1),igkig(1,ik,1),vgkc(1,1,ik,1),apwalmt,evecfvt, &
             evecsvt,pmat)
     else
        ! find the matching coefficients
        call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1), &
             apwalmt)
        ! generate APW expansion coefficients for muffin-tin
        call genapwcmt(lmaxapw,ngk(ik,1),1,nstfv,apwalmt,evecfvt,apwcmt)
        ! generate local orbital expansion coefficients for muffin-tin
        if (nlotot.gt.0) call genlocmt(ngk(ik,1),1,nstfv,evecfvt,locmt)
        call genpmat2(ngk(ik,1),igkig(1,ik,1),vgkc(1,1,ik,1),ripaa,ripalo, &
             riploa,riplolo,apwcmt,locmt,evecfvt,evecsvt,pmat)
     end if
     call putpmat(ik,.false.,trim(fnpmat_t),pmat)
     ! synchronize for common number of k-points to all processes
     if (ik-kpari+1 <= nkpt/procs) call barrier
  end do
  call barrier
10 continue
  if ((procs.gt.1).and.(rank.eq.0)) call pmatgather
  deallocate(evecfvt,evecsvt,pmat)
  if (pmatstrat.eq.0) then
     deallocate(apwalmt)
  else
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
