
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writepmattd
! !INTERFACE:
subroutine writepmattd(lgather)
  ! !USES:
  use modmain
  use modmpi
  use modtddft
  use m_putpmat
  use m_getunit
  use m_filedel
  use m_genfilname
! !DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMAT_TD.OUT}. Derived from
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
  character(*), parameter :: thisnam = 'writepmattd'
  integer recl,ik,iproc,ki,kf,ngridk_save(3)
  integer :: nempty_save
  integer :: j,ist
  integer :: ipar,fu,un
  real(8) :: gmaxvrt
  real(8), parameter :: gmaxvrmax=15.d0
  complex(8), allocatable :: apwalmt(:,:,:,:)
  complex(8), allocatable :: evecfvt(:,:)
  complex(8), allocatable :: evecsvt(:,:)
  complex(8), allocatable :: pmat(:,:,:)
  logical :: existent

  if (pmatira) then
     write(unitout,'(a)') 'Info('//thisnam//'): using an analytic method for &
          &calculation of momentum'
     write(unitout,'(a)') ' matrix elements in the interstitial'
  end if

  ! initialise universal variables
  call init0
  call init1
  call init2td
  ! k-point interval for process
  kpari=firstofset(rank,nkpt)
  kparf=lastofset(rank,nkpt)
!!$  ! jump into next checkpoint
!!$  if (tresume) resumechkpts(1,1)=resumechkpts(1,1)+1
!!$  resumechkpts(1,2)=kpari
!!$  resumechkpts(1,3)=kparf

  allocate(apwalmt(ngkmax,apwordmax,lmmaxapw,natmtot))
  allocate(evecfvt(nmatmax,nstfv))
  allocate(evecsvt(nstsv,nstsv))
  ! allocate the momentum matrix elements array
  allocate(pmat(3,nstsv,nstsv))
  ! get eigenvectors for q=0
  call genfilname(iq=0,setfilext=.true.)

  ! limits for k-point loop
  ki=kpari
  kf=kparf
  ! resume task, first checkpoint is k-point index
  if (lgather) goto 10
!!$  if (tresume) ki=resumechkpts(1,1)
  do ik=ki,kf
     ! get the eigenvectors and values from file
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfvt)
     call getevecsv(vkl(1,ik),evecsvt)
     ! find the matching coefficients
     call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalmt)
     ! calculate the momentum matrix elements
     call genpmat(ngk(ik,1),igkig(1,ik,1),vgkc(1,1,ik,1),apwalmt,evecfvt, &
          evecsvt,pmat)
     call putpmat(ik,.false.,trim(fnpmat_t),pmat)
!!$     resumechkpts(1,1)=ik
!!$     call resupd(un,task,resumechkpts,' : k-point index')
#ifdef MPI     
       if (ik-ki+1 <= nkpt/procs) then
          ! synchronize for common number of k-points to all processes
          call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')
       end if
#endif
  end do
  call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')

10 continue

  ! lgather from processes
  if ((procs.gt.1).and.(rank.eq.0)) call pmatgather()

  ! reset global file extension to default
  call genfilname(setfilext=.true.)
  deallocate(apwalmt,evecfvt,evecsvt,pmat)

  call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')

!!$  if (tresume) tresume=.false.
  write(unitout,'(a)') "Info("//trim(thisnam)//"): momentum matrix elements &
       &finished"

end subroutine writepmattd
!EOC
