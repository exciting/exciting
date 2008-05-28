
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writeevec(vq,voff,filxt)
  use modmain
  use modmpi
  use modxs
  use m_gndstateq
  use m_filedel
  implicit none
  ! arguments
  real(8), intent(in) :: vq(3),voff(3)
  character(*), intent(in) :: filxt
  ! local variables
  integer :: ik,j
  ! read from STATE.OUT exclusively
  isreadstate0=.true.
  ! SCF calculation with one cycle
  call gndstateq(voff,filxt)
  if (allocated(evecfv)) deallocate(evecfv)
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  if (allocated(apwalm)) deallocate(apwalm)
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
  allocate(apwdlm(nstsv,apwordmax,lmmaxapw,natmtot))
  allocate(lodlm(nstsv,nlomax,-lolmax:lolmax,natmtot))
  ! delete existing coefficients files
  if (rank.eq.0) call filedel('APWDLM'//trim(filxt))
  if (rank.eq.0) call filedel('LODLM'//trim(filxt))
  call genparidxran('k',nkpt)
  do ik=kpari,kparf
     apwdlm(:,:,:,:)=zzero
     lodlm(:,:,:,:)=zzero
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
     call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1), &
          apwalm)
     call genapwcmt(lmaxapw,ngk(ik,1),1,nstfv,apwalm,evecfv,apwdlm)
     call genlocmt(ngk(ik,1),1,nstfv,evecfv,lodlm)
     do j=0,procs-1
        if (rank.eq.j) then
           call putapwcmt('APWDLM'//trim(filxt),ik,vkl(1,ik),vq,apwdlm)
           call putlocmt('LODLM'//trim(filxt),ik,vkl(1,ik),vq,lodlm)
        end if
        call barrier
     end do
  end do
  call endloopbarrier(nkpt,procs)
  isreadstate0=.false.
  deallocate(evecfv,apwalm,apwdlm)
end subroutine writeevec
