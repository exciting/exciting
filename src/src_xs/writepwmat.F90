
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writepwmat
! !INTERFACE:
subroutine writepwmat
! !USES:
  use modmain
  use modxs
! !DESCRIPTION:
!   Calculates the matrix elements of the plane wave $exp(-i(\mathbf{G}+
!   mathbf{q})\mathbf{r})$ using routine {\tt genpwmat} and writes them to
!   direct access file {\tt PWMAT.OUT}.
!
! !REVISION HISTORY:
!   Created November 2007 (Sagmeister)
!EOP
!BOC
  implicit none
  ! local variables
  integer, parameter :: iq=1
  integer ik,ikp,recl,isymkp,igq,un
  real(8) :: vpl(3),vkpl(3)
  complex(8), allocatable :: apwalmk(:,:,:,:),apwalmkp(:,:,:,:)
  complex(8), allocatable :: evecfvk(:,:),evecfvkp(:,:)
  complex(8), allocatable :: evecsvk(:,:),evecsvkp(:,:)
  complex(8), allocatable :: pwmat(:,:,:)
  ! initialise universal variables
  call init0
  call init1
  call init2xs
  call findocclims(0,istocc0,istocc,istunocc0,istunocc,isto0,isto,istu0,istu)
  allocate(apwalmk(ngkmax,apwordmax,lmmaxapw,natmtot))
  allocate(apwalmkp(ngkmax,apwordmax,lmmaxapw,natmtot))
  allocate(evecfvk(nmatmax,nstfv))
  allocate(evecfvkp(nmatmax,nstfv))
  allocate(evecsvk(nstsv,nstsv))
  allocate(evecsvkp(nstsv,nstsv))
  ! allocate the momentum matrix elements array
  allocate(pwmat(ngq(iq),nstsv,nstsv))


  ! allocate matrix elements array
  if (allocated(xiou)) deallocate(xiou)
  if (allocated(xiuo)) deallocate(xiuo)
  allocate(xiou(nstval,nstcon,ngq(iq)))
  allocate(xiuo(nstcon,nstval,ngq(iq)))


  ! read in the density and potentials from file
  call readstate
  ! find the new linearisation energies
  call linengy
  ! generate the APW radial functions
  call genapwfr
  ! generate the local-orbital radial functions
  call genlofr
  ! find the record length
  inquire(iolength=recl) pwmat
!!$  open(50,file='PWMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
!!$       status='REPLACE',recl=recl)
  open(50,file='PWMAT_ASC.OUT',action='WRITE',form='FORMATTED',status='REPLACE')
  do ik=1,nkpt
     write(*,*) 'Info(writepwmat): ik',ik
     ! get the eigenvectors from file for k-point
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfvk)
     call getevecsv(vkl(1,ik),evecsvk)
     ! find the matching coefficients for k-point
     call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalmk)
     ! find k-point equivalent to k+q
     vpl(:)=vql(:,iq)+vkl(:,ik)
     call findkpt(vpl,isymkp,ikp)
write(*,*) 'ik,ikp',ik,ikp
     vkpl(:)=vkl(:,ikp)
     ! get the eigenvectors from file for kp-point
     call getevecfv(vkl(1,ikp),vgkl(1,1,ikp,1),evecfvkp)
     call getevecsv(vkl(1,ikp),evecsvkp)
     ! find the matching coefficients for kp-point
     call match(ngk(ikp,1),gkc(1,ikp,1),tpgkc(1,1,ikp,1),sfacgk(1,1,ikp,1), &
          apwalmkp)
     ! calculate the matrix elements of the plane wave
     call genpwmat(vql(1,iq),ngqmax,ngq(iq),gqc(1,iq),igqig(1,iq), &
          ylmgq(1,1,iq),sfacgq(1,1,iq),vkl(1,ik),ngk(ik,1),igkig(1,ik,1), &
          apwalmk,evecfvk,evecsvk,vkl(1,ikp),ngk(ikp,1),igkig(1,ikp,1), &
          apwalmkp,evecfvkp,evecsvkp,pwmat)
     write(*,*) '*** after genpwmat ***'


     do igq=1,ngq(iq)
        xiou(:,:,igq)=pwmat(igq,1:nstval,nstval+1:nstsv)
        xiuo(:,:,igq)=pwmat(igq,nstval+1:nstsv,1:nstval)
     end do

     inquire(iolength=recl) nstval, nstcon, nkpt, ngq(iq), vql(:,iq), &
          vkl(:,ik), xiou,xiuo
     un=51
     open(unit=un,file='EMAT_Q00001.OUT',form='unformatted', &
          action='write',access='direct',recl=recl)
     write(un,rec=ik) nstval, nstcon, nkpt, ngq(iq), vql(:,iq), vkl(:,ik), &
          xiou, xiuo
     close(un)



     ! write the matrix elements to direct-access file
!!$     write(50,rec=ik) pwmat
!!$     write(50,'(i8,3g18.10)') ik,vkl(:,ik)
!!$     do igq=1,ngq(iq)
!!$        do ist=1,nstsv
!!$           do jst=1,nstsv
!!$              write(50,'(4i8,3g18.10)') ik,igq,ist,jst,pwmat(igq,ist,jst), &
!!$                   abs(pwmat(igq,ist,jst))**2
!!$           end do
!!$        end do
!!$     end do
!!$     write(50,*)
  end do
  close(50)
  write(*,*)
  write(*,'("Info(writepwmat):")')
  write(*,'(" matrix elements of the plane wave written to file PWMAT.OUT")')
  write(*,*)
  deallocate(apwalmk,evecfvk,evecsvk,apwalmkp,evecfvkp,evecsvkp,pwmat)
end subroutine writepwmat
!EOC
