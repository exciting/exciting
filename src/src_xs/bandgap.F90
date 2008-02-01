
! Copyright (C) 2002-2006 S. Sagmeister, J. K. Dewhurst, S. Sharma and 
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_bandgap
  implicit none
contains

!BOP
! !ROUTINE: bandgapgrid
! !INTERFACE:
  subroutine bandgap(n,e,ef,egf,ego,ikgf,ikgo,istho)
! !USES:
    use modmain
! !DESCRIPTION:
!   Determines the fundamental and optical band gap if present.
!
! !REVISION HISTORY:
!   Created July 2007 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    integer, intent(in) :: n
    real(8), intent(in) :: e(:,:),ef
    real(8), intent(out) :: egf, ego
    integer, intent(out) :: ikgf(2), ikgo, istho
    ! local variables
    integer ik
    integer :: klu1(1),kho1(1),kluho1(1)
    ! allocatable arrays
    real(8), allocatable :: de(:),eho(:),elu(:)

    allocate(de(nkpt),eho(nkpt),elu(nkpt))
    do ik=1,n
       istho=count(e(:,ik)<=ef)
       if (istho == nstfv) goto 10
       eho(ik)=e(istho,ik)
       elu(ik)=e(istho+1,ik)
       de(ik)=elu(ik)-eho(ik)
    end do
    kho1=maxloc(eho)
    klu1=minloc(elu)
    kluho1=minloc(elu-eho)
    ikgf(1)=kho1(1)
    ikgf(2)=klu1(1)
    ikgo=kluho1(1)
    egf=elu(ikgf(2))-eho(ikgf(1))
    ego=elu(ikgo)-eho(ikgo)
    return
10  continue
    ! all states occupied
    egf=0.0
    ego=0.0
    ikgf=1
    ikgo=1
    istho=nstfv
    deallocate(de,eho,elu)

  end subroutine bandgap
!EOP
end module m_bandgap


subroutine writebandgap
  use modmain
  use modxs
  use m_bandgap
  implicit none
  ! local variables
  real(8) :: egf,ego,eho,elu,de,v1(3),v2(3)
  integer :: ikgf(2),ikgo,istho,iv,ik
  logical, allocatable :: done(:)
  real(8), external :: r3dist

  ! calculate bandgap
  call bandgap(nkpt,evalsv,efermi,egf,ego,ikgf,ikgo,istho)
  
  ! write band gap to file
  if (egf == 0.d0) goto 10
  if (task == 23) then
     open(50,file='BANDGAP_GRID.OUT',form='formatted',action='write', &
          status='replace')
     write(50,'(a)') 'Band gaps determined from energies on grid'
     write(50,'(a,3i6)') ' k-point grid:',ngridk
     write(50,'(a,3f12.6)') ' k-point offset:',vkloff
     write(50,'(a,3f12.6)') ' k-point shift :',vkloff/dble(ngridk)
  else if (task == 20) then
     open(50,file='BANDGAP.OUT',form='formatted',action='write', &
          status='replace')
     write(50,'(a)') 'Band gaps determined from energies on k-point path'
     write(50,'(a)') 'energies and differences on vertex locations'
     write(50,'(a)') 'iv, ik, vkl, homo, lumo, de, de[eV] below'
     allocate(done(nvp1d))
     done(:)=.false.
     do iv=1,nvp1d
        v1(:)=vvlp1d(:,iv)
        do ik=1,nkpt
           v2(:)=vkl(:,ik)
           if(.not.done(iv).and.(r3dist(v1,v2)<epslat)) then
              eho=evalsv(istho,ik)
              elu=evalsv(istho+1,ik)
              de=elu-eho
              write(50,'(2i6,7f12.3)') iv,ik,vkl(:,ik),eho,elu,de,h2ev*de
              done(iv)=.true.
           end if
        end do
     end do
  end if
  write(50,'(a,g16.8,a,g16.8,a)') 'fundamental gap: ',egf, &
       ' (',h2ev*egf,' eV )'
  write(50,'(a,i9,4f12.6)') ' k-point (homo), energy: ',ikgf(1),&
       vkl(:,ikgf(1)),evalsv(istho,ikgf(1))
  write(50,'(a,i9,4f12.6)') ' k-point (lumo), energy: ',ikgf(2),&
       vkl(:,istho+1),evalsv(istho+1,ikgf(2))
  write(50,'(a,g16.8,a,g16.8,a)') 'optical gap    : ',ego, &
       ' (',h2ev*ego,' eV )'
  write(50,'(a,i9,4f12.6)') ' k-point (homo), energy: ',ikgo, &
       vkl(:,ikgo),evalsv(istho,ikgo)
  write(50,'(a,i9,4f12.6)') ' k-point (lumo), energy: ',ikgo, &
       vkl(:,ikgo),evalsv(istho+1,ikgo)
  close(50)
  return

10 continue
  write(50,'(a)') 'No bandgaps found.'

end subroutine writebandgap


subroutine writebandgapgrid
  use modmain
  use modxs
  use m_genfilname
  implicit none
  ! local variables
  integer :: ik
  ! initialise universal variables
  if (calledxs.eq.1) call init0
  ! file extension for q-point
  call genfilname(iqmt=0,setfilext=.true.)
  call init1
  ! read Fermi energy from file
  call readfermi
  do ik=1,nkpt
     ! get the eigenvectors and values from EIGVEC.OUT
     call getevalsv(vkl(1,ik),evalsv(1,ik))
     call getoccsv(vkl(1,ik),occsv(1,ik))
  end do
  call writebandgap
  call genfilname(revertfilext=.true.)
end subroutine writebandgapgrid

