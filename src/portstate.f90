
! Copyright (C) 2002-2005 S. Sagmeister J. K. Dewhurst, S. Sharma and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: portstate
! !INTERFACE:
subroutine portstate(tb2a)
  ! !USES:
  ! !DESCRIPTION:
  !   Toggle file format of {\tt STATE.OUT}. If tb2a is true an ASCII
  !   file with the name {\tt STATE_ASC.OUT} is generated and the data
  !   from {\tt STATE.OUT} is transferred. If tb2a is false the conversion
  !   goes in the other direction. based upon the routines {\tt readstate}
  !   and {\tt writestate}.
  !   
  ! !REVISION HISTORY:
  !   Created 2007 (Sagmeister)
  !EOP
  !BOC
  implicit none
  ! arguments
  logical, intent(in) :: tb2a
  ! local variables
  logical spinpol_
  integer natmtot,is
  integer version_(3),nspecies_,lmmaxvr_,nrmtmax_
  integer natoms_,ngrid_(3)
  integer ngrtot_,ngvec_,ndmag_,nspinor_,ldapu_,lmmaxlu_
  ! allocatable arrays
  integer, allocatable :: nrmt_(:)
  real(8), allocatable :: spr_(:,:)
  real(8), allocatable :: rhomt_(:,:,:)
  real(8), allocatable :: rhoir_(:)
  real(8), allocatable :: vclmt_(:,:,:)
  real(8), allocatable :: vclir_(:)
  real(8), allocatable :: vxcmt_(:,:,:)
  real(8), allocatable :: vxcir_(:)
  real(8), allocatable :: veffmt_(:,:,:)
  real(8), allocatable :: veffir_(:)
  real(8), allocatable :: magmt_(:,:,:,:)
  real(8), allocatable :: magir_(:,:)
  real(8), allocatable :: bxcmt_(:,:,:,:)
  real(8), allocatable :: bxcir_(:,:)
  complex(8), allocatable :: veffig_(:)
  complex(8), allocatable :: vmatlu_(:,:,:,:,:)
  if (tb2a) then
     open(50,file='STATE.OUT',action='READ',form='UNFORMATTED', &
          status='OLD')
     open(51,file='STATE.xml',action='WRITE',form='FORMATTED',&
          status='replace')
     read(50) version_
     read(50) spinpol_
     read(50) nspecies_
     read(50) lmmaxvr_
     read(50) nrmtmax_
     write(51,'(a)') '<?xml version="1.0"?>'
     write(51,'(a)') '<data name="version" type="integer" dimension="1" &
          &shape="3">'
     write(51,*) version_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="spinpol" type="logical" dimension="1" &
          &shape="1">'
     write(51,*) spinpol_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="nspecies" type="integer" dimension="1" &
          &shape="1">'
     write(51,*) nspecies_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="lmmaxvr" type="integer" dimension="1" &
          &shape="1">'
     write(51,*) lmmaxvr_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="nrmtmax" type="integer" dimension="1" &
          &shape="1">'
     write(51,*) nrmtmax_
     write(51,'(a)') '</data>'
  else
     open(50,file='STATE.xml',action='READ',form='FORMATTED', &
          status='OLD')
     open(51,file='STATE.OUT',action='WRITE',form='UNFORMATTED', &
          status='replace')
     read(50,*)
     read(50,*)
     read(50,*) version_
     read(50,*)
     read(50,*)
     read(50,*) spinpol_
     read(50,*)
     read(50,*)
     read(50,*) nspecies_
     read(50,*)
     read(50,*)
     read(50,*) lmmaxvr_
     read(50,*)
     read(50,*)
     read(50,*) nrmtmax_
     read(50,*)
     write(51) version_
     write(51) spinpol_
     write(51) nspecies_
     write(51) lmmaxvr_
     write(51) nrmtmax_
  end if
  allocate(spr_(nrmtmax_,nspecies_))
  allocate(nrmt_(nspecies_))
  if (tb2a) then
     natmtot=0
     do is=1,nspecies_
        read(50) natoms_
        read(50) nrmt_(is)
        read(50) spr_(1:nrmt_(is),is)
        write(51,'(a)') '<data name="natoms" type="integer" dimension="1" &
             &shape="1" index="species" indexval="'//trim(i2str(is))//'">'
        write(51,*) natoms_
        write(51,'(a)') '</data>'
        write(51,'(a)') '<data name="nrmt" type="integer" dimension="1" &
             &shape="1" index="species" indexval="'//trim(i2str(is))//'">'
        write(51,*) nrmt_(is)
        write(51,'(a)') '</data>'
        write(51,'(a)') '<data name="spr" type="real(8)" dimension="1" &
             &shape="'//trim(i2str(nrmt_(is)))//'" index="species" indexval="'&
             //trim(i2str(is))//'">'
        write(51,*) spr_(1:nrmt_(is),is)
        write(51,'(a)') '</data>'
        natmtot=natmtot+natoms_
     end do
     read(50) ngrid_
     read(50) ngvec_
     read(50) ndmag_
     read(50) nspinor_
     read(50) ldapu_
     read(50) lmmaxlu_
     write(51,'(a)') '<data name="ngrid" type="integer" dimension="1" &
          &shape="3">'
     write(51,*) ngrid_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="ngvec" type="integer" dimension="1" &
          &shape="1">'
     write(51,*) ngvec_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="ndmag" type="integer" dimension="1" &
          &shape="1">'
     write(51,*) ndmag_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="ndmag" type="integer" dimension="1" &
          &shape="1">'
     write(51,*) ndmag_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="nspinor" type="integer" dimension="1" &
          &shape="1">'
     write(51,*) nspinor_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="ldapu" type="integer" dimension="1" &
          &shape="1">'
     write(51,*) ldapu_
     write(51,'(a)') '</data>'
          write(51,'(a)') '<data name="lmmaxlu" type="integer" dimension="1" &
          &shape="1">'
     write(51,*) lmmaxlu_
     write(51,'(a)') '</data>'
  else
     natmtot=0
     do is=1,nspecies_
        read(50,*)
        read(50,*) natoms_
        read(50,*)
        read(50,*)
        read(50,*) nrmt_(is)
        read(50,*)
        read(50,*)
        read(50,*) spr_(1:nrmt_(is),is)
        read(50,*)
        write(51) natoms_
        write(51) nrmt_(is)
        write(51) spr_(1:nrmt_(is),is)
        natmtot=natmtot+natoms_
     end do
     read(50,*)
     read(50,*) ngrid_
     read(50,*)
     read(50,*)
     read(50,*) ngvec_
     read(50,*)
     read(50,*)
     read(50,*) ndmag_
     read(50,*)
     read(50,*)
     read(50,*) nspinor_
     read(50,*)
     read(50,*)
     read(50,*) ldapu_
     read(50,*)
     read(50,*)
     read(50,*) lmmaxlu_
     read(50,*)
     write(51) ngrid_
     write(51) ngvec_
     write(51) ndmag_
     write(51) nspinor_
     write(51) ldapu_
     write(51) lmmaxlu_
  end if
  ngrtot_=ngrid_(1)*ngrid_(2)*ngrid_(3)
  allocate(rhomt_(lmmaxvr_,nrmtmax_,natmtot))
  allocate(rhoir_(ngrtot_))
  allocate(vclmt_(lmmaxvr_,nrmtmax_,natmtot))
  allocate(vclir_(ngrtot_))
  allocate(vxcmt_(lmmaxvr_,nrmtmax_,natmtot))
  allocate(vxcir_(ngrtot_))
  allocate(veffmt_(lmmaxvr_,nrmtmax_,natmtot))
  allocate(veffir_(ngrtot_))
  allocate(veffig_(ngvec_))
  if (spinpol_) then
     allocate(magmt_(lmmaxvr_,nrmtmax_,natmtot,ndmag_))
     allocate(magir_(ngrtot_,ndmag_))
     allocate(bxcmt_(lmmaxvr_,nrmtmax_,natmtot,ndmag_))
     allocate(bxcir_(ngrtot_,ndmag_))
  end if
  if (ldapu_.ne.0) then
     allocate(vmatlu_(lmmaxlu_,lmmaxlu_,nspinor_,nspinor_,natmtot))
  end if
  if (tb2a) then
     ! read muffin-tin density
     read(50) rhomt_,rhoir_
     ! read Coulomb potential (spin independent)
     read(50) vclmt_,vclir_
     ! read exchange-correlation potential
     read(50) vxcmt_,vxcir_
     ! read effective potential
     read(50) veffmt_,veffir_,veffig_
     ! write the density
     write(51,'(a)') '<data name="rhomt" type="real(8)" &
          &dimension="3" shape="'//&
          trim(i2str(lmmaxvr_))//','//&
          trim(i2str(nrmtmax_))//','//&
          trim(i2str(natmtot))//'">'
     write(51,*) rhomt_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="rhoir" type="real(8)" &
          &dimension="1" shape="'//&
          trim(i2str(ngrtot_))//'">'
     write(51,*) rhoir_
     write(51,'(a)') '</data>'
     ! write the Coulomb potential
     write(51,'(a)') '<data name="vclmt" type="real(8)" &
          &dimension="3" shape="'//&
          trim(i2str(lmmaxvr_))//','//&
          trim(i2str(nrmtmax_))//','//&
          trim(i2str(natmtot))//'">'
     write(51,*) vclmt_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="vclir" type="real(8)" &
          &dimension="1" shape="'//&
          trim(i2str(ngrtot_))//'">'
     write(51,*) vclir_
     write(51,'(a)') '</data>'
     ! write the exchange-correlation potential
     write(51,'(a)') '<data name="vxcmt" type="real(8)" &
          &dimension="3" shape="'//&
          trim(i2str(lmmaxvr_))//','//&
          trim(i2str(nrmtmax_))//','//&
          trim(i2str(natmtot))//'">'
     write(51,*) vxcmt_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="vxcir" type="real(8)" &
          &dimension="1" shape="'//&
          trim(i2str(ngrtot_))//'">'
     write(51,*) vxcir_
     write(51,'(a)') '</data>'
     ! write the effective potential
     write(51,'(a)') '<data name="veffmt" type="real(8)" &
          &dimension="3" shape="'//&
          trim(i2str(lmmaxvr_))//','//&
          trim(i2str(nrmtmax_))//','//&
          trim(i2str(natmtot))//'">'
     write(51,*) veffmt_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="veffir" type="real(8)" &
          &dimension="1" shape="'//&
          trim(i2str(ngrtot_))//'">'
     write(51,*) veffir_
     write(51,'(a)') '</data>'
     write(51,'(a)') '<data name="veffig" type="complex(8)" &
          &dimension="1" shape="'//&
          trim(i2str(ngvec_))//'">'
     write(51,*) veffig_
     write(51,'(a)') '</data>'
     if (spinpol_) then
        ! read magnetisation and effective field
        read(50) magmt_,magir_
        read(50) bxcmt_,bxcir_
        ! write the magnetisation and effective magnetic fields
        write(51,'(a)') '<data name="magmt" type="real(8)" &
             &dimension="4" shape="'//&
          trim(i2str(lmmaxvr_))//','//&
          trim(i2str(nrmtmax_))//','//&
          trim(i2str(natmtot))//','//&
          trim(i2str(ndmag_))//'">'
        write(51,*) magmt_
        write(51,'(a)') '</data>'
        write(51,'(a)') '<data name="magir" type="real(8)" &
             &dimension="2" shape="'//&
             trim(i2str(ngrtot_))//','//&
             trim(i2str(ndmag_))//'">'
        write(51,*) magir_
        write(51,'(a)') '</data>'
        write(51,'(a)') '<data name="bxcmt" type="real(8)" &
             &dimension="4" shape="'//&
          trim(i2str(lmmaxvr_))//','//&
          trim(i2str(nrmtmax_))//','//&
          trim(i2str(natmtot))//','//&
          trim(i2str(ndmag_))//'">'
        write(51,*) bxcmt_
        write(51,'(a)') '</data>'
        write(51,'(a)') '<data name="bxcir" type="real(8)" &
             &dimension="2" shape="'//&
             trim(i2str(ngrtot_))//','//&
             trim(i2str(ndmag_))//'">'
        write(51,*) bxcir_
        write(51,'(a)') '</data>'
     end if
     if (ldapu_.ne.0) then
        ! read the LDA+U potential matrix elements
        read(50) vmatlu_
        ! write the LDA+U potential matrix elements
        write(51,'(a)') '<data name="vmatlu" type="complex(8)" &
             &dimension="5" shape="'//&
             trim(i2str(lmmaxlu_))//','//&
             trim(i2str(lmmaxlu_))//','//&
             trim(i2str(nspinor_))//','//&
             trim(i2str(nspinor_))//','//&
             trim(i2str(natmtot))//'">'
        write(51,*) vmatlu_
        write(51,'(a)') '</data>'        
     end if
  else
     ! read muffin-tin density
     read(50,*)
     read(50,*) rhomt_
     read(50,*)
     read(50,*)
     read(50,*) rhoir_
     read(50,*)
     ! read Coulomb potential (spin independent)
     read(50,*)
     read(50,*) vclmt_
     read(50,*)
     read(50,*)
     read(50,*) vclir_
     read(50,*)
     ! read exchange-correlation potential
     read(50,*)
     read(50,*) vxcmt_
     read(50,*)
     read(50,*)
     read(50,*) vxcir_
     read(50,*)
     ! read effective potential
     read(50,*)
     read(50,*) veffmt_
     read(50,*)
     read(50,*)
     read(50,*) veffir_
     read(50,*)
     read(50,*)
     read(50,*) veffig_
     read(50,*)
     ! write the density
     write(51) rhomt_,rhoir_
     ! write the Coulomb potential
     write(51) vclmt_,vclir_
     ! write the exchange-correlation potential
     write(51) vxcmt_,vxcir_
     ! write the effective potential
     write(51) veffmt_,veffir_,veffig_
     if (spinpol_) then
        ! read magnetisation and effective field
        read(50,*)
        read(50,*) magmt_
        read(50,*)
        read(50,*)
        read(50,*) magir_
        read(50,*)
        read(50,*)
        read(50,*) bxcmt_
        read(50,*)
        read(50,*)
        read(50,*) bxcir_
        read(50,*)
        ! write the magnetisation and effective magnetic fields
        write(51) magmt_,magir_
        write(51) bxcmt_,bxcir_
     end if
     if (ldapu_.ne.0) then
        ! read the LDA+U potential matrix elements
        read(50,*)
        read(50,*) vmatlu_
        read(50,*)        
        ! write the LDA+U potential matrix elements
        write(51) vmatlu_
     end if
  end if
  close(50)
  close(51)
  deallocate(nrmt_,spr_,rhomt_,rhoir_,vclmt_,vclir_)
  deallocate(vxcmt_,vxcir_,veffmt_,veffir_,veffig_)
  if (spinpol_) deallocate(magmt_,magir_,bxcmt_,bxcir_)
  if (tb2a) then
     write(*,*)
     write(*,'("Info(portstate): generated portable ASCII state file &
          &STATE.xml from STATE.OUT file")')
     write(*,*)
  else
     write(*,*)
     write(*,'("Info(portstate): generated STATE.OUT file from portable &
          &ASCII state file STATE.xml")')
     write(*,*)
  end if
contains
  character(256) function i2str(i)
    ! arguments
    integer, intent(in) :: i
    ! local variables
    character(1024) :: str
    write(str,*) i
    i2str=trim(adjustl(str))
  end function i2str
end subroutine portstate
!EOC
