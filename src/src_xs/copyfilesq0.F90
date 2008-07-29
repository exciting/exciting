
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine copyfilesq0
#ifdef ISO
  use modmain
  use modxs
  use m_getapwcmt
  use m_getlocmt
  implicit none
  ! local variables
  integer, parameter :: iq=1
  integer :: ik
  complex(8), allocatable :: evecfvt(:,:,:)
  complex(8), allocatable :: apwlm(:,:,:,:),lolm(:,:,:,:)
  allocate(evecfvt(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
  allocate(apwlm(nstfv,apwordmax,lmmaxapw,natmtot))
  allocate(lolm(nstfv,nlomax,-lolmax:lolmax,natmtot))
  do ik=1,nkpt
     ! read files
     filext='_QMT001.OUT'
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfvt)
     call getevecsv(vkl(1,ik),evecsv)     
     call getevalsv(vkl(1,ik),evalsv(1,ik))
     call getoccsv(vkl(1,ik),occsv(1,ik))
     call getapwcmt(iq,ik,1,nstfv,lmaxapw,apwlm)
     call getlocmt(iq,ik,1,nstfv,lolm)
     ! write files
     filext='_QMT000.OUT'
     call putevecfv(ik,evecfvt)
     call putevecsv(ik,evecsv)
     call putevalsv(ik,evalsv(1,ik))
     call putoccsv(ik,occsv(1,ik))
     call putapwcmt('APWCMT_QMT000.OUT',ik,vkl(:,ik),vql(:,iq),apwlm)
     call putlocmt('LOCMT_QMT000.OUT',ik,vkl(:,ik),vql(:,iq),lolm)
  end do
  ! read files
  filext='_QMT001.OUT'
  call readfermi
  ! write files
  filext='_QMT000.OUT'
  call writeeval
  call writefermi
  deallocate(evecfvt)
  deallocate(evecsv)
  deallocate(apwlm)
  deallocate(lolm)
#endif
#ifndef ISO
  call system('ln -sf EVECFV_QMT001.OUT  EVECFV_QMT000.OUT')
  call system('ln -sf EVECSV_QMT001.OUT  EVECSV_QMT000.OUT')
  call system('ln -sf EVALSV_QMT001.OUT  EVALSV_QMT000.OUT')
  call system('ln -sf OCCSV_QMT001.OUT   OCCSV_QMT000.OUT')
  call system('ln -sf APWCMT_QMT001.OUT  APWCMT_QMT000.OUT')
  call system('ln -sf LOCMT_QMT001.OUT   LOCMT_QMT000.OUT')
  call system('ln -sf EIGVAL_QMT001.OUT  EIGVAL_QMT000.OUT')
  call system('ln -sf KPOINTS_QMT001.OUT KPOINTS_QMT000.OUT')
  call system('ln -sf EFERMI_QMT001.OUT  EFERMI_QMT000.OUT')
#endif
end subroutine copyfilesq0
