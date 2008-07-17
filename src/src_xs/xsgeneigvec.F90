
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xsgeneigvec
  use modmain
  use modmpi
  use modxs
  use m_writegqpts
  use m_filedel
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='xsgeneigvec'
  real(8) :: vqlt(3)
  integer :: iq,qi,qf
  logical, external :: tqgamma
  ! initialize universal variables
  call init0
  call init1
  call init2xs
  ! SCF allready parallelized for k-point set
  qi=1
  qf=nqpt
  ! add extra q-point for if files for q=0 are to be calculated
  if (tq0ev) qi=0
  ! if first Q-point is Gamma-point we copy files
  if (tqgamma(1)) qi=1
  if (tscreen) then
     qi=0
     qf=0
  end if
  ! write q-points
  if (rank.eq.0) call writeqpts
  ! calculate eigenvectors for each q-point (k+q point set)
  do iq=qi,qf
     if (.not.tscreen) &
          call genfilname(iqmt=max(0,iq),setfilext=.true.)
     vqlt(:)=0.d0
     if ((iq.ne.0).and.(rank.eq.0)) then
        vqlt(:)=vql(:,iq)
        call writegqpts(iq)
     end if
     ! write eigenvectors, -values, occupancies and contracted MT coefficients
     call writeevec(vqlt,qvkloff(1,iq),filext)
     if (.not.tscreen) then
        write(unitout,'(a,i8,3g18.10)') 'Info('//thisnam//'): eigenvectors &
             &generated for Q-point:', iq, vqlt(:)
     else
        write(unitout,'(a)') 'Info('//thisnam//'): eigenvectors &
             &generated for screening/screened interaction:'
     end if
     if (rank.eq.0) then
        ! safely remove unnecessary files
        call filedel('EQATOMS'//trim(filext))
        call filedel('EVALCORE'//trim(filext))
        call filedel('FERMIDOS'//trim(filext))
        call filedel('GEOMETRY'//trim(filext))
        call filedel('LATTICE'//trim(filext))
        call filedel('IADIST'//trim(filext))
        call filedel('LINENGY'//trim(filext))
        call filedel('SYMCRYS'//trim(filext))
        call filedel('SYMLAT'//trim(filext))
        call filedel('SYMSITE'//trim(filext))
        call filedel('TOTENERGY'//trim(filext))
        call filedel('EVALFV'//trim(filext))
        call filedel('RMSDVEFF'//trim(filext))
     end if
     ! end loop over q-points
  end do
  ! generate symbolic links for Gamma-Q-point files (probably not ISO-Fortran)
  ! alternatively one could write these files again one by one
  if ((rank.eq.0).and.tqgamma(1).and.(.not.tscreen)) then
     write(unitout,'(a)') 'Info('//thisnam//'): First Q-point is &
          &Gamma-point'
#ifdef ISO
     call copyfilesq0
#elseif
     call linkfilesq0
#endif

     call system('ln -sf EVECFV_QMT001.OUT  EVECFV_QMT000.OUT')
     call system('ln -sf EVECSV_QMT001.OUT  EVECSV_QMT000.OUT')
     call system('ln -sf EVALSV_QMT001.OUT  EVALSV_QMT000.OUT')
     call system('ln -sf OCCSV_QMT001.OUT   OCCSV_QMT000.OUT')
     call system('ln -sf APWCMT_QMT001.OUT  APWCMT_QMT000.OUT')
     call system('ln -sf LOCMT_QMT001.OUT   LOCMT_QMT000.OUT')
     call system('ln -sf EIGVAL_QMT001.OUT  EIGVAL_QMT000.OUT')
     call system('ln -sf KPOINTS_QMT001.OUT KPOINTS_QMT000.OUT')
     call system('ln -sf EFERMI_QMT001.OUT  EFERMI_QMT000.OUT')
  end if
  call barrier
  write(unitout,'(a)') "Info("//trim(thisnam)//"): generation of &
       &eigenvectors finished"
  call genfilname(setfilext=.true.)
end subroutine xsgeneigvec

subroutine copyfilesq0
  use modmain
  use modxs
  implicit none
  ! local variables
  integer, parameter :: iq=1
  integer :: ik
  real(8), allocatable :: evecfvt(:,:,:)
  complex(8), allocatable :: apwlm(:,:,:,:),lolm(:,:,:,:)
  allocate(evecfvt(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
  allocate(apwlm(nstfv,apwordmax,lmmaxapw,natmtot))
  allocate(lolm(nstfv,nlomax,-lolmax:lolmax,natmtot))
  do ik=1,nkpt
     ! read files
     filext='QMT001.OUT'
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfvt)
     call getevecsv(vkl(1,ik),evecsv)     
     call getevalsv(vkl(1,ik),evalsv(1,ik))
     call getoccsv(vkl(1,ik),occsv(1,ik))
     call getapwcmt(iq,ik,1,nstfv,lmaxapw,apwlm)
     call getlocmt(iq,ik,1,nstfv,lolm)
     ! write files
     filext='QMT000.OUT'
     call putevecfv(ik,evecfvt)
     call putevecsv(ik,evecsv)
     call putevalsv(ik,evalsv(1,ik))
     call putoccsv(ik,occsv(1,ik))
     call putapwcmt('APWCMT_QMT000.OUT',ik,vkl(:,ik),vql(:,iq),apwcmt)
     call putlocmt('LOCMT_QMT000.OUT',ik,vkl(:,ik),vql(:,iq),locmt)
  end do
  ! read files
  filext='QMT001.OUT'
  call readfermi
  ! write files
  filext='QMT000.OUT'
  call writeeval
  call writefermi
  deallocate(evecfvt,evecsv,apwlm,lolm)
end subroutine copyfilesq0

subroutine linkfilesq0
  implicit none
  call system('ln -sf EVECFV_QMT001.OUT  EVECFV_QMT000.OUT')
  call system('ln -sf EVECSV_QMT001.OUT  EVECSV_QMT000.OUT')
  call system('ln -sf EVALSV_QMT001.OUT  EVALSV_QMT000.OUT')
  call system('ln -sf OCCSV_QMT001.OUT   OCCSV_QMT000.OUT')
  call system('ln -sf APWCMT_QMT001.OUT  APWCMT_QMT000.OUT')
  call system('ln -sf LOCMT_QMT001.OUT   LOCMT_QMT000.OUT')
  call system('ln -sf EIGVAL_QMT001.OUT  EIGVAL_QMT000.OUT')
  call system('ln -sf KPOINTS_QMT001.OUT KPOINTS_QMT000.OUT')
  call system('ln -sf EFERMI_QMT001.OUT  EFERMI_QMT000.OUT')
end subroutine linkfilesq0

