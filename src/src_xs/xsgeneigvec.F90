
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
     call system('ln -sf EVECFV_QMT001.OUT  EVECFV_QMT000.OUT')
     call system('ln -sf EVECSV_QMT001.OUT  EVECSV_QMT000.OUT')
     call system('ln -sf EVALSV_QMT001.OUT  EVALSV_QMT000.OUT')
     call system('ln -sf OCCSV_QMT001.OUT   OCCSV_QMT000.OUT')
     call system('ln -sf APWDLM_QMT001.OUT  APWDLM_QMT000.OUT')
     call system('ln -sf LODLM_QMT001.OUT   LODLM_QMT000.OUT')
     call system('ln -sf EIGVAL_QMT001.OUT  EIGVAL_QMT000.OUT')
     call system('ln -sf KPOINTS_QMT001.OUT KPOINTS_QMT000.OUT')
     call system('ln -sf EFERMI_QMT001.OUT  EFERMI_QMT000.OUT')
  end if
  call barrier
  write(unitout,'(a)') "Info("//trim(thisnam)//"): generation of &
       &eigenvectors finished"
  call genfilname(setfilext=.true.)
end subroutine xsgeneigvec

subroutine writeevec(vq,voff,filxt)
  use modmain
  use modxs
  use m_gndstateq
  use m_genapwcmt
  use m_getunit
  implicit none
  ! arguments
  real(8), intent(in) :: vq(3),voff(3)
  character(*), intent(in) :: filxt
  ! local variables
  integer :: ik,reclapw,recllo
  ! read from STATE.OUT exclusively
  isreadstate0=.true.
  call gndstateq(voff,filxt)
  if (allocated(evecfv)) deallocate(evecfv)
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  if (allocated(apwalm)) deallocate(apwalm)
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
  allocate(apwdlm(nstsv,apwordmax,lmmaxapw,natmtot))
  allocate(lodlm(nstsv,nlomax,-lolmax:lolmax,natmtot))
  call getunit(unit1)
  inquire(iolength=reclapw) vq,vkl(:,1),apwdlm
  open(unit1,file='APWDLM'//trim(filxt),action='write',&
       form='unformatted',status='replace',access='direct',recl=reclapw)
  call getunit(unit2)
  inquire(iolength=recllo) vq,vkl(:,1),lodlm
  open(unit2,file='LODLM'//trim(filxt),action='write',&
       form='unformatted',status='replace',access='direct',recl=recllo)
  do ik=1,nkpt
     call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
     call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1), &
          apwalm)
     call genapwcmt(lmaxapw,ngk(ik,1),1,nstsv,apwalm,evecfv,apwdlm)
     write(unit1,rec=ik) vq,vkl(:,ik),apwdlm
     call genlocmt(ngk(ik,1),1,nstsv,evecfv,lodlm)
     write(unit2,rec=ik) vq,vkl(:,ik),lodlm
  end do
  close(unit1)
  close(unit2)
  isreadstate0=.false.
  deallocate(evecfv,apwalm,apwdlm)
end subroutine writeevec
