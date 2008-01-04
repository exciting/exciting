
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tdgeneigvec
  use modmain
  use modxs
  use modmpi
  use m_gndstateq
  use m_genapwcmt
  use m_getunit
  use m_writegqpts
  use m_filedel
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'tdgeneigvec'
  real(8), parameter :: zero3(3)=0.d0
  integer :: iq,qi,qf,un
  integer :: ik,recl


!@@@@@@@@@@@@@
qi=0; qf=1

call init0

!!$  ! initialize universal variables
!!$  call init0
!!$  call init1
!!$
!!$  ! initialize q-point set
!!$  call init2xs
!!$
!!$  ! SCF allready parallelized for k-point set
!!$  qi=1
!!$  ! add extra q-point for if files for q=0 are to be calculated
!!$  if (tq0ev) qi=0
!!$  qf=nqpt
!!$
!!$  ! write q-points
!!$  if (rank.eq.0) call writeqpts
!!$
!!$  ! allocate arrays for APW expansion coefficients
!!$  allocate(apwdlm(nstsv,apwordmax,lmmaxapwtd,natmtot))

  ! read from STATE.OUT exclusively
  isreadstate0=.true.
  ! calculate eigenvectors for each q-point (k+q point set)
  call getunit(unit1)
  do iq=qi,qf

     ! file extension for q-point
!@@     call genfilname(iq=max(0,iq),setfilext=.true.)
write(filext,'("_Q",i5.5,".OUT")') iq

     ! one more iteration for q=0
     if (iq.eq.0) then
        call gndstateq(vkloff,filext)
        write(unitout,'(a)') 'Info('//thisnam//'): eigenvectors generated &
             &for Gamma q-point'        
     else
        ! write G+q vectors
        call writegqpts(iq)
        ! shift k-mesh with q-point
        call gndstateq(qvkloff(:,iq),filext)
        write(unitout,'(a,i8)') 'Info('//thisnam//'): eigenvectors generated &
             &for q-point:', iq
     end if
stop 'stopped'
     ! store product of eigenvectors with matching coefficients
     if (allocated(evecfv)) deallocate(evecfv)
     if (allocated(apwalm)) deallocate(apwalm)
     allocate(evecfv(nmatmax,nstfv,nspnfv))
     allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
     call getunit(unit1)
     inquire(iolength=recl) vql(:,1),vkl(:,1),apwdlm
     open(unit1,file='APWDLM'//trim(filext),action='write',&
          form='unformatted',status='replace',access='direct',recl=recl)
     do ik=1,nkpt
        call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
        call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1), &
             apwalm)
        call genapwcmt(lmaxapwtd,ngk(ik,1),1,nstsv,apwalm,evecfv,apwdlm)
        if (iq.eq.0) then
           write(unit1,rec=ik) zero3,vkl(:,ik),apwdlm
        else
           write(unit1,rec=ik) vql(:,max(iq,1)),vkl(:,ik),apwdlm
        end if
     end do
     close(unit1)
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
  isreadstate0=.false.

  call genfilname(setfilext=.true.)

  deallocate(apwdlm)

  call getunit(un)
  call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')
!!!  call sleepifc(1)

!!$  if (tresume) tresume=.false.

  write(unitout,'(a)') "Info("//trim(thisnam)//"): generation of &
       &eigenvectors finished"

end subroutine tdgeneigvec
