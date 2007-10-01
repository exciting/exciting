
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

program species
implicit none
! local variables
! order of predictor-corrector polynomial
integer, parameter :: np=4
! maximum number of states
integer, parameter :: maxspst=40
! maximum angular momentum allowed
integer, parameter :: lmax=50
! exchange-correlation type
integer, parameter :: xctype=3
integer, parameter :: xcgrad=0
integer nz,spnst,spnr
integer nrmt,nlx,nlorb,i,l
integer ist1,ist2,ir,iostat
real(8), parameter :: pi=3.1415926535897932385d0
! core-valence cut-off energy
real(8), parameter :: ecvcut=-3.5d0
! semi-core cut-off energy
real(8), parameter :: esccut=-0.5d0
! band offset energy
real(8), parameter :: boe=0.15d0
real(8) spmass,rmt,spzn,sprmin,sprmax,t1
character(256) spsymb,spname
! automatic arrays
logical spcore(maxspst)
integer spn(maxspst),spl(maxspst),spk(maxspst)
real(8) spocc(maxspst),eval(maxspst)
! allocatable arrays
real(8), allocatable :: r(:),rho(:),vr(:),rwf(:,:,:)
real(8), allocatable :: fr(:),gr(:),cf(:,:)
open(40,file='species.dat',action='READ',status='OLD',form='FORMATTED')
10 continue
read(40,*,iostat=iostat) nz
if (iostat.ne.0) stop
read(40,*) spsymb
read(40,*) spname
read(40,*) spmass
read(40,*) rmt
read(40,*) spnst
if (spnst.gt.maxspst) then
  write(*,*)
  write(*,'("Error(species): too many states for species ",A)') trim(spname)
  write(*,*)
  stop
end if
do ist1=1,spnst
  read(40,*) spn(ist1),spl(ist1),spk(ist1),i
  if (ist1.ge.2) then
    if (spn(ist1).lt.spn(ist1-1)) then
      write(*,*)
      write(*,'("Error(species): states improperly ordered")')
      write(*,'(" for species ",A)') trim(spname)
      write(*,*)
      stop
    end if
  end if
  spocc(ist1)=dble(i)
end do
read(40,*)
write(*,'("Info(species): running Z = ",I4,", (",A,")")') nz,trim(spname)
! nuclear charge in units of e
spzn=-dble(nz)
! minimum radial mesh point proportional to 1/sqrt(Z)
sprmin=2.d-6/sqrt(abs(dble(spzn)))
! set the number of radial mesh points proportional to number of nodes
nrmt=100*(spn(spnst)+1)
! find the optimal effective infinity
sprmax=80.d0
do i=1,2
  t1=log(sprmax/sprmin)/log(rmt/sprmin)
  spnr=int(t1*dble(nrmt))
  if (allocated(r)) deallocate(r)
  if (allocated(rho)) deallocate(rho)
  if (allocated(vr)) deallocate(vr)
  if (allocated(rwf)) deallocate(rwf)
  if (allocated(fr)) deallocate(fr)
  if (allocated(gr)) deallocate(gr)
  if (allocated(cf)) deallocate(cf)
  allocate(r(spnr))
  allocate(rho(spnr))
  allocate(vr(spnr))
  allocate(rwf(spnr,2,spnst))
  allocate(fr(spnr))
  allocate(gr(spnr))
  allocate(cf(3,spnr))
! generate the radial mesh
  call radmesh(spnr,nrmt,rmt,sprmin,r)
! solve the Kohn-Sham-Dirac equations for the atom
  call atom(spzn,spnst,spn,spl,spk,spocc,xctype,xcgrad,np,spnr,r,eval,rho,vr, &
   rwf)
  do ir=spnr,1,-1
    if (rho(ir).gt.1.d-20) then
      sprmax=1.5d0*r(ir)
      goto 20
    end if
  end do
20 continue
end do
! check total charge is correct
do ir=1,spnr
  fr(ir)=4.d0*pi*rho(ir)*r(ir)**2
end do
call fderiv(-1,spnr,r,fr,gr,cf)
if (abs(gr(spnr)+spzn).gt.1.d-5) then
  write(*,*)
  write(*,'("Error(species): charge mismatch")')
  write(*,*)
  stop
end if
! find which states belong to core
do ist1=1,spnst
  if (eval(ist1).lt.ecvcut) then
    spcore(ist1)=.true.
  else
    spcore(ist1)=.false.
  end if
end do
! check that the state for same n and l but different k is also core
do ist1=1,spnst
  if (spcore(ist1)) then
    do ist2=1,spnst
      if ((spn(ist1).eq.spn(ist2)).and.(spl(ist1).eq.spl(ist2))) &
       spcore(ist2)=.true.
    end do
  end if
end do
! find the total number of local orbitals
nlorb=4
! search from lowest to highest energy
do ist1=1,spnst
  if (.not.spcore(ist1)) then
    if ((spl(ist1).eq.0).or.(spl(ist1).eq.spk(ist1))) then
      if (eval(ist1).lt.esccut) nlorb=nlorb+1
    end if
  end if
end do
nlx=0
! open the atomic data file
open(50,file=trim(spsymb)//'.in',action='WRITE',form='FORMATTED')
write(50,'(" ''",A,"''",T45,": spsymb")') trim(spsymb)
write(50,'(" ''",A,"''",T45,": spname")') trim(spname)
write(50,'(G14.6,T45,": spzn")') spzn
write(50,'(G18.10,T45,": spmass")') spmass
write(50,'(G14.6,2F10.4,I6,T45,": sprmin, rmt, sprmax, nrmt")') sprmin,rmt, &
 sprmax,nrmt
write(50,'(I4,T45,": spnst")') spnst
write(50,'(3I4,G14.6,L1,T45,": spn, spl, spk, spocc, spcore")') spn(1),spl(1), &
 spk(1),spocc(1),spcore(1)
do ist1=2,spnst
  write(50,'(3I4,G14.6,L1)') spn(ist1),spl(ist1),spk(ist1),spocc(ist1), &
   spcore(ist1)
end do
write(50,'(I4,T45,": apword")') 1
write(50,'(F8.4,I4,"  ",L1,T45,": apwe0, apwdm, apwve")') boe,0,.false.
write(50,'(I4,T45,": nlx")') nlx
write(50,'(I4,T45,": nlorb")') nlorb
! write the local-orbitals
do l=0,3
  write(50,'(2I4,T45,": lorbl, lorbord")') l,2
  write(50,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') boe,0,.false.
  write(50,'(F8.4,I4,"  ",L1)') boe,1,.false.
end do
do ist1=1,spnst
  if (.not.spcore(ist1)) then
    if ((spl(ist1).eq.0).or.(spl(ist1).eq.spk(ist1))) then
      if (eval(ist1).lt.esccut) then
        write(50,'(2I4,T45,": lorbl, lorbord")') spl(ist1),3
        write(50,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') boe,0, &
         .false.
        write(50,'(F8.4,I4,"  ",L1)') boe,1,.false.
        write(50,'(F8.4,I4,"  ",L1)') eval(ist1)+boe,0,.true.
      end if
    end if
  end if
end do
close(50)
! read another element from file
goto 10
end program
