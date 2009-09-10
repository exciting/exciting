
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

program species
use modspecies
use  mod_muffin_tin , only:epsedirac,epspotatom
implicit none
!
! Parameterset corresponding to original source code for species-program.
! Note that for some of the shipped species files the cutoff parameters are
! different (epsedirac=1.d-10 and epspotatom=eps=1.d-7)
epsedirac=1.d-11
epspotatom=1.d-6
! set default values
inspecies=''
ecvcut=-3.5d0
esccut=-0.35d0
apwdescr=''
suffix=''
apword=1
apwdm(:)=0
apwve(:)=.false.
apwordx=0
apwdmx(:)=0
apwve(:)=.false.
locorb=.true.
locorbsc=.true.
searchlocorb=.false.
fullsearchlocorbsc=.false.

! get generation strategy from input file 'species.input'
open(50,file='species.input',action='READ',status='OLD',form='FORMATTED', &
 iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(species): error opening species.input")')
  write(*,*)
  stop
end if
100 continue
read(50,*,end=300) bname
! check for a comment
if ((scan(trim(bname),'!').eq.1).or.(scan(trim(bname),'#').eq.1)) goto 100
select case(trim(bname))
case('epsedirac')
  read(50,*,err=200) epsedirac
  if (epsedirac.lt.0.d0) then
    write(*,*)
    write(*,'("Error(species): epsedirac < 0 : ",g18.10)') epsedirac
    write(*,*)
    stop
  end if
case('epspotatom')
  read(50,*,err=200) epspotatom
  if (epspotatom.lt.0.d0) then
    write(*,*)
    write(*,'("Error(species): epspotatom < 0 : ",g18.10)') epspotatom
    write(*,*)
    stop
  end if
case('inspecies')
  read(50,*,err=200) inspecies
case('ecvcut')
  read(50,*,err=200) ecvcut
case('esccut')
  read(50,*,err=200) esccut
case('apwdescr')
  read(50,*,err=200) apwdescr
case('suffix')
  read(50,*,err=200) suffix
case('apw')
  read(50,*,err=200) apword
  if ((apword.lt.1).or.(apword.gt.maxapword)) then
    write(*,*)
    write(*,'("Error(species): apword < 1 : ",A)') apword
    write(*,*)
    stop
  end if
  if (apword.gt.maxapword) then
    write(*,*)
    write(*,'("Error(species): apword too large : ",I8)') apword
    write(*,*)
    stop
  end if
  do io=1,apword
    read(50,*) apwdm(io),apwve(io)
    if (apwdm(io).lt.0) then
      write(*,*)
      write(*,'("Error(species): apwdm < 0 : ",I8)') apwdm(io)
      write(*,'(" for order ",I4)') io
      write(*,*)
      stop
    end if
  end do
case('apwx')
  read(50,*,err=200) apwordx
  if ((apwordx.lt.1).or.(apwordx.gt.maxapword)) then
    write(*,*)
    write(*,'("Error(species): apwordx < 1 : ",A)') apwordx
    write(*,*)
    stop
  end if
  if (apwordx.gt.maxapword) then
    write(*,*)
    write(*,'("Error(species): apwordx too large : ",I8)') apwordx
    write(*,*)
    stop
  end if
  do io=1,apwordx
    read(50,*) apwdmx(io),apwvex(io)
    if (apwdmx(io).lt.0) then
      write(*,*)
      write(*,'("Error(species): apwdmx < 0 : ",I8)') apwdmx(io)
      write(*,'(" for order ",I4)') io
      write(*,*)
      stop
    end if
  end do
case('locorb')
  read(50,*,err=200) locorb
case('locorbsc')
  read(50,*,err=200) locorbsc
case('searchlocorb')
  read(50,*,err=200) searchlocorb
case('fullsearchlocorbsc')
  read(50,*,err=200) fullsearchlocorbsc
case('')
  goto 100
case default
  write(*,*)
  write(*,'("Error(species): invalid block name : ",A)') trim(bname)
  write(*,*)
  stop
end select
goto 100
200 continue
write(*,*)
write(*,'("Error(species): error reading from species.in")')
write(*,'("Problem occurred in ''",A,"'' block")') trim(bname)
write(*,*)
stop
300 continue
close(50)



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
do ist=1,spnst
  read(40,*) spn(ist),spl(ist),spk(ist),i
  if (ist.ge.2) then
    if (spn(ist).lt.spn(ist-1)) then
      write(*,*)
      write(*,'("Error(species): states improperly ordered")')
      write(*,'(" for species ",A)') trim(spname)
      write(*,*)
      stop
    end if
  end if
  spocc(ist)=dble(i)
end do
read(40,*)
if ((adjustl(trim(inspecies)).ne.'').and.(adjustl(trim(inspecies)).ne. &
	adjustl(trim(spsymb)))) goto 10
write(*,'("Info(species): method: ",a,"; running Z = ",I4,", (",A,")")') trim(apwdescr),nz, &
  trim(spname)
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
  call atom(.true.,spzn,spnst,spn,spl,spk,spocc,xctype,xcgrad,np,spnr,r,eval, &
   rho,vr,rwf)
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
do ist=1,spnst
  if (eval(ist).lt.ecvcut) then
    spcore(ist)=.true.
  else
    spcore(ist)=.false.
  end if
end do
! check that the state for same n and l but different k is also core
do ist=1,spnst
  if (spcore(ist)) then
    do jst=1,spnst
      if ((spn(ist).eq.spn(jst)).and.(spl(ist).eq.spl(jst))) spcore(jst)=.true.
    end do
  end if
end do
! find the total number of local orbitals
nlorb=0
maxl=0
do ist=1,spnst
  if (.not.spcore(ist)) then
    if ((spl(ist).eq.0).or.(spl(ist).eq.spk(ist))) then
      l=spl(ist)
      if (eval(ist).lt.esccut) nlorb=nlorb+1
    end if
    if (spl(ist).gt.maxl) maxl=spl(ist)
  end if
end do
! save number of semi-core LO's
nlorbsc=nlorb
maxl=maxl+1
if (maxl.gt.3) maxl=3
! case for no local orbitals
if ((.not.locorb).and.(.not.locorbsc)) nlorb=0
! case for small lo's and semi-core lo's
if (locorb.and.locorbsc) nlorb=nlorbsc+maxl+1
! case for only semi-core LO's
if ((.not.locorb).and.locorbsc) nlorb=nlorbsc
! case for only small lo's
if (locorb.and.(.not.locorbsc)) nlorb=maxl+1

! open the atomic data file
open(50,file=trim(spsymb)//trim(suffix)//'.in',action='WRITE',form='FORMATTED')
write(50,'(" ''",A,"''",T45,": spsymb")') trim(spsymb)
write(50,'(" ''",A,"''",T45,": spname")') trim(spname)
write(50,'(G14.6,T45,": spzn")') spzn
write(50,'(G18.10,T45,": spmass")') spmass
write(50,'(G14.6,2F10.4,I6,T45,": sprmin, rmt, sprmax, nrmt")') sprmin,rmt, &
 sprmax,nrmt
write(50,'(I4,T45,": spnst")') spnst
write(50,'(3I4,G14.6,L1,T45,": spn, spl, spk, spocc, spcore")') spn(1),spl(1), &
 spk(1),spocc(1),spcore(1)
do ist=2,spnst
  write(50,'(3I4,G14.6,L1)') spn(ist),spl(ist),spk(ist),spocc(ist),spcore(ist)
end do


! overall APW order
write(50,'(I4,T45,": apword")') apword
do io=1,apword
   write(50,'(F8.4,I4,"  ",L1,T45,": apwe0, apwdm, apwve")') boe,apwdm(io),apwve(io)
end do



if (apwordx.gt.0) then
  ! number of exceptions corresponds to number of l-values
  nlx=maxl+1
  write(50,'(I4,T45,": nlx")') nlx
  ! write the exceptions
  do l=0,maxl
    write(50,'(2I4,T45,": lorbl, lorbord")') l,apwordx
    do io=1,apwordx
      if (io.eq.1) then
        write(50,'(F8.4,I4,"  ",L1,T45,": apwe0, apwdm, apwve")') boe,apwdmx(io),apwvex(io)
      else
        write(50,'(F8.4,I4,"  ",L1,T45)') boe,apwdmx(io),apwvex(io)
      end if
    end do
  end do
else
  nlx=0
  write(50,'(I4,T45,": nlx")') nlx
end if

write(50,'(I4,T45,": nlorb")') nlorb

if (locorb) then
  ! write the local-orbitals
  do l=0,maxl
    write(50,'(2I4,T45,": lorbl, lorbord")') l,2
    if (searchlocorb) then
       write(50,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') boe,0,.true.
       write(50,'(F8.4,I4,"  ",L1)') boe,1,.true.
    else
       write(50,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') boe,0,.false.
       write(50,'(F8.4,I4,"  ",L1)') boe,1,.false.
    end if
  end do
end if

if (locorbsc) then
  do ist=1,spnst
    if (.not.spcore(ist)) then
      if ((spl(ist).eq.0).or.(spl(ist).eq.spk(ist))) then
        if (eval(ist).lt.esccut) then
          write(50,'(2I4,T45,": lorbl, lorbord")') spl(ist),3
          write(50,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') boe,0, &
           fullsearchlocorbsc
          write(50,'(F8.4,I4,"  ",L1)') boe,1,fullsearchlocorbsc
          write(50,'(F8.4,I4,"  ",L1)') eval(ist)+0.5d0*boe,0,.true.
        end if
      end if
    end if
  end do
end if

write(50,*)
write(50,'("# Exciting code version : ",a)') version
write(50,'("# Description of method : ",a)') trim(apwdescr)

close(50)
 call writexmlspecies()
! read another element from file
goto 10
end program
