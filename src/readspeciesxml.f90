! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine readspeciesxml
use modmain
use modinput,only: input
use modsp
use  FoX_dom
implicit none
! local variables
integer::is, ist, iostat
integer::io, nlx, ilx, lx, ilo
type(Node), pointer :: speciesnp ,speciesdbnp

allocate(speziesdeflist(nspecies))
config => newDOMConfig()
parseerror=.false.
! parse xml and create derived type for species definitions  speziesdeflist
do is=1, nspecies
  if(trim(input%structure%speciesarray(is)%species%speciesfile).eq."")then
    input%structure%speciesarray(is)%species%speciesfile = &
    trim (input%structure%speciesarray(is)%species%chemicalSymbol)//".xml"
  endif
    write(*,*)"open ", trim(input%structure%speciespath)//"/"//trim(input%structure%speciesarray(is)%species%speciesfile)
  doc => parseFile(trim(input%structure%speciespath)//"/"//&
  trim(input%structure%speciesarray(is)%species%speciesfile),config)
  speciesdbnp=>getDocumentElement(doc)
  speciesnp=>item(getElementsByTagname(speciesdbnp,"sp"),0)
  parseerror=.false.
  speziesdeflist(is)%sp=>getstructsp(speciesnp)
  call destroy(doc)
end do

do is=1, nspecies
   spsymb(is)=speziesdeflist(is)%sp%chemicalSymbol
   input%structure%speciesarray(is)%species%chemicalSymbol= spsymb(is)
   spname(is)=speziesdeflist(is)%sp%name
   spzn(is)=speziesdeflist(is)%sp%z
   spmass(is)=speziesdeflist(is)%sp%mass
   sprmin(is)=speziesdeflist(is)%sp%muffinTin%rmin
   rmt(is)=speziesdeflist(is)%sp%muffinTin%radius
   if (input%structure%speciesarray(is)%species%rmt .gt. 0) then
    rmt(is)=input%structure%speciesarray(is)%species%rmt
   end if
   sprmax(is)=speziesdeflist(is)%sp%muffinTin%rinf
   nrmt(is)=speziesdeflist(is)%sp%muffinTin%radialmeshPoints
   if (sprmin(is).le.0.d0) then
     write(*, *)
     write(*, '("Error(readinput): sprmin <= 0 : ", G18.10)') sprmin(is)
     write(*, '(" for species ", I4)') is
     write(*, *)
     stop
   end if
   if (rmt(is).le.sprmin(is)) then
     write(*, *)
    write(*, '("Error(readinput): rmt <= sprmin : ", 2G18.10)') rmt(is), sprmin(is)
    write(*, '(" for species ", I4)') is
    write(*, *)
    stop
  end if
  if (sprmax(is).lt.rmt(is)) then
    write(*, *)
    write(*, '("Error(readinput): sprmax < rmt : ", 2G18.10)') sprmax(is), rmt(is)
    write(*, *)
    stop
  end if
  if (nrmt(is).lt.20) then
    write(*, *)
    write(*, '("Error(readinput): nrmt too small : ", I8)') nrmt(is)
    write(*, '(" for species ", I4)') is
    write(*, *)
    stop
  end if
  spnst(is)=size(speziesdeflist(is)%sp%AtomicStatearray)
  if (spnst(is).le.0) then
    write(*, *)
    write(*, '("Error(readinput): invalid spnst : ", I8)') spnst(is)
    write(*, '(" for species ", I4)') is
    write(*, *)
    stop
  end if
  if (spnst(is).gt.maxspst) then
    write(*, *)
    write(*, '("Error(readinput): too many states for species ", I8)') is
    write(*, *)
    stop
  end if
  do ist=1, spnst(is)

    spn(ist, is)=speziesdeflist(is)%sp%atomicStatearray(ist)%atomicState%n
    spl(ist, is)=speziesdeflist(is)%sp%atomicStatearray(ist)%atomicState%l
    spk(ist, is)=speziesdeflist(is)%sp%atomicStatearray(ist)%atomicState%kappa
    spocc(ist, is)=speziesdeflist(is)%sp%atomicStatearray(ist)%atomicState%occ
    spcore(ist, is)=speziesdeflist(is)%sp%atomicStatearray(ist)%atomicState%core
    if (spn(ist, is).lt.1) then
      write(*, *)
      write(*, '("Error(readinput): spn < 1 : ", I8)') spn(ist, is)
      write(*, '(" for species ", I4)') is
      write(*, '(" and state ", I4)') ist
      write(*, *)
      stop
    end if
    if (spl(ist, is).lt.0) then
      write(*, *)
      write(*, '("Error(readinput): spl < 0 : ", I8)') spl(ist, is)
      write(*, '(" for species ", I4)') is
      write(*, '(" and state ", I4)') ist
      write(*, *)
      stop
    end if
    if (spk(ist, is).lt.1) then
      write(*, *)
      write(*, '("Error(readinput): spk < 1 : ", I8)') spk(ist, is)
      write(*, '(" for species ", I4)') is
      write(*, '(" and state ", I4)') ist
      write(*, *)
      stop
    end if
    if (spocc(ist, is).lt.0.d0) then
      write(*, *)
      write(*, '("Error(readinput): spocc < 0 : ", G18.10)') spocc(ist, is)
      write(*, '(" for species ", I4)') is
      write(*, '(" and state ", I4)') ist
      write(*, *)
      stop
    end if
  end do
  apword(0, is)=speziesdeflist(is)%sp%basis%order
  if (apword(0, is).le.0) then
    write(*, *)
    write(*, '("Error(readinput): apword <= 0 : ", I8)') apword(0, is)
    write(*, '(" for species ", I4)') is
    write(*, *)
    stop
  end if
  if (apword(0, is).gt.maxapword) then
    write(*, *)
    write(*, '("Error(readinput): apword too large : ", I8)') apword(0, is)
    write(*, '(" for species ", I4)') is
    write(*, '("Adjust maxapword in modmain and recompile code")')
    write(*, *)
    stop
  end if
! set the APW orders for l>0
  apword(1:input%groundstate%lmaxapw, is)=apword(0, is)
  do io=1, apword(0, is)
    apwe0(io, 0, is)=speziesdeflist(is)%sp%basis%wfarray(io)%wf%trialEnergy
    apwdm(io, 0, is)=speziesdeflist(is)%sp%basis%wfarray(io)%wf%matchingOrder
    apwve(io, 0, is)=speziesdeflist(is)%sp%basis%wfarray(io)%wf%searchE
    if (apwdm(io, 0, is).lt.0) then
      write(*, *)
      write(*, '("Error(readinput): apwdm < 0 : ", I8)') apwdm(io, 0, is)
      write(*, '(" for species ", I4)') is
      write(*, '(" and order ", I4)') io
      write(*, *)
      stop
    end if
! set the APW linearisation energies, derivative orders and variability for l>0
    apwe0(io, 1:input%groundstate%lmaxapw, is)=apwe0(io, 0, is)
    apwdm(io, 1:input%groundstate%lmaxapw, is)=apwdm(io, 0, is)
    apwve(io, 1:input%groundstate%lmaxapw, is)=apwve(io, 0, is)
  end do

  nlx=size(speziesdeflist(is)%sp%basis%exceptionarray)
  if (nlx.lt.0 .and. associated(speziesdeflist(is)%sp%basis%exceptionarray)) then
    write(*, *)
    write(*, '("Error(readinput): nlx < 0 : ", I8)') nlx
    write(*, '(" for species ", I4)') is
    write(*, *)
    stop
  end if
  do ilx=1, nlx
     lx=speziesdeflist(is)%sp%basis%exceptionarray(ilx)%exception%l
     io=size(speziesdeflist(is)%sp%basis%exceptionarray(ilx)%exception%wfarray)
    if (lx.lt.0) then
      write(*, *)
      write(*, '("Error(readinput): lx < 0 : ", I8)') lx
      write(*, '(" for species ", I4)') is
      write(*, '(" and exception number ", I4)') ilx
      write(*, *)
      stop
    end if
    if (lx.gt.input%groundstate%lmaxapw) then
      write(*, *)
      write(*, '("Error(readinput): lx > lmaxapw : ", I8)') lx
      write(*, '(" for species ", I4)') is
      write(*, '(" and exception number ", I4)') ilx
      write(*, *)
      stop
    end if
    apword(lx, is)=io
    if (apword(lx, is).le.0) then
      write(*, *)
      write(*, '("Error(readinput): apword <= 0 : ", I8)') apword(lx, is)
      write(*, '(" for species ", I4)') is
      write(*, '(" and exception number ", I4)') ilx
      write(*, *)
      stop
    end if
    if (apword(lx, is).gt.maxapword) then
      write(*, *)
      write(*, '("Error(readinput): apword too large : ", I8)') apword(lx, is)
      write(*, '(" for species ", I4)') is
      write(*, '(" and exception number ", I4)') ilx
      write(*, '("Adjust maxapword in modmain and recompile code")')
      write(*, *)
      stop
    end if
    do io=1, apword(lx, is)

       apwe0(io, lx, is)=speziesdeflist(is)%sp%basis%exceptionarray(ilx)%exception%wfarray(io)%wf%trialEnergy
       apwdm(io, lx, is)=speziesdeflist(is)%sp%basis%exceptionarray(ilx)%exception%wfarray(io)%wf%matchingOrder
       apwve(io, lx, is)=speziesdeflist(is)%sp%basis%exceptionarray(ilx)%exception%wfarray(io)%wf%searchE
      if (apwdm(io, lx, is).lt.0) then
	write(*, *)
	write(*, '("Error(readinput): apwdm < 0 : ", I8)') apwdm(io, lx, is)
	write(*, '(" for species ", I4)') is
	write(*, '(" exception number ", I4)') ilx
	write(*, '(" and order ", I4)') io
	write(*, *)
	stop
      end if
    end do
  end do
 nlorb(is)=size(speziesdeflist(is)%sp%lorbarray)
  if (nlorb(is).lt.0) then
    write(*, *)
    write(*, '("Error(readinput): nlorb < 0 : ", I8)') nlorb(is)
    write(*, '(" for species ", I4)') is
    write(*, *)
    stop
  end if
  if (nlorb(is).gt.maxlorb) then
    write(*, *)
    write(*, '("Error(readinput): nlorb too large : ", I8)') nlorb(is)
    write(*, '(" for species ", I4)') is
    write(*, '("Adjust maxlorb in modmain and recompile code")')
    write(*, *)
    stop
  end if
  do ilo=1, nlorb(is)
    lorbl(ilo, is)=speziesdeflist(is)%sp%lorbarray(ilo)%lorb%l
    lorbord(ilo, is)=size(speziesdeflist(is)%sp%lorbarray(ilo)%lorb%wfarray)
    if (lorbl(ilo, is).lt.0) then
      write(*, *)
      write(*, '("Error(readinput): lorbl < 0 : ", I8)') lorbl(ilo, is)
      write(*, '(" for species ", I4)') is
      write(*, '(" and local-orbital ", I4)') ilo
      write(*, *)
      stop
    end if
    if (lorbl(ilo, is).gt.input%groundstate%lmaxmat) then
      write(*, *)
      write( * , '("Error(readinput): lorbl > lmaxmat : ", 2I8)') lorbl(ilo, is), &
       input%groundstate%lmaxmat
      write(*, '(" for species ", I4)') is
      write(*, '(" and local-orbital ", I4)') ilo
      write(*, *)
      stop
    end if
    if (lorbord(ilo, is).lt.2) then
      write(*, *)
      write(*, '("Error(readinput): lorbord < 2 : ", I8)') lorbord(ilo, is)
      write(*, '(" for species ", I4)') is
      write(*, '(" and local-orbital ", I4)') ilo
      write(*, *)
      stop
    end if
    if (lorbord(ilo, is).gt.maxlorbord) then
      write(*, *)
      write(*, '("Error(readinput): lorbord too large : ", I8)') lorbord(ilo, is)
      write(*, '(" for species ", I4)') is
      write(*, '(" and local-orbital ", I4)') ilo
      write(*, '("Adjust maxlorbord in modmain and recompile code")')
      write(*, *)
      stop
    end if
    do io=1, lorbord(ilo, is)
     lorbe0(io, ilo, is)=speziesdeflist(is)%sp%lorbarray(ilo)%lorb%wfarray(io)%wf%trialEnergy
     lorbdm(io, ilo, is)=speziesdeflist(is)%sp%lorbarray(ilo)%lorb%wfarray(io)%wf%matchingOrder
     lorbve(io, ilo, is)=speziesdeflist(is)%sp%lorbarray(ilo)%lorb%wfarray(io)%wf%searchE
      if (lorbdm(io, ilo, is).lt.0) then
	write(*, *)
	write(*, '("Error(readinput): lorbdm < 0 : ", I8)') lorbdm(io, ilo, is)
	write(*, '(" for species ", I4)') is
	write(*, '(" local-orbital ", I4)') ilo
	write(*, '(" and order ", I4)') io
	write(*, *)
	stop
      end if
    end do
  end do
  close(50)
end do
return
end subroutine
