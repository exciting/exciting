
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readspecies
use modmain
implicit none
! local variables
integer is,ist,iostat
integer io,nlx,ilx,lx,ilo
do is=1,nspecies
  open(50,file=trim(sppath)//trim(spfname(is)),action='READ',status='OLD', &
   form='FORMATTED',iostat=iostat)
  if (iostat.ne.0) then
    write(*,*)
    write(*,'("Error(readinput): error opening species file ",A)') &
     trim(sppath)//trim(spfname(is))
    write(*,*)
    stop
  end if
  read(50,*) spsymb(is)
  read(50,*) spname(is)
  read(50,*) spzn(is)
  read(50,*) spmass(is)
  read(50,*) sprmin(is),rmt(is),sprmax(is),nrmt(is)
  if (sprmin(is).le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): sprmin <= 0 : ",G18.10)') sprmin(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (rmt(is).le.sprmin(is)) then
    write(*,*)
    write(*,'("Error(readinput): rmt <= sprmin : ",2G18.10)') rmt(is),sprmin(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (sprmax(is).lt.rmt(is)) then
    write(*,*)
    write(*,'("Error(readinput): sprmax < rmt : ",2G18.10)') sprmax(is),rmt(is)
    write(*,*)
    stop
  end if
  if (nrmt(is).lt.20) then
    write(*,*)
    write(*,'("Error(readinput): nrmt too small : ",I8)') nrmt(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  read(50,*) spnst(is)
  if (spnst(is).le.0) then
    write(*,*)
    write(*,'("Error(readinput): invalid spnst : ",I8)') spnst(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (spnst(is).gt.maxspst) then
    write(*,*)
    write(*,'("Error(readinput): too many states for species ",I8)') is
    write(*,*)
    stop
  end if
  do ist=1,spnst(is)
    read(50,*) spn(ist,is),spl(ist,is),spk(ist,is),spocc(ist,is),spcore(ist,is)
    if (spn(ist,is).lt.1) then
      write(*,*)
      write(*,'("Error(readinput): spn < 1 : ",I8)') spn(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
    if (spl(ist,is).lt.0) then
      write(*,*)
      write(*,'("Error(readinput): spl < 0 : ",I8)') spl(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
    if (spk(ist,is).lt.1) then
      write(*,*)
      write(*,'("Error(readinput): spk < 1 : ",I8)') spk(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
    if (spocc(ist,is).lt.0.d0) then
      write(*,*)
      write(*,'("Error(readinput): spocc < 0 : ",G18.10)') spocc(ist,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and state ",I4)') ist
      write(*,*)
      stop
    end if
  end do
  read(50,*) apword(0,is)
  if (apword(0,is).le.0) then
    write(*,*)
    write(*,'("Error(readinput): apword <= 0 : ",I8)') apword(0,is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (apword(0,is).gt.maxapword) then
    write(*,*)
    write(*,'("Error(readinput): apword too large : ",I8)') apword(0,is)
    write(*,'(" for species ",I4)') is
    write(*,'("Adjust maxapword in modmain and recompile code")')
    write(*,*)
    stop
  end if
! set the APW orders for l>0
  apword(1:lmaxapw,is)=apword(0,is)
  do io=1,apword(0,is)
    read(50,*) apwe0(io,0,is),apwdm(io,0,is),apwve(io,0,is)
    if (apwdm(io,0,is).lt.0) then
      write(*,*)
      write(*,'("Error(readinput): apwdm < 0 : ",I8)') apwdm(io,0,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and order ",I4)') io
      write(*,*)
      stop
    end if
! set the APW linearisation energies, derivative orders and variability for l>0
    apwe0(io,1:lmaxapw,is)=apwe0(io,0,is)
    apwdm(io,1:lmaxapw,is)=apwdm(io,0,is)
    apwve(io,1:lmaxapw,is)=apwve(io,0,is)
  end do
  read(50,*) nlx
  if (nlx.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nlx < 0 : ",I8)') nlx
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  do ilx=1,nlx
    read(50,*) lx,io
    if (lx.lt.0) then
      write(*,*)
      write(*,'("Error(readinput): lx < 0 : ",I8)') lx
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    end if
    if (lx.gt.lmaxapw) then
      write(*,*)
      write(*,'("Error(readinput): lx > lmaxapw : ",I8)') lx
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    end if
    apword(lx,is)=io
    if (apword(lx,is).le.0) then
      write(*,*)
      write(*,'("Error(readinput): apword <= 0 : ",I8)') apword(lx,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,*)
      stop
    end if
    if (apword(lx,is).gt.maxapword) then
      write(*,*)
      write(*,'("Error(readinput): apword too large : ",I8)') apword(lx,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and exception number ",I4)') ilx
      write(*,'("Adjust maxapword in modmain and recompile code")')
      write(*,*)
      stop
    end if
    do io=1,apword(lx,is)
      read(50,*) apwe0(io,lx,is),apwdm(io,lx,is),apwve(io,lx,is)
      if (apwdm(io,lx,is).lt.0) then
        write(*,*)
        write(*,'("Error(readinput): apwdm < 0 : ",I8)') apwdm(io,lx,is)
        write(*,'(" for species ",I4)') is
        write(*,'(" exception number ",I4)') ilx
        write(*,'(" and order ",I4)') io
        write(*,*)
        stop
      end if
    end do
  end do
  read(50,*) nlorb(is)
  if (nlorb(is).lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nlorb < 0 : ",I8)') nlorb(is)
    write(*,'(" for species ",I4)') is
    write(*,*)
    stop
  end if
  if (nlorb(is).gt.maxlorb) then
    write(*,*)
    write(*,'("Error(readinput): nlorb too large : ",I8)') nlorb(is)
    write(*,'(" for species ",I4)') is
    write(*,'("Adjust maxlorb in modmain and recompile code")')
    write(*,*)
    stop
  end if
  do ilo=1,nlorb(is)
    read(50,*) lorbl(ilo,is),lorbord(ilo,is)
    if (lorbl(ilo,is).lt.0) then
      write(*,*)
      write(*,'("Error(readinput): lorbl < 0 : ",I8)') lorbl(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    end if
    if (lorbl(ilo,is).gt.lmaxmat) then
      write(*,*)
      write(*,'("Error(readinput): lorbl > lmaxmat : ",2I8)') lorbl(ilo,is), &
       lmaxmat
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    end if
    if (lorbord(ilo,is).lt.2) then
      write(*,*)
      write(*,'("Error(readinput): lorbord < 2 : ",I8)') lorbord(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,*)
      stop
    end if
    if (lorbord(ilo,is).gt.maxlorbord) then
      write(*,*)
      write(*,'("Error(readinput): lorbord too large : ",I8)') lorbord(ilo,is)
      write(*,'(" for species ",I4)') is
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,'("Adjust maxlorbord in modmain and recompile code")')
      write(*,*)
      stop
    end if
    do io=1,lorbord(ilo,is)
      read(50,*) lorbe0(io,ilo,is),lorbdm(io,ilo,is),lorbve(io,ilo,is)
      if (lorbdm(io,ilo,is).lt.0) then
        write(*,*)
        write(*,'("Error(readinput): lorbdm < 0 : ",I8)') lorbdm(io,ilo,is)
        write(*,'(" for species ",I4)') is
        write(*,'(" local-orbital ",I4)') ilo
        write(*,'(" and order ",I4)') io
        write(*,*)
        stop
      end if
    end do
  end do
  close(50)
end do
return
end subroutine

