
! Copyright (C) 2005-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine rdmenergy
! calculate the total energy for RDMFT
use modmain
implicit none
! local variables
integer is,ia,ias
integer ik,ist,ir
real(8) vn,t1
complex(8) zt1
! allocatable arrays
real(8), allocatable :: rfmt(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: c(:,:)
! external functions
real(8) rfmtinp
complex(8) zdotc
external rfmtinp,zdotc
allocate(rfmt(lmmaxvr,nrmtmax))
allocate(evecsv(nstsv,nstsv))
allocate(c(nstsv,nstsv))
! Coulomb energy from core states
engyvcl=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      rfmt(1,ir)=rhocr(ir,ias)/y00
    end do
    engyvcl=engyvcl+rfmtinp(1,0,nrmt(is),spr(:,is),lmmaxvr,rfmt,vclmt(:,:,ias))
  end do
end do
engykn=engykncr
do ik=1,nkpt
  call getevecsv(vkl(:,ik),evecsv)
  do ist=1,nstsv
    t1=wkpt(ik)*occsv(ist,ik)
! Coulomb energy from valence states
    engyvcl=engyvcl+t1*dble(vclmat(ist,ist,ik))
! kinetic energy from valence states
    zt1=zdotc(nstsv,evecsv(:,ist),1,dkdc(:,ist,ik),1)
    engykn=engykn+t1*dble(zt1)
  end do
end do
! Madelung term
engymad=0.d0
do is=1,nspecies
! compute the bare nucleus potential at the origin
  call potnucl(ptnucl,1,spr(:,is),spzn(is),vn)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    engymad=engymad+0.5d0*spzn(is)*(vclmt(1,1,ias)*y00-vn)
  end do
end do
! exchange-correlation energy
call rdmengyxc
! total energy
engytot=0.5d0*engyvcl+engymad+engykn+engyx
if (rdmtemp.gt.0.d0) then
  call rdmentropy
  engytot=engytot-rdmtemp*rdmentrpy
end if
write(*,*)
write(*,'("Info(rdmengy): Energies :")')
write(*,'(" electronic kinetic",T30,": ",G18.10)') engykn
write(*,'(" core electron kinetic",T30,": ",G18.10)') engykncr
write(*,'(" Coulomb",T30,": ",G18.10)') engyvcl
write(*,'(" Madelung",T30,": ",G18.10)') engymad
write(*,'(" exchange-correlation",T30,": ",G18.10)') engyx
if (rdmtemp.gt.0.d0) then
  write(*,'(" entropy",T30,": ",G18.10)') rdmtemp*rdmentrpy
end if
write(*,'(" total",T30,": ",G18.10)') engytot
write(*,*)
deallocate(evecsv,rfmt,c)
return
end subroutine

