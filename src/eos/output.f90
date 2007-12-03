
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine output
use modmain
implicit none
! local variables
integer ip,ipt,iplt
real(8) v
! external functions
real(8) eveos,pveos
external eveos,pveos
! output parameters
open(60,file='PARAM.OUT')
write(60,*)
write(60,'(A)') trim(cname)
write(60,*)
write(60,'(A)') trim(ename(1))
write(60,'(A)') trim(ename(2))
write(60,*)
write(60,'("(Default units are atomic: Hartree, Bohr etc.) ")')
write(60,*)
do ip=1,nparam
  write(60,'(" ",A,T20,"=",T30,G18.10)') trim(pname(ip)),popt(ip)
end do
write(60,*)
do ip=1,nparam
  if (trim(pname(ip)).eq."B0") then
    write(60,'(" B0 (GPa)",T20,"=",T30,G18.10)') popt(ip)*aupress_gpa
  end if
  if (trim(pname(ip)).eq."B0''") then
    write(60,'(A4," (/GPa)",T20,"=",T30,G18.10)') "B0''",popt(ip)/aupress_gpa
  end if
end do
write(60,*)
close(60)
! output energy vs volume per atom at data points
open(60,file='EVPAP.OUT')
do ipt=1,nevpt
  write(60,*) vpt(ipt)/dble(natoms),ept(ipt)/dble(natoms)
end do
close(60)
! output energy vs volume per atom over volume interval
open(60,file='EVPAI.OUT')
do iplt=1,nvplt
  v=(vplt2-vplt1)*dble(iplt)/dble(nvplt)+vplt1
  write(60,*) v/dble(natoms),eveos(etype,popt,v)/dble(natoms)
end do
close(60)
! output pressure vs volume per atom at data points
open(60,file='PVPAP.OUT')
do ipt=1,nevpt
  write(60,*) vpt(ipt)/dble(natoms),pveos(etype,popt,vpt(ipt))*aupress_gpa
end do
close(60)
! output pressure vs volume per atom over volume interval
open(60,file='PVPAI.OUT')
do iplt=1,nvplt
  v=(vplt2-vplt1)*dble(iplt)/dble(nvplt)+vplt1
  write(60,*) v/dble(natoms),pveos(etype,popt,v)*aupress_gpa
end do
close(60)
! output enthalpy vs pressure per atom over volume interval
open(60,file='HPPAI.OUT')
do iplt=1,nvplt
  v=(vplt2-vplt1)*dble(iplt)/dble(nvplt)+vplt1
  write(60,*) pveos(etype,popt,v)*aupress_gpa, &
   (eveos(etype,popt,v)+pveos(etype,popt,v)*v)/dble(natoms)
end do
close(60)
write(*,*)
write(*,'("All units are atomic unless otherwise stated")')
write(*,'("EOS parameters written to PARAM.OUT")')
write(*,'("Energy-volume per atom at data points written to EVPAP.OUT")')
write(*,'("Energy-volume per atom over interval written to EVPAI.OUT")')
write(*,'("Pressure(GPa)-volume per atom at data points written to PVPAP.OUT")')
write(*,'("Pressure(GPa)-volume per atom over interval written to PVPAI.OUT")')
write(*,'("Enthalpy-pressure(GPa) per atom over interval written to HPPAI.OUT")')
write(*,*)
return
end subroutine
