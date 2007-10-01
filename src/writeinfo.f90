
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeinfo
! !INTERFACE:
subroutine writeinfo(fnum)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   fnum : unit specifier for INFO.OUT file (in,integer)
! !DESCRIPTION:
!   Outputs basic information about the run to the file {\tt INFO.OUT}. Does not
!   close the file afterwards.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer fnum
! local variables
integer i,is,ia
character(10) dat,tim
write(fnum,'("+----------------------------------+")')
write(fnum,'("| EXCITING version ",I1.1,".",I1.1,".",I3.3," started |")') &
 version
write(fnum,'("+----------------------------------+")')
if (notelns.gt.0) then
  write(fnum,*)
  write(fnum,'("Notes :")')
  do i=1,notelns
    write(fnum,'(A)') notes(i)
  end do
end if
call date_and_time(date=dat,time=tim)
write(fnum,*)
write(fnum,'("Date (YYYY-MM-DD) : ",A4,"-",A2,"-",A2)') dat(1:4),dat(5:6), &
 dat(7:8)
write(fnum,'("Time (hh:mm:ss)   : ",A2,":",A2,":",A2)') tim(1:2),tim(3:4), &
 tim(5:6)
select case(task)
case(0)
  write(fnum,*)
  write(fnum,'("+-------------------------------------------------+")')
  write(fnum,'("| Ground-state run starting from atomic densities |")')
  write(fnum,'("+-------------------------------------------------+")')
case(1)
  write(fnum,*)
  write(fnum,'("+------------------------------------------+")')
  write(fnum,'("| Ground-state run resuming from STATE.OUT |")')
  write(fnum,'("+------------------------------------------+")')
case(2)
  write(fnum,*)
  write(fnum,'("+--------------------------------------------------------+")')
  write(fnum,'("| Structural optimisation starting from atomic densities |")')
  write(fnum,'("+--------------------------------------------------------+")')
case(3)
  write(fnum,*)
  write(fnum,'("+-----------------------------------------------------+")')
  write(fnum,'("| Structural optimisation run resuming from STATE.OUT |")')
  write(fnum,'("+-----------------------------------------------------+")')
case(5)
  write(fnum,*)
  write(fnum,'("+-------------------------------+")')
  write(fnum,'("| Ground-state Hartree-Fock run |")')
  write(fnum,'("+-------------------------------+")')
case(300)
  write(fnum,*)
  write(fnum,'("+---------------------------------------+")')
  write(fnum,'("| Reduced density matrix functional run |")')
  write(fnum,'("+---------------------------------------+")')
case default
  write(*,*)
  write(*,'("Error(writeinfo): task not defined : ",I8)') task
  write(*,*)
  stop
end select
write(fnum,*)
write(fnum,'("Lattice vectors :")')
write(fnum,'(3G18.10)') avec(1,1),avec(2,1),avec(3,1)
write(fnum,'(3G18.10)') avec(1,2),avec(2,2),avec(3,2)
write(fnum,'(3G18.10)') avec(1,3),avec(2,3),avec(3,3)
write(fnum,*)
write(fnum,'("Reciprocal lattice vectors :")')
write(fnum,'(3G18.10)') bvec(1,1),bvec(2,1),bvec(3,1)
write(fnum,'(3G18.10)') bvec(1,2),bvec(2,2),bvec(3,2)
write(fnum,'(3G18.10)') bvec(1,3),bvec(2,3),bvec(3,3)
write(fnum,*)
write(fnum,'("Unit cell volume      : ",G18.10)') omega
write(fnum,'("Brillouin zone volume : ",G18.10)') (twopi**3)/omega
if (autormt) then
  write(fnum,*)
  write(fnum,'("Automatic determination of muffin-tin radii")')
  write(fnum,'(" parameters : ",2G18.10)') rmtapm
end if
do is=1,nspecies
  write(fnum,*)
  write(fnum,'("Species : ",I4," (",A,")")') is,trim(spsymb(is))
  write(fnum,'(" parameters loaded from : ",A)') trim(spfname(is))
  write(fnum,'(" name : ",A)') trim(spname(is))
  write(fnum,'(" nuclear charge    : ",G18.10)') spzn(is)
  write(fnum,'(" electronic charge : ",G18.10)') spze(is)
  write(fnum,'(" atomic mass : ",G18.10)') spmass(is)
  write(fnum,'(" muffin-tin radius : ",G18.10)') rmt(is)
  write(fnum,'(" number of radial points in muffin-tin : ",I6)') nrmt(is)
  write(fnum,'(" atomic position (lattice), magnetic field (Cartesian) :")')
  do ia=1,natoms(is)
    write(fnum,'(I4,3F12.8,"  ",3F12.8)') ia,atposl(:,ia,is),bfcmt(:,ia,is)
  end do
end do
write(fnum,*)
write(fnum,'("Total number of atoms per unit cell : ",I4)') natmtot
write(fnum,*)
write(fnum,'("Spin treatment :")')
if (spinpol) then
  write(fnum,'(" spin-polarised")')
else
  write(fnum,'(" spin-unpolarised")')
end if
if (spinorb) then
  write(fnum,'(" spin-orbit coupling")')
end if
if (fixspin) then
  write(fnum,'(" fixed spin moment (Cartesian) : ",3G18.10)') momfix
end if
if (spinpol) then
  write(fnum,'(" global magnetic field (Cartesian) : ",3G18.10)') bfieldc
end if
if (ndmag.eq.1) then
  write(fnum,'(" collinear magnetisation in z-direction")')
else if (ndmag.eq.3) then
  write(fnum,'(" non-collinear magnetisation")')
end if
if (spinsprl) then
  write(fnum,'(" spin-spiral state assumed")')
  write(fnum,'("  q-vector (lattice)   : ",3G18.10)') vqlss
  write(fnum,'("  q-vector (Cartesian) : ",3G18.10)') vqcss
  write(fnum,'("  q-vector length      : ",G18.10)') sqrt(vqcss(1)**2 &
   +vqcss(2)**2+vqcss(3)**2)
end if
write(fnum,*)
write(fnum,'("Number of Bravais lattice symmetries : ",I4)') nsymlat
write(fnum,'("Number of crystal symmetries         : ",I4)') nsymcrys
write(fnum,*)
if (autokpt) then
  write(fnum,'("de Broglie wavelength used to determine k-point grid : ",&
   &G18.10)') rlambda
end if
write(fnum,'("k-point grid : ",3I6)') ngridk
write(fnum,'("k-point offset : ",3G18.10)') vkloff
if (reducek) then
  write(fnum,'("k-point set is reduced with crystal symmetries")')
else
  write(fnum,'("k-point set is not reduced")')
end if
write(fnum,'("Total number of k-points : ",I8)') nkpt
write(fnum,*)
write(fnum,'("Smallest muffin-tin radius times maximum |G+k| : ",G18.10)') &
 rgkmax
write(fnum,'("Maximum |G+k| for APW functions       : ",G18.10)') gkmax
write(fnum,'("Maximum |G| for potential and density : ",G18.10)') gmaxvr
write(fnum,'("Polynomial order for pseudocharge density : ",I4)') npsden
write(fnum,*)
write(fnum,'("G-vector grid sizes : ",3I6)') ngrid(1),ngrid(2),ngrid(3)
write(fnum,'("Total number of G-vectors : ",I8)') ngvec
write(fnum,*)
write(fnum,'("Maximum angular momentum used for")')
write(fnum,'(" APW functions                     : ",I4)') lmaxapw
write(fnum,'(" computing H and O matrix elements : ",I4)') lmaxmat
write(fnum,'(" potential and density             : ",I4)') lmaxvr
write(fnum,'(" inner part of muffin-tin          : ",I4)') lmaxinr
write(fnum,*)
write(fnum,'("Total nuclear charge    : ",G18.10)') chgzn
write(fnum,'("Total core charge       : ",G18.10)') chgcr
write(fnum,'("Total valence charge    : ",G18.10)') chgval
write(fnum,'("Total excess charge     : ",G18.10)') chgexs
write(fnum,'("Total electronic charge : ",G18.10)') chgtot
write(fnum,*)
write(fnum,'("Number of empty states         : ",I4)') nempty
write(fnum,'("Total number of valence states : ",I4)') nstsv
write(fnum,*)
write(fnum,'("Total number of local-orbitals : ",I4)') nlotot
write(fnum,*)
if (task.eq.5) write(fnum,'("Hartree-Fock calculation using Kohn-Sham states")')
if (task.eq.300) then
  write(fnum,'("RDMFT calculation using Kohn-Sham states")')
  if (rdmxctype.eq.1) then
    write(fnum,'(" RDMFT exchange-correlation type : ",I4)') rdmxctype
    write(fnum,'(" Hartree-Fock functional")')
  else if (rdmxctype.eq.2) then
    write(fnum,'(" RDMFT exchange-correlation type : ",I4)') rdmxctype
    write(fnum,'(" Mueller functional")')
  else if (rdmxctype.eq.3) then
    write(fnum,'(" RDMFT exchange-correlation type : ",I4)') rdmxctype
    write(fnum,'(" Mueller-type functional, exponent : ",G18.10)') rdmalpha
  endif
end if
if (xctype.lt.0) then
  write(fnum,'("Optimised effective potential (OEP) and exact exchange (EXX)")')
  write(fnum,'(" Phys. Rev. B 53, 7024 (1996)")')
  write(fnum,'("Correlation type : ",I4)') abs(xctype)
  write(fnum,'(" ",A)') trim(xcdescr)
else
  write(fnum,'("Exchange-correlation type : ",I4)') xctype
  write(fnum,'(" ",A)') trim(xcdescr)
end if
if (xcgrad.eq.1) write(fnum,'(" Generalised gradient approximation (GGA)")')
write(fnum,*)
write(fnum,'("Smearing scheme :")')
write(fnum,'(" ",A)') trim(sdescr)
write(fnum,'("Smearing width : ",G18.10)') swidth
write(fnum,*)
write(fnum,'("Radial integration step length : ",I4)') lradstp
return
end subroutine
!EOC

