


! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writeinfo
! !INTERFACE:


subroutine writeinfo(fnum)
! !USES:
use modinput
use modmain
#ifdef TETRA
use modtetra
#endif
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
integer::fnum
! local variables
integer::i, is, ia
character(10)::dat, tim
write(fnum, '(" +---------------------------+")')
write(fnum, '(" | EXCITING hydrogen started |")')
#include "version.inc"
write(fnum, *)"| git hash id ", GITHASH, " |"
write(fnum, '(" +---------------------------+")')
if (notelns.gt.0) then
  write(fnum, *)
  write(fnum, '("Notes :")')
  do i=1, notelns
    write(fnum, '(A)') notes(i)
  end do
end if
call date_and_time(date=dat, time=tim)
write(fnum, *)
write(fnum, '("Date (YYYY-MM-DD) : ", A4, "-", A2, "-", A2)') dat(1:4), dat(5:6), &
 dat(7:8)
write(fnum, '("Time (hh:mm:ss)	 : ", A2, ":", A2, ":", A2)') tim(1:2), tim(3:4), &
 tim(5:6)
write(fnum, *)
write(fnum, '("All units are atomic (Hartree, Bohr, etc.)")')
select case(task)
case(0)
  write(fnum, *)
  write(fnum, '("+-------------------------------------------------+")')
  write(fnum, '("| Ground-state run starting from atomic densities |")')
  write(fnum, '("+-------------------------------------------------+")')
case(1, 200)
  write(fnum, *)
  write(fnum, '("+------------------------------------------+")')
  write(fnum, '("| Ground-state run resuming from STATE.OUT |")')
  write(fnum, '("+------------------------------------------+")')
case(2)
  write(fnum, *)
  write(fnum, '("+--------------------------------------------------------+")')
  write(fnum, '("| Structural optimisation starting from atomic densities |")')
  write(fnum, '("+--------------------------------------------------------+")')
case(3)
  write(fnum, *)
  write(fnum, '("+-----------------------------------------------------+")')
  write(fnum, '("| Structural optimisation run resuming from STATE.OUT |")')
  write(fnum, '("+-----------------------------------------------------+")')
case(5, 6)
  write(fnum, *)
  write(fnum, '("+-------------------------------+")')
  write(fnum, '("| Ground-state Hartree-Fock run |")')
  write(fnum, '("+-------------------------------+")')
case(300)
  write(fnum, *)
  write(fnum, '("+----------------------------------------------+")')
  write(fnum, '("| Reduced density matrix functional theory run |")')
  write(fnum, '("+----------------------------------------------+")')
case default
  write(*, *)
  write(*, '("Error(writeinfo): task not defined : ", I8)') task
  write(*, *)
  stop
end select
write(fnum, *)
write(fnum, '("Lattice vectors :")')
write(fnum, '(3G18.10)') input%structure%crystal%basevect(1, 1), input%structure%crystal%basevect(2, 1), &
    &input%structure%crystal%basevect(3, 1)
write(fnum, '(3G18.10)') input%structure%crystal%basevect(1, 2), input%structure%crystal%basevect(2, 2), &
    &input%structure%crystal%basevect(3, 2)
write(fnum, '(3G18.10)') input%structure%crystal%basevect(1, 3), input%structure%crystal%basevect(2, 3), &
    &input%structure%crystal%basevect(3, 3)
write(fnum, *)
write(fnum, '("Reciprocal lattice vectors :")')
write(fnum, '(3G18.10)') bvec(1, 1), bvec(2, 1), bvec(3, 1)
write(fnum, '(3G18.10)') bvec(1, 2), bvec(2, 2), bvec(3, 2)
write(fnum, '(3G18.10)') bvec(1, 3), bvec(2, 3), bvec(3, 3)
write(fnum, *)
write(fnum, '("Unit cell volume      : ", G18.10)') omega
write(fnum, '("Brillouin zone volume : ", G18.10)') (twopi**3)/omega
if (input%structure%autormt) then
  write(fnum, *)
  write(fnum, '("Automatic determination of muffin-tin radii")')
  write(fnum, '(" parameters : ", 2G18.10)') input%groundstate%rmtapm
end if
if (input%groundstate%frozencore) then
  write(fnum, *)
  write(fnum, '("A frozen-core calculation is performed")')
end if
do is=1, nspecies
  write(fnum, *)
  write(fnum, '("Species : ", I4, " (", A, ")")') is, trim(spsymb(is))
  write(fnum, '(" parameters loaded from : ", A)') trim(input%structure%speciesarray(is)%species%speciesfile)
  write(fnum, '(" name : ", A)') trim(spname(is))
  write(fnum, '(" nuclear charge    : ", G18.10)') spzn(is)
  write(fnum, '(" electronic charge : ", G18.10)') spze(is)
  write(fnum, '(" atomic mass : ", G18.10)') spmass(is)
  write(fnum, '(" muffin-tin radius : ", G18.10)') rmt(is)
  write(fnum, '(" number of radial points in muffin-tin : ", I6)') nrmt(is)
  write(fnum, '(" atomic positions (lattice), magnetic fields (Cartesian) :")')
  do ia=1, natoms(is)
    write(fnum, '(I4, " : ", 3F12.8, "	", 3F12.8)') ia, &
    &input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), &
     input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:)
  end do
end do
write(fnum, *)
write(fnum, '("Total number of atoms per unit cell : ", I4)') natmtot
write(fnum, *)
write(fnum, '("Spin treatment :")')
if (associated(input%groundstate%spin)) then
  write(fnum, '(" spin-polarised")')
else
  write(fnum, '(" spin-unpolarised")')
end if
if (isspinorb()) then
  write(fnum, '(" spin-orbit coupling")')
end if
if (associated(input%groundstate%spin)) then
  write(fnum, '(" global magnetic field (Cartesian) : ", 3G18.10)') input%groundstate%spin%bfieldc
  if (ncmag) then
    write(fnum, '(" non-collinear magnetisation")')
  else
    write(fnum, '(" collinear magnetisation in z-direction")')
  end if
end if
if (isspinspiral()) then
  write(fnum, '(" spin-spiral state assumed")')
  write(fnum, '("  q-vector (lattice)	: ", 3G18.10)') input%groundstate%spin%vqlss
  write(fnum, '("  q-vector (Cartesian) : ", 3G18.10)') vqcss
  write(fnum, '("  q-vector length	: ", G18.10)') sqrt(vqcss(1)**2 &
   +vqcss(2)**2+vqcss(3)**2)
end if
if (getfixspinnumber().ne.0) then
  write(fnum, '(" fixed spin moment (FSM) calculation")')
end if
if ((getfixspinnumber().eq.1).or.(getfixspinnumber().eq.3)) then
  write(fnum, '("  fixing total moment to (Cartesian) :")')
  write(fnum, '("  ", 3G18.10)') input%groundstate%spin%momfix
end if
if ((getfixspinnumber().eq.2).or.(getfixspinnumber().eq.3)) then
  write(fnum, '("  fixing local muffin-tin moments to (Cartesian) :")')
  do is=1, nspecies
    write(fnum, '("  species : ", I4, " (", A, ")")') is, trim(spsymb(is))
    do ia=1, natoms(is)
      write(fnum, '("	", I4, 3G18.10)') ia, input%structure%speciesarray(is)%species%atomarray(ia)%atom%mommtfix(:)
    end do
  end do
end if
write(fnum, *)
write(fnum, '("Number of Bravais lattice symmetries : ", I4)') nsymlat
write(fnum, '("Number of crystal symmetries	    : ", I4)') nsymcrys
write(fnum, *)
if (input%groundstate%autokpt) then
  write(fnum, '("radius of sphere used to determine k-point grid density : ", &
   &G18.10)') input%groundstate%radkpt
end if
write(fnum, '("k-point grid : ", 3I6)') input%groundstate%ngridk
write(fnum, '("k-point offset : ", 3G18.10)') input%groundstate%vkloff
if (input%groundstate%reducek) then
  write(fnum, '("k-point set is reduced with crystal symmetries")')
else
  write(fnum, '("k-point set is not reduced")')
end if
write(fnum, '("Total number of k-points : ", I8)') nkpt
write(fnum, *)
write(fnum, '("Smallest muffin-tin radius times maximum |G+k| : ", G18.10)') &
 input%groundstate%rgkmax
if ((input%groundstate%isgkmax.ge.1).and.(input%groundstate%isgkmax.le.nspecies)) then
  write(fnum, '("Species with smallest (or selected) muffin-tin radius : ", &
   &I4, " (", A, ")")') input%groundstate%isgkmax, &
    &trim(spsymb(input%groundstate%isgkmax))
end if
write(fnum, '("Maximum |G+k| for APW functions	     : ", G18.10)') gkmax
write(fnum, '("Maximum |G| for potential and density : ", G18.10)') input%groundstate%gmaxvr
write(fnum, '("Polynomial order for pseudocharge density : ", I4)') input%groundstate%npsden
write(fnum, *)
write(fnum, '("G-vector grid sizes : ", 3I6)') ngrid(1), ngrid(2), ngrid(3)
write(fnum, '("Total number of G-vectors : ", I8)') ngvec
write(fnum, *)
write(fnum, '("Maximum angular momentum used for")')
write(fnum, '(" APW functions			  : ", I4)') input%groundstate%lmaxapw
write(fnum, '(" computing H and O matrix elements : ", I4)') input%groundstate%lmaxmat
write(fnum, '(" potential and density		  : ", I4)') input%groundstate%lmaxvr
write(fnum, '(" inner part of muffin-tin	  : ", I4)') input%groundstate%lmaxinr
write(fnum, *)
write(fnum, '("Total nuclear charge    : ", G18.10)') chgzn
write(fnum, '("Total core charge       : ", G18.10)') chgcr
write(fnum, '("Total valence charge    : ", G18.10)') chgval
write(fnum, '("Total excess charge     : ", G18.10)') input%groundstate%chgexs
write(fnum, '("Total electronic charge : ", G18.10)') chgtot
write(fnum, *)
write(fnum, '("Effective Wigner radius, r_s : ", G18.10)') rwigner
write(fnum, *)
write(fnum, '("Number of empty states	      : ", I4)') input%groundstate%nempty
write(fnum, '("Total number of valence states : ", I4)') nstsv
write(fnum, *)
write(fnum, '("Total number of local-orbitals : ", I4)') nlotot
write(fnum, *)
if ((task.eq.5).or.(task.eq.6)) &
 write(fnum, '("Hartree-Fock calculation using Kohn-Sham states")')
if (input%groundstate%xctypenumber.lt.0) then
  write(fnum, '("Optimised effective potential (OEP) and exact exchange (EXX)")')
  write(fnum, '(" Phys. Rev. B 53, 7024 (1996)")')
  write(fnum, '("Correlation type : ", I4)') abs(input%groundstate%xctypenumber)
  write(fnum, '(" ", A)') trim(xcdescr)
else
  write(fnum, '("Exchange-correlation type : ", I4)') input%groundstate%xctypenumber
  write(fnum, '(" ", A)') trim(xcdescr)
end if
if (xcgrad.eq.1) write(fnum, '(" Generalised gradient approximation (GGA)")')
if (ldapu.ne.0) then
  write(fnum, *)
  write(fnum, '("LDA+U calculation")')
  if (ldapu.eq.1) then
    write(fnum, '(" fully localised limit (FLL)")')
  else if (ldapu.eq.2) then
    write(fnum, '(" around mean field (AFM)")')
  else if (ldapu.eq.3) then
    write(fnum, '(" interpolation between FLL and AFM")')
  else
    write(*, *)
    write(*, '("Error(writeinfo): ldapu not defined : ", I8)') ldapu
    write(*, *)
    stop
  end if
  write(fnum, '(" see PRB 67, 153106 (2003) and PRB 52, R5467 (1995)")')
  do is=1, nspecies
    if (llu(is).ge.0) then
      write(fnum, '(" species : ", I4, " (", A, ")", ", l = ", I2, ", U = ", F12.8, &
       &", J = ", F12.8)') is, trim(spsymb(is)), llu(is), ujlu(1, is), &
    &ujlu(2, is)
    end if
  end do
end if
if (task.eq.300) then
  write(fnum, *)
  write(fnum, '("RDMFT calculation")')
  write(fnum, '(" see arXiv:0801.3787v1 [cond-mat.mtrl-sci]")')
  write(fnum, '(" RDMFT exchange-correlation type : ", I4)') input%groundstate%RDMFT%rdmxctype
  if (input%groundstate%RDMFT%rdmxctype.eq.1) then
    write(fnum, '("  Hartree-Fock functional")')
  else if (input%groundstate%RDMFT%rdmxctype.eq.2) then
    write(fnum, '("  SDLG functional, exponent : ", G18.10)') input%groundstate%RDMFT%rdmalpha
  endif
end if
write(fnum, *)
write(fnum, '("Smearing scheme :")')
#ifdef XS
  if(associated(input%xs)) then
  if(associated(input%xs%tetra)) then
if (.not.input%xs%tetra%tetraocc) then
#endif
   write(fnum, '(" ", A)') trim(sdescr)
   write(fnum, '("Smearing width : ", G18.10)') input%groundstate%swidth
#ifdef XS
else
   write(fnum, '(" ", A)') 'No smearing - using the linear tetrahedron method'
   write(fnum, '(" ", A)') 'for occupation numbers and Fermi energy'
end if
endif
endif
#endif
write(fnum, *)
write(fnum, '("Radial integration step length : ", I4)') input%groundstate%lradstep
call flushifc(fnum)
return
end subroutine
!EOC
