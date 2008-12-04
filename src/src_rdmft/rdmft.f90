
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine rdmft
! 1-reduced density matrix functional theory
use modmain
implicit none
! local variables
integer ik
call init0
call init1
! generate q-point set and wiq2 array
call init2
! read density and potentials from file
call readstate
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! compute the overlap radial integrals
call olprad
! compute the Hamiltonian radial integrals
call hmlrad
! compute the kinetic energy of the core
call energykncr
! generate the kinetic matrix elements
call genkinmat
! read in the occupancies
do ik=1,nkpt
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
! calculate Coulomb potential matrix elements
call genvmat(vclmt,vclir,vclmat)
! derivative of kinetic energy w.r.t. evecsv
call rdmdkdc
! open information files
open(60,file='RDM_INFO.OUT',action='WRITE',form='FORMATTED')
! write out general information to RDM_INFO.OUT
call writeinfo(60)
! begin main self-consistent loop
do iscl=1,rdmmaxscl
  write(60,*)
  write(60,'("+-------------------------+")')
  write(60,'("| Iteration number : ",I4," |")') iscl
  write(60,'("+-------------------------+")')
  call flushifc(60)
! minimisation over natural orbitals
  if (maxitc.ge.1) then
    call rdmminc
    write(60,*)
    write(60,'("Natural orbital minimisation done")')
    call rdmwriteengy(60)
  end if
! minimisation over occupation number
  if (maxitn.ge.1) then
    call rdmminn
    write(60,*)
    write(60,'("Occupation number minimisation done")')
    call rdmwriteengy(60)
  end if
! end loop over iscl
end do
! write density to STATE.OUT
call writestate
! write occupation numbers for restart
do ik=1,nkpt
  call putoccsv(ik,occsv(:,ik))
end do
! close RDM_INFO.OUT file
close(60)
return
end subroutine

