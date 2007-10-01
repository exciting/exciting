
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hartfock
use modmain
implicit none
! local variables
logical exist
integer ik,idm
real(8) etp,de
! allocatable arrays
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
! initialise universal variables
call init0
call init1
call init2
! open INFO.OUT file
open(60,file='INFO'//trim(filext),action='WRITE',form='FORMATTED')
! open TOTENERGY.OUT
open(61,file='TOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
! open FERMIDOS.OUT
open(62,file='FERMIDOS'//trim(filext),action='WRITE',form='FORMATTED')
! open MOMENT.OUT if required
if (spinpol) open(63,file='MOMENT'//trim(filext),action='WRITE', &
 form='FORMATTED')
! open DENERGY.OUT
open(65,file='DENERGY'//trim(filext),action='WRITE',form='FORMATTED')
! write out general information to INFO.OUT
call writeinfo(60)
! read the charge density and potentials from file
call readstate
! compute the effective potential
call poteff
! Fourier transform effective potential to G-space
call genveffig
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
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv)
!$OMP DO
do ik=1,nkpt
  allocate(evalfv(nstfv,nspnfv))
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
! solve the first- and second-variational secular equations
  call seceqn(ik,evalfv,evecfv,evecsv)
! write the eigenvalues/vectors to file
  call putevalfv(ik,evalfv)
  call putevalsv(ik,evalsv(1,ik))
  call putevecfv(ik,evecfv)
  call putevecsv(ik,evecsv)
! generate the kinetic matrix elements
  call genkinmat(ik,evecfv,evecsv)
  deallocate(evalfv,evecfv,evecsv)
end do
!$OMP END DO
!$OMP END PARALLEL
! find the occupation numbers and Fermi energy
call occupy
! set last iteration flag
tlast=.false.
etp=0.d0
! begin the self-consistent loop
write(60,*)
write(60,'("+------------------------------+")')
write(60,'("| Self-consistent loop started |")')
write(60,'("+------------------------------+")')
do iscl=1,maxscl
  write(60,*)
  write(60,'("+--------------------------+")')
  write(60,'("| Iteration number : ",I5," |")') iscl
  write(60,'("+--------------------------+")')
  if (iscl.ge.maxscl) then
    write(60,*)
    write(60,'("Reached self-consistent loops maximum")')
    tlast=.true.
  end if
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecsv)
!$OMP DO
  do ik=1,nkpt
    allocate(evecsv(nstsv,nstsv))
    call getevecsv(vkl(1,ik),evecsv)
! solve the Hartree-Fock secular equation
    call seceqnhf(ik,evecsv)
! write the eigenvalues/vectors to file
    call putevalsv(ik,evalsv(1,ik))
    call putevecsv(ik,evecsv)
    deallocate(evecsv)
  end do
!$OMP END DO
!$OMP END PARALLEL
! find the occupation numbers and Fermi energy
  call occupy
! write out the eigenvalues and occupation numbers
  call writeeval
! write the Fermi energy to file
  call writefermi
! set the charge density and magnetisation to zero
  rhomt(:,:,:)=0.d0
  rhoir(:)=0.d0
  if (spinpol) then
    magmt(:,:,:,:)=0.d0
    magir(:,:)=0.d0
  end if
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecfv,evecsv)
!$OMP DO
  do ik=1,nkpt
    allocate(evecfv(nmatmax,nstfv,nspnfv))
    allocate(evecsv(nstsv,nstsv))
! write the occupancies to file
    call putoccsv(ik,occsv(1,ik))
! get the eigenvectors from file
    call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
    call getevecsv(vkl(1,ik),evecsv)
! add to the density and magnetisation
    call rhovalk(ik,evecfv,evecsv)
    deallocate(evecfv,evecsv)
  end do
!$OMP END DO
!$OMP END PARALLEL
! symmetrise the density
  call symrf(lradstp,rhomt,rhoir)
! symmetrise the magnetisation
  if (spinpol) call symrvf(lradstp,magmt,magir)
! convert the density from a coarse to a fine radial mesh
  call rfmtctof(rhomt)
! convert the magnetisation from a coarse to a fine radial mesh
  do idm=1,ndmag
    call rfmtctof(magmt(1,1,1,idm))
  end do
! add the core density to the total density
  call addrhocr
! calculate the charges
  call charge
! calculate the moments
  if (spinpol) call moment
! normalise the density
  call rhonorm
! compute the Coulomb potential
  call potcoul
! compute the energy components
  call energy
! output energy components
  call writeengy(60)
  write(60,*)
  write(60,'("Density of states at Fermi energy : ",G18.10)') fermidos
  write(60,'(" (states/Hartree/spin/unit cell)")')
! write total energy to TOTENERGY.OUT and flush
  write(61,'(G18.10)') engytot
  call flushifc(61)
! write DOS at Fermi energy to FERMIDOS.OUT and flush
  write(62,'(G18.10)') fermidos
  call flushifc(62)
! output charges and moments
  call writechg(60)
! write total moment to MOMENT.OUT and flush
  if (spinpol) then
    write(63,'(3G18.10)') momtot(1:ndmag)
    call flushifc(63)
  end if
  if (tlast) goto 10
! compute the change in total energy and check for convergence
  if (iscl.ge.2) then
    de=abs(engytot-etp)/(abs(engytot)+1.d0)
    write(60,*)
    write(60,'("Relative change in total energy (target) : ",G18.10,&
     &" (",G18.10,")")') de,epsengy
    if (de.lt.epsengy) then
      write(60,*)
      write(60,'("Energy convergence target achieved")')
      tlast=.true.
    end if
    write(65,'(G18.10)') de
    call flushifc(65)
  end if
  etp=engytot
! check for STOP file
  inquire(file='STOP',exist=exist)
  if (exist) then
    write(60,*)
    write(60,'("STOP file exists - stopping self-consistent loop")')
    tlast=.true.
    open(50,file='STOP')
    close(50,status='DELETE')
  end if
end do
10 continue
write(60,*)
write(60,'("+------------------------------+")')
write(60,'("| Self-consistent loop stopped |")')
write(60,'("+------------------------------+")')
if (maxscl.gt.1) then
  call writestate
  write(60,*)
  write(60,'("Wrote STATE.OUT")')
end if
write(60,*)
write(60,'("+----------------------------------+")')
write(60,'("| EXCITING version ",I1.1,".",I1.1,".",I3.3," stopped |")') version
write(60,'("+----------------------------------+")')
! close the INFO.OUT file
close(60)
! close the TOTENERGY.OUT file
close(61)
! close the FERMIDOS.OUT file
close(62)
! close the MOMENT.OUT file
if (spinpol) close(63)
! close the DENERGY.OUT file
close(65)
return
end subroutine

