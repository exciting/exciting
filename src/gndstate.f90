
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gndstate
! !INTERFACE:
subroutine gndstate
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the self-consistent Kohn-Sham ground-state. General information is
!   written to the file {\tt INFO.OUT}. First- and second-variational
!   eigenvalues, eigenvectors and occupancies are written to the unformatted
!   files {\tt EVALFV.OUT}, {\tt EVALSV.OUT}, {\tt EVECFV.OUT}, {\tt EVECSV.OUT}
!   and {\tt OCCSV.OUT}.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
logical exist
integer ik,is,ia,idm,n
real(8) dv,timetot
! allocatable arrays
real(8), allocatable :: nu(:)
real(8), allocatable :: mu(:)
real(8), allocatable :: beta(:)
real(8), allocatable :: f(:)
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
! require forces for structural optimisation
if ((task.eq.2).or.(task.eq.3)) tforce=.true.
! initialise universal variables
call init0
call init1
! initialise OEP variables if required
if (xctype.lt.0) call init2
! write the real and reciprocal lattice vectors to file
call writelat
! write inter-atomic distances to file
call writeiad
! write symmetry matrices to file
call writesym
! output the k-point set to file
call writekpts
! write lattice vectors and atomic positions to file
call writegeom(.false.)
! open INFO.OUT file
open(60,file='INFO'//trim(filext),action='WRITE',form='FORMATTED')
! open TOTENERGY.OUT
open(61,file='TOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
! open FERMIDOS.OUT
open(62,file='FERMIDOS'//trim(filext),action='WRITE',form='FORMATTED')
! open MOMENT.OUT if required
if (spinpol) open(63,file='MOMENT'//trim(filext),action='WRITE', &
 form='FORMATTED')
! open FORCEMAX.OUT if required
if (tforce) open(64,file='FORCEMAX'//trim(filext),action='WRITE', &
 form='FORMATTED')
! open RMSDVEFF.OUT
open(65,file='RMSDVEFF'//trim(filext),action='WRITE',form='FORMATTED')
! write out general information to INFO.OUT
call writeinfo(60)
! initialise or read the charge density and potentials from file
iscl=0
write(60,*)
if ((task.eq.1).or.(task.eq.3)) then
  call readstate
  write(60,'("Density and potential read in from STATE.OUT")')
else
  call rhoinit
  call poteff
  call genveffig
  write(60,'("Density and potential initialised from atomic data")')
end if
call flushifc(60)
! size of mixing vector
n=lmmaxvr*nrmtmax*natmtot+ngrtot
if (spinpol) n=n*(1+ndmag)
! allocate mixing arrays
allocate(nu(n))
allocate(mu(n))
allocate(beta(n))
allocate(f(n))
! set stop flag
tstop=.false.
10 continue
! set last iteration flag
tlast=.false.
! initialise the mixer
call packeff(.true.,n,nu)
call mixer(.true.,beta0,betamax,n,nu,mu,beta,f,dv)
call packeff(.false.,n,nu)
! delete any existing eigenvector files
call delevec
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
  call flushifc(60)
! generate the core wavefunctions and densities
  call gencore
! find the new linearisation energies
  call linengy
! write out the linearisation energies
  call writelinen
! generate the APW radial functions
  call genapwfr
! generate the local-orbital radial functions
  call genlofr
! compute the overlap radial integrals
  call olprad
! compute the Hamiltonian radial integrals
  call hmlrad
! begin parallel loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv)
!$OMP DO
  do ik=1,nkpt
! every thread should allocate its own arrays
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
    deallocate(evalfv,evecfv,evecsv)
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
! compute the effective potential
  call poteff
! pack interstitial and muffin-tin effective potential and field into one array
  call packeff(.true.,n,nu)
! mix in the old potential and field with the new
  call mixer(.false.,beta0,betamax,n,nu,mu,beta,f,dv)
! unpack potential and field
  call packeff(.false.,n,nu)
! add the fixed spin moment effective field
  if (fixspin) call fsmfield
! Fourier transform effective potential to G-space
  call genveffig
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
! output effective field for fixed spin moment calculations
  if (fixspin) then
    write(60,*)
    write(60,'("FSM effective field      : ",3G18.10)') bfsmc(1:ndmag)
  end if
! check for WRITE file
  inquire(file='WRITE',exist=exist)
  if (exist) then
    write(60,*)
    write(60,'("WRITE file exists - writing STATE.OUT")')
    call writestate
    open(50,file='WRITE')
    close(50,status='DELETE')
  end if
! write STATE.OUT file if required
  if (nwrite.ge.1) then
    if (mod(iscl,nwrite).eq.0) then
      call writestate
      write(60,*)
      write(60,'("Wrote STATE.OUT")')
    end if
  end if
! exit self-consistent loop if last iteration is complete
  if (tlast) goto 20
! check for convergence
  if (iscl.ge.2) then
    write(60,*)
    write(60,'("RMS change in effective potential (target) : ",G18.10,&
     &" (",G18.10,")")') dv,epspot
    if (dv.lt.epspot) then
      write(60,*)
      write(60,'("Potential convergence target achieved")')
      tlast=.true.
    end if
    write(65,'(G18.10)') dv
    call flushifc(65)
    if (xctype.lt.0) then
      write(60,'("Magnitude of OEP residue : ",G18.10)') resoep
    end if
  end if
! check for STOP file
  inquire(file='STOP',exist=exist)
  if (exist) then
    write(60,*)
    write(60,'("STOP file exists - stopping self-consistent loop")')
    tstop=.true.
    tlast=.true.
    open(50,file='STOP')
    close(50,status='DELETE')
  end if
! output the current total CPU time
  timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
  write(60,*)
  write(60,'("Time (CPU seconds) : ",F12.2)') timetot
! end the self-consistent loop
end do
20 continue
write(60,*)
write(60,'("+------------------------------+")')
write(60,'("| Self-consistent loop stopped |")')
write(60,'("+------------------------------+")')
! write density and potentials to file only if maxscl > 1
if (maxscl.gt.1) then
  call writestate
  write(60,*)
  write(60,'("Wrote STATE.OUT")')
end if
!-----------------------!
!     compute forces    !
!-----------------------!
if ((.not.tstop).and.(tforce)) then
  call force
! output forces to INFO.OUT
  call writeforce(60)
! write maximum force magnitude to FORCEMAX.OUT
  write(64,'(G18.10)') forcemax
  call flushifc(64)
end if
!---------------------------------------!
!     perform structural relaxation     !
!---------------------------------------!
if ((.not.tstop).and.((task.eq.2).or.(task.eq.3))) then
  write(60,*)
  write(60,'("Maximum force magnitude (target) : ",G18.10," (",G18.10,")")') &
   forcemax,epsforce
  call flushifc(60)
! check force convergence
  if (forcemax.le.epsforce) then
    write(60,*)
    write(60,'("Force convergence target achieved")')
    goto 30
  end if
! update the atomic positions if forces are not converged
  call updatpos
  write(60,*)
  write(60,'("+--------------------------+")')
  write(60,'("| Updated atomic positions |")')
  write(60,'("+--------------------------+")')
  do is=1,nspecies
    write(60,*)
    write(60,'("Species : ",I4," (",A,")")') is,trim(spsymb(is))
    write(60,'(" atomic positions (lattice) :")')
    do ia=1,natoms(is)
      write(60,'(I4,3F14.8)') ia,atposl(:,ia,is)
    end do
  end do
! add blank line to TOTENERGY.OUT, FERMIDOS.OUT and MOMENT.OUT
  write(61,*)
  write(62,*)
  if (spinpol) write (63,*)
! begin new self-consistent loop with updated positions
  goto 10
end if
30 continue
! output timing information
write(60,*)
write(60,'("Timings (CPU seconds) :")')
write(60,'(" initialisation                        : ",F12.2)') timeinit
write(60,'(" Hamiltonian and overlap matrix set up : ",F12.2)') timemat
write(60,'(" first-variational secular equation    : ",F12.2)') timefv
if (spinpol) then
  write(60,'(" second-variational calculation        : ",F12.2)') timesv
end if
write(60,'(" charge density calculation            : ",F12.2)') timerho
write(60,'(" potential calculation                 : ",F12.2)') timepot
if (tforce) then
  write(60,'(" force calculation                     : ",F12.2)') timefor
end if
timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
write(60,'(" total                                 : ",F12.2)') timetot
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
! close the FORCEMAX.OUT file
if (tforce) close(64)
! close the RMSDVEFF.OUT file
close(65)
deallocate(nu,mu,beta,f)
return
end subroutine
!EOC
