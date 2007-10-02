
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gndstate
! !INTERFACE:
subroutine gndstate
  ! !USES:
  use modmain
  use modmpi
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
  logical init,exist,redoscl
  integer ik,is,ia,idm,n,occsvfile
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
  ! initialise OEP/Hartree-Fock variables if required
  if ((xctype.lt.0).or.hartfock) call init2
  ! write the real and reciprocal lattice vectors to file
  call gndstate_write_latticevectors
  if ((task.eq.1).or.(task.eq.3)) then
     call readstate
     if(rank.eq.0) write(60,'("Density and potential read in from STATE.OUT")')
  else
     call rhoinit
     call poteff
     call genveffig
     if(rank.eq.0) write(60,'("Density and potential initialised from atomic data")')
  end if
  if(rank.eq.0)call flushifc(60)
  iscl=0
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
  redoscl=.false.
10 continue
  ! set last iteration flag
  tlast=.false.
  ! initialise the mixer
  init=.true.
  call packeff(.true.,n,nu)
  call mixer(init,beta0,betamax,n,nu,mu,beta,f,dv)
  call packeff(.false.,n,nu)
  ! delete any existing eigenvector files
 if (rank.eq.0) call delevec
  ! begin the self-consistent loop
  if(rank.eq.0) then
     write(60,*)
     write(60,'("+------------------------------+")')
     write(60,'("| Self-consistent loop started |")')
     write(60,'("+------------------------------+")')
     if (procs.gt.1) write(60,'("MPI parallelisation on",I5,"processes")')procs
  endif
  do iscl=1,maxscl
     if(rank.eq.0) then
        write(60,*)
        write(60,'("+--------------------------+")')
        write(60,'("| Iteration number : ",I5," |")') iscl
        write(60,'("+--------------------------+")')
		call flushifc(60)
     endif
     if (iscl.ge.maxscl) then
        if(rank.eq.0) then
           write(60,*)
           write(60,'("Reached self-consistent loops maximum")')
        endif
        tlast=.true.
     end if
     if(rank.eq.0) call flushifc(60)
     ! generate the core wavefunctions and densities
     call gndstate_gencore_wf_density
     ! begin parallel loop over k-points
     allocate(evalfv(nstfv,nspnfv))
     allocate(evecfv(nmatmax,nstfv,nspnfv))
     allocate(evecsv(nstsv,nstsv))
   
#ifdef MPI
	call MPI_barrier(MPI_COMM_WORLD,ierr)
     if (rank.eq.0) call delevec()
#endif
#ifdef MPISEC
     splittfile=.true.
       !$OMP PARALLEL DEFAULT(SHARED) &
     !$OMP PRIVATE(evalfv,evecfv,evecsv)
     !$OMP DO
     do ik=firstk(rank),lastk(rank)
     
#endif
#ifndef MPISEC
     splittfile=.false.
       !$OMP PARALLEL DEFAULT(SHARED) &
     !$OMP PRIVATE(evalfv,evecfv,evecsv)
     !$OMP DO
     do ik=1,nkpt
#endif
        ! solve the first- and second-variational secular equations
        call seceqn(ik,evalfv,evecfv,evecsv)
        ! write the eigenvalues/vectors to file
        call putevalfv(ik,evalfv)
        call putevalsv(ik,evalsv(1,ik))
        call putevecfv(ik,evecfv)
        call putevecsv(ik,evecsv)
     end do
     !$OMP END DO
     !$OMP END PARALLEL
     deallocate(evalfv,evecfv,evecsv)

     ! perform Hartree-Fock calculation if required
     if (hartfock) then
        ! initialise the occupancies
        if (iscl.le.1) call occupy
        !$OMP PARALLEL DEFAULT(SHARED) &
        !$OMP PRIVATE(evecsv)
        !$OMP DO
        do ik=firstk(rank),lastk(rank)
           allocate(evecsv(nstsv,nstsv))
           ! get the eigenvectors from file
           call getevecsv(vkl(1,ik),evecsv)
           ! solve the Hartree-Fock equations
           call seceqnhf(ik,evecsv)
           ! write the eigenvalues/vectors to file
           call putevalsv(ik,evalsv(1,ik))
           call putevecsv(ik,evecsv)
           deallocate(evecsv)
        end do
        !$OMP END DO
        !$OMP END PARALLEL
     end if
#ifdef MPISEC
     call mpisync_evalsv_spnchr
#endif
     ! find the occupation numbers and Fermi energy
     if (rank.eq.0) call occupy

#ifdef MPISEC
     call MPI_BCAST(occsv,nstsv*nkpt, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(fermidos,1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(efermi,1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

     if(rank.eq.0) then
        ! write out the eigenvalues and occupation numbers
        call writeeval
        ! write the Fermi energy to file
        call writefermi
        ! set the charge density and magnetisation to zero
     endif

     rhomt(:,:,:)=0.d0
     rhoir(:)=0.d0
     if (spinpol) then
        magmt(:,:,:,:)=0.d0
        magir(:,:)=0.d0
     end if
  
     
#ifdef MPIRHO	 

     do ik=firstk(rank),lastk(rank)
	call putoccsv(ik,occsv(1,ik))
	end do
 ! begin parallel loop over k-points
     !$OMP PARALLEL DEFAULT(SHARED) &
     !$OMP PRIVATE(evecfv,evecsv)
     !$OMP DO
     do ik=firstk(rank),lastk(rank)
#endif
#ifndef MPIRHO	
 if (rank.eq.0)then
    do ik=1,nkpt
	call putoccsv(ik,occsv(1,ik))
	end do
 	endif
 	! begin parallel loop over k-points
     !$OMP PARALLEL DEFAULT(SHARED) &
     !$OMP PRIVATE(evecfv,evecsv)
     !$OMP DO
	 do ik=1,nkpt	
#endif
        allocate(evecfv(nmatmax,nstfv,nspnfv))
        allocate(evecsv(nstsv,nstsv))
        ! write the occupancies to file
        
        ! get the eigenvectors from file
        call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
        call getevecsv(vkl(1,ik),evecsv)
        ! add to the density
        call rhovalk(ik,evecfv,evecsv)
        deallocate(evecfv,evecsv)
     end do
     !$OMP END DO
     !$OMP END PARALLEL
#ifdef MPIRHO    
	 call mpisumrhoandmag
#endif    
call mpiresumeevecfiles()
	 call gndstate_solvepotential(init,mu,nu,beta,f,n,dv)
     ! compute the energy components
     call energy
     call gndstate_output_state

     ! exit self-consistent loop if last iteration is complete
	 	  if (tlast) goto 20
     call gndstate_check_for_convergence(tlast,tstop,dv)
  !end of SCL
  end do
20 continue
  call gndstate_force_relax(redoscl)
  if(redoscl)goto 10
  call gndstate_output_timing_information
  if (rank.eq.0) then
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
  endif

  deallocate(nu,mu,beta,f)

  return
end subroutine gndstate
!EOC
