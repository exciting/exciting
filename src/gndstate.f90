
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
  logical exist
  integer ik,is,ia,idm
  integer n,nwork
  real(8) timetot
  ! allocatable arrays
  real(8), allocatable :: v(:)
  real(8), allocatable :: work(:)
  real(8), allocatable :: evalfv(:,:)
  complex(8), allocatable :: evecfv(:,:,:)
  complex(8), allocatable :: evecsv(:,:)
  logical::force_converged, redoscl
  ! require forces for structural optimisation
  if ((task.eq.2).or.(task.eq.3)) tforce=.true.
  ! initialise global variables
  call init0
  call init1
  ! initialise OEP variables if required
  if (xctype.lt.0) call init2
  if(rank.eq.0) then
     ! write the real and reciprocal lattice vectors to file
     call writelat
     ! write interatomic distances to file
     call writeiad(.false.)
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
  end if
  iscl=0
  write(60,*)
  if ((task.eq.1).or.(task.eq.3)) then
     call readstate
     if(rank.eq.0) write(60,'("Potential read in from STATE.OUT")')
  else if (task.eq.200) then
     call phveff
     if(rank.eq.0) write(60,'("Supercell potential constructed from STATE.OUT")')
  else
     call rhoinit
     call poteff
     call genveffig
     if(rank.eq.0)  write(60,'("Density and potential initialised from atomic data")')
  end if
  call flushifc(60)
  ! size of mixing vector
  n=lmmaxvr*nrmtmax*natmtot+ngrtot
  if (spinpol) n=n*(1+ndmag)
  if (ldapu.ne.0) n=n+2*lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot
  ! allocate mixing arrays
  allocate(v(n))
  nwork=-1
allocate(work(1))
  call mixerifc(mixtype,n,v,currentconvergence,nwork,work)
deallocate(work)
  allocate(work(nwork))
  ! set stop flag
  tstop=.false.
10 continue
  ! set last iteration flag
  tlast=.false.
  ! delete any existing eigenvector files
  if((splittfile.or.rank.eq.0).and.(task.eq.0).or.(task.eq.2)) call delevec
  ! begin the self-consistent loop
  if(rank.eq.0) then
     write(60,*)
     write(60,'("+------------------------------+")')
     write(60,'("| Self-consistent loop started |")')
     write(60,'("+------------------------------+")')
  endif
  do iscl=1,maxscl
     if(rank.eq.0) then
        write(60,*)
        write(60,'("+-------------------------+")')
        write(60,'("| Iteration number : ",I4," |")') iscl
        write(60,'("+-------------------------+")')
     endif
     if (iscl.ge.maxscl) then
        if(rank.eq.0) then
           write(60,*)
           write(60,'("Reached self-consistent loops maximum")')
        endif
        tlast=.true.
#ifdef MPI
        call MPI_bcast(tlast,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif

     end if
     if(rank.eq.0) call flushifc(60)
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

#ifdef MPI
     call MPI_barrier(MPI_COMM_WORLD,ierr)
     if (rank.eq.0) call delevec()
#endif
#ifdef MPISEC
     splittfile=.true.
     do ik=firstk(rank),lastk(rank)

#endif
#ifdef NEVERDEFINED
     end do
#endif
#ifndef MPISEC
     splittfile=.false.
     ! begin parallel loop over k-points
#ifdef KSMP
     !$OMP PARALLEL DEFAULT(SHARED) &
     !$OMP PRIVATE(evalfv,evecfv,evecsv)
     !$OMP DO
#endif
     do ik=1,nkpt
#endif
        ! every thread should allocate its own arrays
        allocate(evalfv(nstfv,nspnfv))
        allocate(evecfv(nmatmax,nstfv,nspnfv))
        allocate(evecsv(nstsv,nstsv))
        ! solve the first- and second-variational secular equations
        call seceqn(ik,evalfv,evecfv,evecsv)
        ! write the eigenvalues/vectors to file
        call putevalfv(ik,evalfv)
        call putevalsv(ik,evalsv(:,ik))
        call putevecfv(ik,evecfv)
        call putevecsv(ik,evecsv)
        deallocate(evalfv,evecfv,evecsv)
     end do
#ifdef KSMP
     !$OMP END DO
     !$OMP END PARALLEL
#endif
#ifdef MPISEC
     call mpisync_evalsv_spnchr
#endif
     ! find the occupation numbers and Fermi energy
     call occupy
#ifdef MPISEC
     call MPI_BCAST(occsv,nstsv*nkpt, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(fermidos,1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(efermi,1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
     if (rank.eq.0) then
        ! write out the eigenvalues and occupation numbers
        call writeeval
        ! write the Fermi energy to file
        call writefermi
     endif
     ! set the charge density and magnetisation to zero
     rhomt(:,:,:)=0.d0
     rhoir(:)=0.d0
     if (spinpol) then
        magmt(:,:,:,:)=0.d0
        magir(:,:)=0.d0
     end if

#ifdef MPIRHO

     do ik=firstk(rank),lastk(rank)
        !write the occupancies to file
        call putoccsv(ik,occsv(:,ik))
     end do
     do ik=firstk(rank),lastk(rank)
#endif

#ifndef MPIRHO
        if (rank.eq.0)then
           do ik=1,nkpt
              !write the occupancies to file
              call putoccsv(ik,occsv(:,ik))
           end do
        endif
#ifdef KSMP
        ! begin parallel loop over k-points
        !$OMP PARALLEL DEFAULT(SHARED) &
        !$OMP PRIVATE(evecfv,evecsv)
        !$OMP DO
#endif
        do ik=1,nkpt
#endif
           allocate(evecfv(nmatmax,nstfv,nspnfv))
           allocate(evecsv(nstsv,nstsv))
           !
           ! get the eigenvectors from file
           call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
           call getevecsv(vkl(:,ik),evecsv)
           ! add to the density and magnetisation
           call rhovalk(ik,evecfv,evecsv)
           deallocate(evecfv,evecsv)
        end do
#ifndef MPIRHO
#ifdef KSMP
        !$OMP END DO
        !$OMP END PARALLEL
#endif
#endif
#ifdef MPIRHO
        call mpisumrhoandmag
#endif
#ifdef MPI
        if(xctype.lt.0) call mpiresumeevecfiles()
#endif
        ! symmetrise the density
        call symrf(lradstp,rhomt,rhoir)
        ! symmetrise the magnetisation
        if (spinpol) call symrvf(lradstp,magmt,magir)
        ! convert the density from a coarse to a fine radial mesh
        call rfmtctof(rhomt)
        ! convert the magnetisation from a coarse to a fine radial mesh
        do idm=1,ndmag
           call rfmtctof(magmt(:,:,:,idm))
        end do
        ! add the core density to the total density
        call addrhocr
        ! calculate the charges
        call charge
        ! calculate the moments
        if (spinpol) call moment
        ! normalise the density
        call rhonorm
        ! LDA+U
        if (ldapu.ne.0) then
           ! generate the LDA+U density matrix
           call gendmatlu
           ! generate the LDA+U potential matrix
           call genvmatlu
           ! write the LDA+U matrices to file
           call writeldapu
        end if
        ! compute the effective potential
        call poteff
        ! pack interstitial and muffin-tin effective potential and field into one array
        call packeff(.true.,n,v)
        ! mix in the old potential and field with the new

        if(rank.eq.0) call mixerifc(mixtype,n,v,currentconvergence,nwork,work)

#ifdef MPI
        call  MPI_BCAST(v(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call  MPI_BCAST(nwork, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call  MPI_BCAST(work(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)


!	call  MPI_BCAST(nu(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!        call  MPI_BCAST(mu(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!        call  MPI_BCAST(f(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
!        call  MPI_BCAST(beta(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
        ! unpack potential and field
        call packeff(.false.,n,v)
        ! add the fixed spin moment effect field
        if (fixspin.ne.0) call fsmfield
        ! Fourier transform effective potential to G-space
        call genveffig
        ! reduce the external magnetic fields if required
        if (reducebf.lt.1.d0) then
           bfieldc(:)=bfieldc(:)*reducebf
           bfcmt(:,:,:)=bfcmt(:,:,:)*reducebf
        end if
        ! compute the energy components
        call energy
        ! output energy components
        if(rank.eq.0) then
           call writeengy(60)
           write(60,*)
           write(60,'("Density of states at Fermi energy : ",G18.10)') fermidos
           write(60,'(" (states/Hartree/unit cell)")')
           ! write total energy to TOTENERGY.OUT and flush
           write(61,'(G22.12)') engytot
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
           ! output effective fields for fixed spin moment calculations
           if (fixspin.ne.0) call writefsm(60)
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
        endif
        if (tlast) goto 20
        if(rank.eq.0)then
           ! check for convergence
           if (iscl.ge.2) then
              write(60,*)
              write(60,'("RMS change in effective potential (target) : ",G18.10,&
                   &" (",G18.10,")")') currentconvergence,epspot
              if (currentconvergence.lt.epspot) then
                 write(60,*)
                 write(60,'("Potential convergence target achieved")')
                 tlast=.true.
              end if
              write(65,'(G18.10)') currentconvergence
              call flushifc(65)
           end if
           if (xctype.lt.0) then
              write(60,*)
              write(60,'("Magnitude of OEP residual : ",G18.10)') resoep
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
        endif
#ifdef MPI
        call MPI_bcast(tstop,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_bcast(tlast,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif

     end do


20   continue
     if (rank.eq.0)then
        redoscl=.false.
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
           if(rank.eq.0) then
              ! output forces to INFO.OUT
              call writeforce(60)
              ! write maximum force magnitude to FORCEMAX.OUT
              write(64,'(G18.10)') forcemax
              call flushifc(64)
           endif
        end if
     endif
     !---------------------------------------!
     !     perform structural relaxation     !
     !---------------------------------------!
     if ((.not.tstop).and.((task.eq.2).or.(task.eq.3))) then
        if(rank.eq.0) then
           write(60,*)
           write(60,'("Maximum force magnitude (target) : ",G18.10," (",G18.10,")")') &
                forcemax,epsforce
           call flushifc(60)
        endif
        ! check force convergence
        if(rank.eq.0)	then
           force_converged=.false.
           if (forcemax.le.epsforce) then
              write(60,*)
              write(60,'("Force convergence target achieved")')
              force_converged=.true.
           end if
        endif
#ifdef MPI
        call mpi_bcast(force_converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif
        if(force_converged) goto 30
        ! update the atomic positions if forces are not converged
        call updatpos
        if(rank.eq.0)	then
           write(60,*)
           write(60,'("+--------------------------+")')
           write(60,'("| Updated atomic positions |")')
           write(60,'("+--------------------------+")')
           do is=1,nspecies
              write(60,*)
              write(60,'("Species : ",I4," (",A,")")') is,trim(spsymb(is))
              write(60,'(" atomic positions (lattice) :")')
              do ia=1,natoms(is)
                 write(60,'(I4," : ",3F14.8)') ia,atposl(:,ia,is)
              end do
           end do
! add blank line to TOTENERGY.OUT, FERMIDOS.OUT, MOMENT.OUT and RMSDVEFF.OUT
           write(61,*)
           write(62,*)
           if (spinpol) write (63,*)
           write(65,*)
           ! begin new self-consistent loop with updated positions
           redoscl=.true.
        endif

#ifdef MPI
        call MPI_bcast(redoscl,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

        call MPI_Bcast(forcemax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(atposc,size(atposc),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(atposl,size(atposl),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(vkl,size(vkl),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(vkc,size(vkc),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ngvec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(sfacg,size(sfacg),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ngk,size(ngk),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(vgkc,size(vgkc),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(sfacgk,size(sfacgk),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ngkmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(engynn,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

#endif
        if(redoscl)  goto 10
     end if
30   continue
     ! output timing information
     if(rank.eq.0)then
        write(60,*)
        write(60,'("Timings (CPU seconds) :")')
        write(60,'(" initialisation",T40,": ",F12.2)') timeinit
        write(60,'(" Hamiltonian and overlap matrix set up",T40,": ",F12.2)') timemat
        write(60,'(" first-variational secular equation",T40,": ",F12.2)') timefv
        if (spinpol) then
           write(60,'(" second-variational calculation",T40,": ",F12.2)') timesv
        end if
        write(60,'(" charge density calculation",T40,": ",F12.2)') timerho
        write(60,'(" potential calculation",T40,": ",F12.2)') timepot
        if (tforce) then
           write(60,'(" force calculation",T40,": ",F12.2)') timefor
        end if
        timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
        write(60,'(" total",T40,": ",F12.2)') timetot
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
     endif
     deallocate(v,work)
     call mpiresumeevecfiles()
     return
   end subroutine gndstate
!EOC
