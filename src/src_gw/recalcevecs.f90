subroutine recalcevecs
    
    use modinput
    use modmain
    use modgw
    implicit none
    
! local variables
    Integer :: ik 
    real(8) :: tstart,tend

    Real (8), Allocatable :: evalfv (:, :)
    Complex (8), Allocatable :: evecfv (:, :, :)
    Complex (8), Allocatable :: evecsv (:, :)
    
    Integer(8) :: Recl
    Integer(4) :: nmatmax_, nstfv_, nspnfv_
    Real(8) :: vkl_(3)
    integer :: iostat, n
    Logical :: exist

!
!   Check first if the eigenvalues and eigenvectors are needed to be updated
!    
    Inquire (IoLength=Recl) vkl_, nmatmax_, nstfv_, nspnfv_
    Open (70, File='EVECFV.OUT', Action='READ', Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
    Read(70,Rec=1) vkl_, nmatmax_, nstfv_, nspnfv_
    Close (70)
    
    n=int(chgval/2.d0)+input%gw%nempty+1
    if ((nmatmax_.eq.nmatmax).and.(nstfv_.eq.n)) then
      allocate(evecfv(nmatmax_,nstfv_,nspnfv))
      Inquire (IoLength=Recl) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
      Open (70, File='EVECFV.OUT', Action='READ', Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
      Read(70,REC=nkptnr,IOSTAT=iostat) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
      close (70)
      deallocate(evecfv)
      if (iostat.eq.0) Return
    end if
!
!   Else one has to recalculate them accordingly
!    
    call cpu_time(tstart)
    
    input%groundstate%reducek=.false.
    input%groundstate%nempty=input%gw%nempty

    call init1
    call writekpts
    
    write(fgw,*)
    write(fgw,*) 'WARNING(genevecs):'
    write(fgw,*) '  KS eigenvalues and eigenvectors have been recalculated'
    write(fgw,*) '  for the complete k-grid and/or different number of the unoccupied states.'
    write(fgw,*)

! initialise the charge density and potentials from file
    Call readstate

! generate the core wavefunctions and densities
    Call gencore

! find the new linearisation energies
    Call linengy

! write out the linearisation energies
    Call writelinen

! generate the APW radial functions
    Call genapwfr

! generate the local-orbital radial functions
    Call genlofr

! compute the overlap radial integrals
    Call olprad

! compute the Hamiltonian radial integrals
    Call hmlint

! compute "relativistic mass" on the G-grid
      Call genmeffig

! delete any existing eigenvector files
    Call delevec
     
    Do ik = 1, nkpt

      Allocate(evalfv(nstfv, nspnfv))
      Allocate(evecfv(nmatmax, nstfv, nspnfv))
      Allocate(evecsv(nstsv, nstsv))

! solve the first- and second-variational secular equations
      Call seceqn(ik, evalfv, evecfv, evecsv)

! write the eigenvalues/vectors to file
      Call putevalfv(ik, evalfv)
      Call putevalsv(ik, evalsv(:, ik))
      Call putevecfv(ik, evecfv)
      Call putevecsv(ik, evecsv)

      Deallocate(evalfv, evecfv, evecsv)
    
    End Do ! ik

    if (allocated(meffig)) deallocate(meffig)
    if (allocated(m2effig)) deallocate(m2effig)
    input%groundstate%reducek=.true.
    call init1
    
    call cpu_time(tend)
    if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
    call write_cputime(fgw,tend-tstart,'GENEVECS')

10  Return

end subroutine
