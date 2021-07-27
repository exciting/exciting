! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! Created Jan 2021 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module implementing general initializations for RT-TDDFT
module rttddft_init
  implicit none

  private

  public :: initialize_rttddft

contains
!> This subroutine initializes many global variables in a RT-TDDFT calculation.  
!> TODO (Ronaldo, issue #82): split `initialize_rttddft` into subroutines.
subroutine initialize_rttddft
  use modmpi
  use modinput
  use rttddft_GlobalVariables
  use mod_kpoint, only: vkl
  use mod_gkvector, only: ngk, ngkmax, vgkl, gkc, tpgkc, sfacgk
  use mod_eigenvalue_occupancy, only: occsv, nstfv, nstsv
  use mod_muffin_tin, only: lmmaxapw
  use mod_eigensystem, only: nmatmax
  use mod_APW_LO, only: apwordmax
  use mod_kpoint, only: nkpt
  use mod_misc, only: filext
  use mod_atoms, only: natmtot
  use rttddft_pmat, only: Obtain_Pmat_LAPWLOBasis
  use m_getunit, only: getunit
  use rttddft_Density, only: updatedensity
  use rttddft_HamiltonianOverlap, only: UpdateHam
  use rttddft_CurrentDensity, only: UpdateCurrentDensity
  use modxs, only: isreadstate0
  use m_gndstateq, only: gndstateq
  use precision, only: dp
  use constants, only: zzero
  use physical_constants, only: c

  integer                 :: ik, recl, first_kpt, last_kpt
#ifdef MPI
  integer                 :: count
#endif
  integer, parameter      :: MB = 1048576
  character(50)           :: string
  character(50),parameter :: formatMemory = '(A40,F12.1)'
  logical                 :: filepmatexists
  real(dp)                :: voff(3)


  ! Backup groundstate variables
  call backup0
  call backup1

  !--------------------------------------------!
  !     map xs parameters associated to gs     !
  !--------------------------------------------!
  if( input%xs%rgkmax == 0.d0 ) input%xs%rgkmax = input%groundstate%rgkmax
  call mapxsparameters

  ! Initialize universal variables
  call init0
  call init1
  call init2

#ifdef MPI
  first_kpt = firstk(rank)
  last_kpt = lastk(rank)
#else
  first_kpt = 1
  last_kpt = nkpt
#endif

  ! Print to RTTDDFT_INFO that we will start the single-shot GS calculation
  if ( rank == 0 ) then
    call printline(fileinfortddft, "=")
    write( fileinfortddft, '("One-shot GS runs for TDDFT calculations")' )
    write( fileinfortddft, * )
  end if

  ! Read from STATE.OUT exclusively
  isreadstate0 = .true.

  ! One-shot GS calculation
  voff(1:3) = input%xs%vkloff(1:3)
  call gndstateq (voff, '_RTTDDFT.OUT')

  ! Print to RTTDDFT_INFO that we finished the single-shot GS calculation
  if ( rank == 0 ) then
    write(fileinfortddft, '("One-shot GS runs for TDDFT calculations - finished")')
    call printline(fileinfortddft, "=")
  end if

  ! Allocate the matching coefficients, and print the amount of memory required
  ! to store them
  if ( rank == 0 ) then
    write(fileinfortddft,*) 'Allocating memory (MiB per MPI process)'
  end if
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,first_kpt:last_kpt))
  if (rank == 0) then
    write(fileinfortddft,formatMemory) &
      & 'Coefficients to match LAPW functions:', dble( sizeof( apwalm ) / MB )
  end if

  ! Allocate the coefficients of the KS wavefunctions (in terms of our basis)
  ! and print the amount of memory required to store them
  allocate(evecfv_gnd(nmatmax,nstfv,first_kpt:last_kpt))
  allocate(evecfv_time(nmatmax,nstfv,first_kpt:last_kpt))
  allocate(evecsv(nstsv,nstsv,first_kpt:last_kpt))
  if ( rank == 0 ) then
    write(fileinfortddft,formatMemory) 'Wavefunctions:', &
    & dble((sizeof(evecfv_gnd)+sizeof(evecfv_time)+sizeof(evecsv))/MB)
  end if

  ! Allocate the overlap and hamiltonian matrices (for time t, for time t-\Delta t)
  ! and print the amount of memory required to store them
  allocate(overlap(nmatmax,nmatmax,first_kpt:last_kpt))
  allocate(ham_time(nmatmax,nmatmax,first_kpt:last_kpt))
  allocate(ham_past(nmatmax,nmatmax,first_kpt:last_kpt))
  if ( rank == 0 ) then
    write(fileinfortddft,formatMemory) 'Hamiltonian and Overlap matrices:', &
    & dble((sizeof(overlap)+sizeof(ham_time)+sizeof(ham_past))/MB)
  end if

  ! Initiliaze
  evecfv_gnd(:,:,:) = zzero
  overlap(:,:,:) = zzero
  ham_time(:,:,:) = zzero
  ham_past(:,:,:) = zzero

  ! If we are required to follow the predictorCorrector scheme
  if ( associated(input%xs%realTimeTDDFT%predictorCorrector) ) then
    allocate(ham_predcorr(nmatmax,nmatmax,nkpt))
    allocate(evecfv_save(nmatmax,nstfv,nkpt))
    ham_predcorr(:,:,:) = zzero
    predictorCorrector = .True.
    tolPredCorr = input%xs%realTimeTDDFT%predictorCorrector%tol
    maxstepsPredictorCorrector = input%xs%realTimeTDDFT%predictorCorrector%maxIterations
    if ( rank == 0 ) then
      write(fileinfortddft,formatMemory) 'Extra storage (predictor-corrector):', &
      & dble((sizeof(ham_predcorr)+sizeof(evecfv_save))/MB)
    end if
  else
    predictorCorrector = .False.
  end if

  ! Allocate the momentum matrix and print the amount of memory required to store it
  allocate(pmat(nmatmax,nmatmax,3,first_kpt:last_kpt))
  if ( rank == 0 ) then
    write(fileinfortddft,formatMemory) 'Momentum matrix:', &
    & dble((sizeof(pmat))/MB)
  end if

  ! Some general info to be printed to RTTDDFT_INFO
  if ( rank == 0 ) then
    call printline( fileinfortddft, "=" )
    write( fileinfortddft, * ) 'Important output files: AVEC.OUT, PVEC.OUT, JIND.OUT.'
    write( fileinfortddft, * ) 'JIND.OUT contains the x, y, and z components of the current density.'
    write( fileinfortddft, * ) 'PVEC.OUT contains the x, y, and z components of the polarization vector.'
    write( fileinfortddft, * ) 'AVEC.OUT contains in each line 6 elements:'
    write( fileinfortddft, * ) ': the x components of the induced and the total vector potential.'
    write( fileinfortddft, * ) ': the y components of the induced and the total vector potential.'
    write( fileinfortddft, * ) ': the z components of the induced and the total vector potential.'
  end if

  ! Initiliazing global variables
  method = input%xs%realTimeTDDFT%propagator
  printTimesGeneral = input%xs%realTimeTDDFT%printTimingGeneral
  printTimesDetailed = (printTimesGeneral .and. input%xs%realTimeTDDFT%printTimingDetailed)
  calculateTotalEnergy = input%xs%realTimeTDDFT%calculateTotalEnergy
  calculateNexc = input%xs%realTimeTDDFT%calculateNExcitedElectrons


  string = filext
  filext = '_RTTDDFT.OUT'

  call readstate        ! read the density and potentials from file
  call gencore          ! generate the core wavefunctions and densities
  call genmeffig
  call linengy          ! find the new linearization energies
  call genapwfr         ! generate the APW radial functions
  call genlofr          ! generate the local-orbital radial functions
  call olprad           ! compute the overlap radial integrals



  ! Get the eigenvectors and eigenvalues from file
  do ik = first_kpt, last_kpt
    ! Eigenvectors (first and second-variational components)
    call getevecfv(vkl(:,ik), vgkl(:,:,:,ik), evecfv_gnd(:,:,ik))
    evecfv_time(1:nmatmax,1:nstfv,ik) = evecfv_gnd(1:nmatmax,1:nstfv,ik)
    call getevecsv(vkl(:,ik), evecsv(:,:,ik))
    call getoccsv (vkl(:,ik), occsv(:,ik))
  end do


  filext = string

  do ik = first_kpt, last_kpt
  ! Matching coefficients (apwalm)
    call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm(:,:,:,:,ik))
  end do


  ! Check if the momentum matrix elements have already been calculated
  inquire(file='PMATBASIS.OUT', exist=filepmatexists)
  ! Calculate or read the momentum matrix elements
  if ( filepmatexists .and. input%xs%realTimeTDDFT%readPmatbasis ) then
#ifdef MPI
    call getunit(filepmat)
    do count = 1, procs
      if (rank .eq. count-1) then
        do ik = first_kpt, last_kpt
          inquire (IoLength=recl) pmat(:,:,:,ik)
          open(filepmat,file='PMATBASIS'//trim(filext),action='READ', &
            & form='UNFORMATTED',access='DIRECT',recl=recl)
          read(filepmat,rec=ik) pmat(:,:,:,ik)
        end do
        close(filepmat)
      end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end do
#else
    call getunit(filepmat)
    do ik = 1, nkpt
      inquire (IoLength=recl) pmat(:,:,:,ik)
      open(filepmat,file='PMATBASIS'//trim(filext),action='READ', &
        & form='UNFORMATTED',access='DIRECT',recl=recl)
      read(filepmat,rec=ik) pmat(:,:,:,ik)
    end do
    close(filepmat)
#endif
  else !if ( filepmatexists .and. input%xs%rt_tddft%readpmatbasis ) then
    call Obtain_Pmat_LAPWLOBasis
#ifdef MPI
    call getunit(filepmat)
    do count = 1, procs
      if (rank .eq. count-1) then
        do ik = first_kpt, last_kpt
          inquire (IoLength=recl) pmat(:,:,:,ik)
          open(filepmat,file='PMATBASIS'//trim(filext),action='WRITE',&
            & form='UNFORMATTED',access='DIRECT',recl=recl)
          write(filepmat,rec=ik) pmat(:,:,:,ik)
        end do
        close(filepmat)
      end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    end do
#else
    call getunit(filepmat)
    do ik = 1, nkpt
      inquire (IoLength=recl) pmat(:,:,:,ik)
      open(filepmat,file='PMATBASIS'//trim(filext),action='WRITE',&
        & form='UNFORMATTED',access='DIRECT',recl=recl)
      write(filepmat,rec=ik) pmat(:,:,:,ik)
    end do
    close(filepmat)
#endif
  end if !if ( filepmatexists .and. input%xs%rt_tddft%readpmatbasis ) then

  ! Setup time variables
  tstep = input%xs%realTimeTDDFT%timeStep
  tend = input%xs%realTimeTDDFT%endTime
  nsteps = int(tend/tstep)
  time = 0._dp

  ! Here, we initialize the most important variables about the vector potential
  ! applied by an external laser
  if ( associated(input%xs%realTimeTDDFT%laser) ) then
    nkicks = size( input%xs%realTimeTDDFT%laser%kickarray )
    if ( nkicks .ge. 1 ) then
      allocate(wkick(nkicks))
      allocate(dirkick(nkicks))
      allocate(amplkick(nkicks))
      allocate(t0kick(nkicks))
      do ik = 1, nkicks
        wkick(ik) = input%xs%realTimeTDDFT%laser%kickarray(ik)%kick%width
        dirkick(ik) = input%xs%realTimeTDDFT%laser%kickarray(ik)%kick%direction
        amplkick(ik) = -c*(input%xs%realTimeTDDFT%laser%kickarray(ik)%kick%amplitude)
        t0kick(ik) = input%xs%realTimeTDDFT%laser%kickarray(ik)%kick%t0
      end do
    end if
    ntrapcos = size( input%xs%realTimeTDDFT%laser%trapCosarray )
    if ( ntrapcos .ge. 1 ) then
      allocate(dirtrapcos(ntrapcos))
      allocate(ampltrapcos(ntrapcos))
      allocate(omegatrapcos(ntrapcos))
      allocate(phasetrapcos(ntrapcos))
      allocate(t0trapcos(ntrapcos))
      allocate(trtrapcos(ntrapcos))
      allocate(wtrapcos(ntrapcos))
      do ik = 1, ntrapcos
        dirtrapcos(ik) = input%xs%realTimeTDDFT%laser%trapCosarray(ik)%trapCos%direction
        ampltrapcos(ik) = input%xs%realTimeTDDFT%laser%trapCosarray(ik)%trapCos%amplitude
        omegatrapcos(ik) = input%xs%realTimeTDDFT%laser%trapCosarray(ik)%trapCos%omega
        phasetrapcos(ik) = input%xs%realTimeTDDFT%laser%trapCosarray(ik)%trapCos%phase
        t0trapcos(ik) = input%xs%realTimeTDDFT%laser%trapCosarray(ik)%trapCos%t0
        trtrapcos(ik) = input%xs%realTimeTDDFT%laser%trapCosarray(ik)%trapCos%riseTime
        wtrapcos(ik) = input%xs%realTimeTDDFT%laser%trapCosarray(ik)%trapCos%width
      end do
    end if
    nsinsq = size( input%xs%realTimeTDDFT%laser%sinSqarray )
    if ( nsinsq .ge. 1 ) then
      allocate(dirsinsq(nsinsq))
      allocate(amplsinsq(nsinsq))
      allocate(omegasinsq(nsinsq))
      allocate(phasesinsq(nsinsq))
      allocate(t0sinsq(nsinsq))
      allocate(tpulsesinsq(nsinsq))
      do ik = 1, nsinsq
        dirsinsq(ik) = input%xs%realTimeTDDFT%laser%sinSqarray(ik)%sinSq%direction
        amplsinsq(ik) = input%xs%realTimeTDDFT%laser%sinSqarray(ik)%sinSq%amplitude
        omegasinsq(ik) = input%xs%realTimeTDDFT%laser%sinSqarray(ik)%sinSq%omega
        phasesinsq(ik) = input%xs%realTimeTDDFT%laser%sinSqarray(ik)%sinSq%phase
        t0sinsq(ik) = input%xs%realTimeTDDFT%laser%sinSqarray(ik)%sinSq%t0
        tpulsesinsq(ik) = input%xs%realTimeTDDFT%laser%sinSqarray(ik)%sinSq%pulseLength
      end do
    end if
  end if

  ! Initialize fields
  pvec(:) = 0._dp
  jpara(:) = 0._dp
  jparaold(:) = 0._dp
  jdia(:) = 0._dp
  jind(:) = 0._dp
  aext(:) = 0._dp
  aind(:) = 0._dp
  atot(:) = 0._dp
  ! Open files: AVEC (where the vector potential is written), JIND (current density)
  ! PVEC (polarization vector), and, if it is the case, ETOT_RTTDDFT (total energy),
  ! NEXC (number of excited electrons)
  if ( rank == 0 ) then
    call getunit(fileavec)
    open(fileavec,file='AVEC'//trim(filext),status='replace')
    call getunit(filejind)
    open(filejind,file='JIND'//trim(filext),status='replace')
    call getunit(filepvec)
    open(filepvec,file='PVEC'//trim(filext),status='replace')
    if( calculateTotalEnergy ) then
      call getunit(fileetot)
      open(fileetot,file='ETOT_RTTDDFT'//trim(filext),status='replace')
    end if
    if( calculateNexc ) then
      call getunit(filenexc)
      open(filenexc,file='NEXC'//trim(filext),status='replace')
    end if
  end if


  ! Hamiltonian at time t=0
  call UpdateHam(.False.,.True.)
  ham_past(:,:,:) = ham_time(:,:,:)

  ! Spurious current
  if ( input%xs%realTimeTDDFT%subtractJ0 ) then
    call UpdateCurrentDensity(first_kpt,last_kpt,evecfv_gnd(:,:,:),jparaspurious(:))
  else
    jparaspurious(:) = 0.d0
  end if

end subroutine initialize_rttddft
end module rttddft_init
