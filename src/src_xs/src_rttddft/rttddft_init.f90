! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! Created Jan 2021 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module implementing general initializations for RT-TDDFT
module rttddft_init
  use constants, only: zzero
  use m_gndstateq, only: gndstateq
  use MD, only: force, MD_input_keys
  use mod_APW_LO, only: apwordmax
  use mod_atoms, only: natmtot, nspecies, idxas, natoms
  use mod_eigensystem, only: nmatmax
  use mod_eigenvalue_occupancy, only: occsv, nstfv, nstsv
  use mod_gkvector, only: ngk, ngkmax, vgkl, gkc, tpgkc, sfacgk
  use mod_kpoint, only: vkl, nkpt
  use mod_muffin_tin, only: lmmaxapw
  use mod_misc, only: filext
  use modinput, only: input
  use modmpi
  use modxs, only: isreadstate0
  use physical_constants, only: c
  use precision, only: dp, i32
  use rttddft_CurrentDensity, only: UpdateCurrentDensity
  use rttddft_Density, only: updatedensity
  use rttddft_GlobalVariables
  use rttddft_HamiltonianOverlap, only: UpdateHam
  use rttddft_io, only: file_pmat_exists, read_pmat, write_pmat, &
    file_pmat_mt_exists, read_pmat_mt, write_pmat_mt, write_file_info, &
    write_file_info_fill_line_with_char
  use rttddft_pmat, only: Obtain_Pmat_LAPWLOBasis

  implicit none

  private

  public :: initialize_rttddft

contains
!> This subroutine initializes many global variables in a RT-TDDFT calculation.  
subroutine initialize_rttddft( molecular_dynamics )
  !> variable that is an interface to the input keys defined in `input.xml` inside the `MD` block
  type(MD_input_keys), intent(in)         :: molecular_dynamics

  integer                     :: ik, first_kpt, last_kpt
  character(len=50)           :: string
  character(len=*), parameter :: new_line = achar(13) // achar(10) 
  logical                     :: file_exists, readPmatBasis, forcePmatHermitian
  real(dp)                    :: voff(3)


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

  call distribute_loop(mpi_env_k, nkpt, first_kpt, last_kpt)

  ! Interface with input variables
  voff(1:3) = input%xs%vkloff(1:3)
  readPmatBasis = input%xs%realTimeTDDFT%readPmatbasis
  forcePmatHermitian = input%xs%realTimeTDDFT%forcePmatHermitian
  method = input%xs%realTimeTDDFT%propagator
  printTimesGeneral = input%xs%realTimeTDDFT%printTimingGeneral
  printTimesDetailed = (printTimesGeneral .and. input%xs%realTimeTDDFT%printTimingDetailed)
  calculateTotalEnergy = input%xs%realTimeTDDFT%calculateTotalEnergy
  calculateNexc = input%xs%realTimeTDDFT%calculateNExcitedElectrons
  predictorCorrector = associated(input%xs%realTimeTDDFT%predictorCorrector)
  if( predictorCorrector ) then
    tolPredCorr = input%xs%realTimeTDDFT%predictorCorrector%tol
    maxstepsPredictorCorrector = input%xs%realTimeTDDFT%predictorCorrector%maxIterations
  end if
  tstep = input%xs%realTimeTDDFT%timeStep
  tend = input%xs%realTimeTDDFT%endTime
  nsteps = int(tend/tstep)
  time = 0._dp

  ! Print to RTTDDFT_INFO that we will start the single-shot GS calculation
  if ( rank == 0 ) then
    call write_file_info_fill_line_with_char('=')
    call write_file_info( 'Non-self-consistent GS for TDDFT calculations - started' // new_line)
  end if

  ! Read from STATE.OUT exclusively
  isreadstate0 = .true.

  ! One-shot GS calculation
  call gndstateq (voff, '_RTTDDFT.OUT')

  call allocate_globals( first_kpt, last_kpt, ionDynamics=molecular_dynamics%on,&
    allocate_mathcalH=molecular_dynamics%valence_corrections, &
    allocate_mathcalB=molecular_dynamics%valence_corrections .or. molecular_dynamics%basis_derivative,&
    allocate_pmatmt=molecular_dynamics%valence_corrections .or. molecular_dynamics%basis_derivative, &
    allocate_B=molecular_dynamics%basis_derivative )
  
  if ( rank == 0 ) call write_to_info( molecular_dynamics%on )

  string = filext
  filext = '_RTTDDFT.OUT'

  call readstate        ! read the density and potentials from file
  call gencore          ! generate the core wavefunctions and densities
  call genmeffig
  call linengy          ! find the new linearization energies
  call genapwfr         ! generate the APW radial functions
  call genlofr          ! generate the local-orbital radial functions
  call olprad           ! compute the overlap radial integrals



  ! Get the eigenvectors and occupations from file
  do ik = first_kpt, last_kpt
    ! Eigenvectors (first and second-variational components)
    call getevecfv(vkl(:,ik), vgkl(:,:,:,ik), evecfv_gnd(:,:,ik))
    evecfv_time(1:nmatmax,1:nstfv,ik) = evecfv_gnd(1:nmatmax,1:nstfv,ik)
    call getevecsv(vkl(:,ik), evecsv(:,:,ik))
    call getoccsv(vkl(:,ik), occsv(:,ik))
  end do


  filext = string

  do ik = first_kpt, last_kpt
  ! Matching coefficients (apwalm)
    call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm(:,:,:,:,ik))
  end do


  ! Check if the momentum matrix elements have already been calculated
  file_exists = file_pmat_exists()
  if ( molecular_dynamics%on ) file_exists = file_exists .and. file_pmat_mt_exists()
  ! Calculate or read the momentum matrix elements
  if ( file_exists .and. readPmatBasis ) then
    call read_pmat( first_kpt, pmat, mpi_env_k )
    if ( molecular_dynamics%on ) call read_pmat_mt( first_kpt, pmatmt, mpi_env_k )
  else
    call Obtain_Pmat_LAPWLOBasis( forcePmatHermitian, molecular_dynamics%on )
    call write_pmat( first_kpt, pmat, mpi_env_k )
    if ( molecular_dynamics%on ) call write_pmat_mt( first_kpt, pmatmt, mpi_env_k )
  end if

  call init_laser

  ! Initialize fields
  pvec(:) = 0._dp
  jpara(:) = 0._dp
  jparaold(:) = 0._dp
  jdia(:) = 0._dp
  jind(:) = 0._dp
  aext(:) = 0._dp
  aind(:) = 0._dp
  atot(:) = 0._dp
  ! Hamiltonian at time t=0
  call UpdateHam( predcorr=.False., calculateOverlap=.True., &
    & update_mathcalH=allocated(mathcalH), update_mathcalB=allocated(mathcalB), update_pmat=.False. )
  ham_past(:,:,:) = ham_time(:,:,:)

  ! Spurious current
  if ( input%xs%realTimeTDDFT%subtractJ0 ) then
    call UpdateCurrentDensity(first_kpt,last_kpt,evecfv_gnd(:,:,:),jparaspurious(:))
  else
    jparaspurious(:) = 0._dp
  end if

end subroutine 

!> Allocate global arrays
subroutine allocate_globals( first_kpt, last_kpt, ionDynamics, allocate_mathcalH, &
    allocate_mathcalB, allocate_pmatmt, allocate_B )
  !> index of the first `k-point` to be considered in the sum
  integer(i32),intent(in)        :: first_kpt
  !> index of the last `k-point` considered
  integer(i32),intent(in)        :: last_kpt
  !> if `.True`, we need to allocate arrays for Ehrenfest molecular dynamics
  logical, intent(in) :: ionDynamics
  !> if `.True`, we need to allocate the global array `mathcalH`
  logical, intent(in) :: allocate_mathcalH
  !> if `.True`, we need to allocate the global array `mathcalB`
  logical, intent(in) :: allocate_mathcalB
  !> if `.True`, we need to allocate the global array `pmatmt`
  logical, intent(in) :: allocate_pmatmt
  !> if `.True`, we need to allocate the global arrays `B_time` and `B_past`
  logical, intent(in) :: allocate_B

  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,first_kpt:last_kpt))
  allocate(evecfv_gnd(nmatmax,nstfv,first_kpt:last_kpt), source=zzero)
  allocate(evecfv_time(nmatmax,nstfv,first_kpt:last_kpt))
  allocate(evecsv(nstsv,nstsv,first_kpt:last_kpt))
  allocate(overlap(nmatmax,nmatmax,first_kpt:last_kpt), source=zzero)
  allocate(ham_time(nmatmax,nmatmax,first_kpt:last_kpt), source=zzero)
  allocate(ham_past(nmatmax,nmatmax,first_kpt:last_kpt), source=zzero)
  allocate(pmat(nmatmax,nmatmax,3,first_kpt:last_kpt))
  if ( predictorCorrector ) then 
    allocate(ham_predcorr(nmatmax,nmatmax,nkpt), source=zzero)
    allocate(evecfv_save(nmatmax,nstfv,nkpt))
  end if
  if ( ionDynamics ) then
    if (allocate_mathcalH) allocate(mathcalH(nmatmax,nmatmax,3,natmtot,last_kpt))
    if (allocate_mathcalB) allocate(mathcalB(nmatmax,nmatmax,3,natmtot,first_kpt:last_kpt))
    if (allocate_pmatmt) allocate(pmatmt(nmatmax,nmatmax,3,natmtot,first_kpt:last_kpt))
    if (allocate_B) then
      allocate(B_time(nmatmax,nmatmax,first_kpt:last_kpt), source=zzero)
      allocate(B_past(nmatmax,nmatmax,first_kpt:last_kpt), source=zzero)
    end if
  end if

end subroutine

!> Output general information to `RTTDDFT_INFO.OUT`
subroutine write_to_info( ionDynamics )
  !> Are we performing an MD calculation?
  logical, intent(in)         :: ionDynamics

  character(len=100)          :: string
  character(len=*), parameter :: formatMemory = '(A40,F12.1)'
  integer(i32), parameter     :: MB = 1048576

  call write_file_info( 'Non-self-consistent GS for TDDFT calculations - finished' )
  call write_file_info_fill_line_with_char('=')
  call write_file_info( 'Allocated memory (MiB per MPI process)' )
  write( string, formatMemory ) 'Coefficients to match LAPW functions:', dble( sizeof( apwalm ) / MB )
  call write_file_info( string )
  write( string, formatMemory ) 'Wavefunctions:', &
    dble((sizeof(evecfv_gnd)+sizeof(evecfv_time)+sizeof(evecsv))/MB)
  call write_file_info( string )
  write( string, formatMemory ) 'Hamiltonian and Overlap matrices:', &
    dble((sizeof(overlap)+sizeof(ham_time)+sizeof(ham_past))/MB)
  call write_file_info( string )
  if( predictorCorrector ) then
    write( string, formatMemory ) 'Extra storage (predictor-corrector):', &
      dble((sizeof(ham_predcorr)+sizeof(evecfv_save))/MB)
    call write_file_info( string )
  end if
  write( string, formatMemory ) 'Momentum matrix:', dble((sizeof(pmat))/MB)
  call write_file_info( string )
  if( ionDynamics ) then
    call write_file_info( string )
    write( string, formatMemory ) 'Molecular Dynamics - Muffin-tin aux. matrices:', &
      dble( (sizeof(pmatmt)+sizeof(mathcalH)+sizeof(mathcalB)+sizeof(B_time)+sizeof(B_past))/MB )
    call write_file_info( string )
  end if
  call write_file_info_fill_line_with_char('=')
  ! General info to be printed to RTTDDFT_INFO
  call write_file_info( 'Important output files: AVEC.OUT, PVEC.OUT, JIND.OUT.' )
  call write_file_info( 'JIND.OUT contains the x, y, and z components of the current density.' )
  call write_file_info( 'PVEC.OUT contains the x, y, and z components of the polarization vector.' )
  call write_file_info( 'AVEC.OUT contains in each line 6 elements:' )
  call write_file_info( ': the x components of the induced and the total vector potential.' )
  call write_file_info( ': the y components of the induced and the total vector potential.' )
  call write_file_info( ': the z components of the induced and the total vector potential.' )
end subroutine

!> Initialize the most important variables related to the vector potential
!> applied by an external laser
subroutine init_laser

  integer(i32) :: ik

  if ( associated(input%xs%realTimeTDDFT%laser) ) then
    nkicks = size( input%xs%realTimeTDDFT%laser%kickarray )
    if ( nkicks >= 1 ) then
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
    if ( ntrapcos >= 1 ) then
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
    if ( nsinsq >= 1 ) then
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
end subroutine

end module rttddft_init
