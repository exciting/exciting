! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! History
! Created by Ronaldo Rodrigues Pela, July 2019
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module to manage "screenshots" of desired properties during a RT-TDDFT
!> propagation
module rttddft_screenshot
  implicit none
  private

  public :: screenshot

contains


  !> subroutine that takes "screenshots" of the Eigenvalues, and the Wavefunctions
  !> during a RT-TDDFT propagation
  !> @param[in]   it          number of the current iteration (to name output files)
  !> @param[in]   first_kpt   first k-point to be considered in the average
  !> @param[in]   last_kpt    last k-point
  !> @param[in]   overlap     overlap matrix
  !>                          Dimensions: nmatmax, nmatmax, first_kpt:last_kpt
  !> @param[in]   evec_gnd    coefficients of the KS-wavefunctions at t=0
  !>                          Dimensions: nmatmax, nstfv, first_kpt:last_kpt
  !> @param[in]   evec_time   coefficients of the KS-wavefunctions at current time
  !>                          Dimensions: nmatmax, nstfv, first_kpt:last_kpt
  !> @param[in]   ham_time    Hamiltonian matrix at time \( t \). Dimensions
  !>                          assumed for it: nmatmax, nmatmax, first_kpt:last_kpt
  subroutine screenshot( it, first_kpt, last_kpt, overlap, evecfv_gnd,&
      & evecfv_time, ham_time )
    use modmpi
    Use modinput, only: input
    use mod_misc, only: filext
    use mod_eigensystem, only: nmat, nmatmax
    use mod_eigenvalue_occupancy, only: nstfv
    use modinput, only: input
    use m_getunit, only: getunit
    use precision, only: dp
    use constants, only: zzero, zone

    implicit none

    !> number of the current iteration (to name output files)
    integer, intent(in)       :: it
    !> index of the first `k-point` to be considered in the sum
    integer,intent(in)        :: first_kpt
    !> index of the last `k-point` considered
    integer,intent(in)        :: last_kpt
    !> overlap matrix, Dimensions: `nmatmax`, `nmatmax`, `first_kpt:last_kpt`
    complex(dp), intent(in)   :: overlap(:, :, first_kpt:)
    !> Basis-expansion coefficients of the KS-wavefunctions at \( t=0 \).
    !> Dimensions: `nmatmax`, `nstfv`, `first_kpt:last_kpt`
    complex(dp), intent(in)   :: evecfv_gnd(:, :,first_kpt:)
    !> Basis-expansion coefficients of the KS-wavefunctions at current time.
    !> Dimensions: `nmatmax`, `nstfv`, `first_kpt:last_kpt`
    complex(dp), intent(in)   :: evecfv_time(:, :, first_kpt:)
    !> Hamiltonian matrix at time \( t \). 
    !> Dimensions: `nmatmax`, `nmatmax`, `first_kpt:last_kpt`
    complex(dp), intent(in)   :: ham_time(:, :, first_kpt:)


    logical                   :: print_abs ! Print just the abs**2 of the projection
    integer                   :: ik,ist,m,lwork,info,nmatp
    integer                   :: count
    integer,allocatable       :: ifail(:),iwork(:)
    integer                   :: fileout
    character(20)             :: strout
    character(20)             :: frmt
    character(20)             :: file_status
    real (dp)                 :: vl,vu
    real (dp),allocatable     :: w(:,:)
    complex(dp)               :: rwork(7*nmatmax)
    complex(dp)               :: scratch(nmatmax,nstfv)
    complex(dp)               :: proj_time(nstfv, nstfv, first_kpt:last_kpt )
    complex(dp), allocatable  :: work(:),evecham(:,:),hamcopy(:,:)
    complex(dp), allocatable  :: overlcopy(:,:)

    ! Print format
    print_abs = input%xs%realTimeTDDFT%screenshots%printAbsProjCoeffs
    if( print_abs ) then
      write(frmt,*)nstfv
    else
      write(frmt,*)2*nstfv
    end if
    frmt = adjustl(frmt)
    frmt = '('//trim(frmt)//trim('F10.5)')

    ! Adjustments about the file name
    write(strout,*) it
    strout = adjustl(strout)

#ifdef MPI
    if ( rank == 0 ) call getunit(fileout)
    call MPI_BCAST (fileout, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#else
    call getunit(fileout)
#endif



    ! Project the current WFs onto the ground-state ones
    do ik = first_kpt, last_kpt
      ! Matrix multiplication C := alpha*AB+beta*C
      call ZGEMM( 'N', 'N', nmatmax, nstfv, nmatmax, zone, overlap(:,:,ik), &
        & nmatmax, evecfv_time(:,:,ik), nmatmax, zzero, scratch, nmatmax )
      call ZGEMM( 'C', 'N', nstfv, nstfv, nmatmax, zone, evecfv_gnd(:,:,ik), &
        & nmatmax, scratch, nmatmax, zzero, proj_time(:,:,ik), nstfv )
    end do
    ! Write to file
    do count = 1, procs
      if ( rank == count - 1 ) then
        ! Trick: for merge to work, we need the same length for both strings
        file_status = merge('REPLACE', 'OLD    ', rank == 0)
        open( fileout, file = 'PROJ_'//trim( strout )//trim( filext ), &
          & action = 'WRITE', position = 'APPEND', status = trim( file_status ) )
        do ik = first_kpt, last_kpt
          write(fileout,'(A5,I10)') 'ik: ', ik
          do ist = 1, nstfv
            if( print_abs ) then
              write(fileout,frmt) abs(proj_time(:,ist,ik))**2
            else
              write(fileout,frmt) proj_time(:,ist,ik)
            end if
          end do
        end do
        close(fileout)
      end if
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
    end do ! do count = 1, procs

    ! Obtain the eigenvalues/eigenvectors of the hamiltonian at current time t
    allocate( w(nmatmax, first_kpt:last_kpt) )
    do ik = first_kpt, last_kpt
      nmatp = nmat(1,ik)
      allocate( overlcopy(nmatp, nmatp), hamcopy(nmatp, nmatp) )
      allocate( evecham(nmatp, nmatp) )
      allocate( ifail(nmatp), iwork(5*nmatp) )

      hamcopy(1:nmatp, 1:nmatp) = ham_time(1:nmatp, 1:nmatp, ik)
      overlcopy(1:nmatp, 1:nmatp) = overlap(1:nmatp, 1:nmatp, ik)

      ! Obtain the optimum lwork
      vl = 0._dp
      vu = 0._dp
      lwork = -1
      allocate( work(2) )
      call ZHEGVX( 1, 'V', 'I', 'U', nmatp, hamcopy, nmatp, overlcopy, nmatp, &
        & vl, vu, 1, nmatp, input%groundstate%solver%evaltol, m, w(:,ik), evecham, nmatp,&
        & work, lwork, rwork, iwork, ifail, info )
      lwork = int( work(1) )
      deallocate( work )
      allocate( work(lwork) )

      ! Solves A*x = (lambda)*B*x
      call ZHEGVX( 1, 'V', 'I', 'U', nmatp, hamcopy, nmatp, overlcopy, nmatp, &
        & vl, vu, 1, nmatp, input%groundstate%solver%evaltol, m, w(:,ik), evecham, nmatp,&
        & work, lwork, rwork, iwork, ifail, info )

      deallocate( work, overlcopy, hamcopy )
      deallocate( evecham, ifail, iwork )
    end do

  ! Write the eigenvalues onto EIGVAL
#ifdef MPI
    do count = 1, procs
      if ( rank == count - 1 ) then
#endif
        ! Trick: for merge to work, we need the same length for both strings
        file_status = merge('REPLACE', 'OLD    ', rank == 0)
        open( fileout, file = 'EIGVAL_'//trim( strout )//trim( filext ), &
          & action = 'WRITE', position = 'APPEND', status = trim(file_status) )
        do ik = first_kpt, last_kpt
          nmatp = nmat(1, ik)
          write( fileout, '(A5,I7)' ) 'ik = ', ik
          do ist = 1, nmatp
            write( fileout, '(I5,F20.12)' ) ist, w(ist, ik)
          end do
          write( fileout, * ) ''
        end do
        close(fileout)
#ifdef MPI
      end if ! if (rank .eq. count-1) then
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
    end do ! do count = 1, procs
#endif

  end subroutine screenshot

end module rttddft_screenshot
