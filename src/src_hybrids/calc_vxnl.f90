!
!BOP
! !ROUTINE: calc_vxnl
! !INTERFACE:
!
subroutine calc_vxnl()
! !USES:
    use modmain
    use mod_hybrids, only: vxnl, exnl, evalfv, kset, nomax, numin, ikvbm, ikcbm, ikvcm

    use modfvsystem
    use modmpi
!
! !DESCRIPTION:
!   Calculates the non-local exchange potential
!   and the non-local exchange energy for Hartree-Fock based hybrid functionals.
!
!EOP
!BOC
    implicit none

    integer(4) :: ik, ib
    integer(4) :: ikfirst, iklast, COMM_LEVEL_2
    real(8)    :: tstart, tend

    call cpu_time(tstart)

    !----------------------------------------
    ! Read KS eigenvalues from file EVALSV.OUT
    !----------------------------------------
    if (allocated(evalfv)) deallocate(evalfv)
    allocate(evalfv(nstfv,kset%nkpt))
    evalfv(:,:) = 0.d0
    do ik = 1, nkpt
      call getevalfv(kset%vkl(:,ik), evalfv(:,ik))
    end do

    ! VB / CB state index
    call find_vbm_cbm(1, nstfv, kset%nkpt, evalfv, efermi, nomax, numin, ikvbm, ikcbm, ikvcm)

    ! BZ integration weights
    call kintw()
    deallocate(evalfv)

#ifdef MPI
    if (rank == 0) write(*,*) "Computing non-local potential"
    ikfirst = firstofset(mod(rank, nkpt), nkpt)
    iklast = lastofset(mod(rank, nkpt), nkpt)
    write(*,*) "On proc", rank, "do ik", ikfirst," to ", iklast
    ! create communicator object for each k-point
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, ikfirst,  rank/nkpt, COMM_LEVEL_2, ierr)
    write(*,*) "proc ", rank, " does color ", ikfirst, " with key(rank)", rank/nkpt
#else
    ikfirst = 1
    iklast = nkpt
#endif

    ! Matrix elements of non-local potential   !
    if (allocated(vxnl)) deallocate(vxnl)
    allocate(vxnl(nstfv,nstfv,ikfirst:iklast))
    vxnl(:,:,:) = zzero

    !---------------------------------------
    ! Loop over k-points
    !---------------------------------------
    ! each process does a subset
    exnl = 0.d0
    do ik = ikfirst, iklast
      ! Calculate the non-local potential
      call calc_vxnl_k(ik, COMM_LEVEL_2)
      ! Calculate the non-local exchange energy
      do ib = 1, nomax
        exnl = exnl + kset%wkpt(ik)*vxnl(ib,ib,ik)
      end do
    end do ! ikp

#ifdef MPI
    call MPI_AllReduce(MPI_IN_PLACE, exnl, 1, &
    &                  MPI_DOUBLE_PRECISION, MPI_SUM, &
    &                  MPI_COMM_WORLD, ierr)
    call barrier()
#endif

    call cpu_time(tend)

    return
end subroutine
