!> Read quasiparticle energies from a direct-access binary file.
subroutine readevalqp(fname, kset, ib, nb, eks, efks, eqp, efqp)

    use mod_kpointset
    use modmpi, only: terminate_if_false
    use mod_large_io, only: inquire_large, open_direct_unformatted_large
    use precision, only: i32, dp, long_int
    implicit none
    !> File name to read.
    character(*), intent(in)  :: fname
    !> k-point set associated with the stored eigenvalues.
    type(k_set),  intent(in)  :: kset
    !> First band index to read.
    integer(i32), intent(in)  :: ib
    !> Last band index to read.
    integer(i32), intent(in)  :: nb
    !> Kohn-Sham eigenvalues.
    real(dp),     intent(out) :: eks(ib:nb,kset%nkpt)
    !> Kohn-Sham Fermi energy.
    real(dp),     intent(out) :: efks
    !> Quasiparticle energies.
    real(dp),     intent(out) :: eqp(ib:nb,kset%nkpt)
    !> Quasiparticle Fermi energy.
    real(dp),     intent(out) :: efqp
    ! local
    integer(i32)  :: ik, nk0, ib0, nb0, unit
    real(dp)      :: vkl(3)
    integer(long_int) :: recl
    logical       :: exist

    !-----------------------------------------------------------------------------
    ! Read the file
    !-----------------------------------------------------------------------------
    inquire( File=trim(fname), Exist=exist )
    call terminate_if_false(exist, 'ERROR(readevalqp): File ' // trim(fname) // ' does not exist!')

    call inquire_large( recl, [nk0, ib0, nb0] )
    call open_direct_unformatted_large( unit, trim(fname), "read", recl, "old" )
    read(unit, Rec=1) nk0, ib0, nb0
    close(unit)

    ! Consistency check
    call terminate_if_false(nk0 == kset%nkpt, 'ERROR(readevalqp): Inconsistent number of k-points!')
    call terminate_if_false(ib0 == ib, 'ERROR(readevalqp): Inconsistent first GW state in checkpoint!')
    call terminate_if_false(nb0 == nb, 'ERROR(readevalqp): Inconsistent last GW state in checkpoint!')

    call inquire_large( recl, [nk0, ib0, nb0], vkl, &
                        eqp(ib:nb,1), eks(ib:nb,1), [efqp, efks] )

    call open_direct_unformatted_large( unit, trim(fname), "read", recl, "old" )

    do ik = 1, kset%nkpt
        read(unit, Rec=ik) nk0, ib0, nb0, vkl, &
             eqp(ib:nb,ik), eks(ib:nb,ik), &
             efqp, efks
        call terminate_if_false(abs(sum(vkl(:)-kset%vkl(:,ik))) <= 1.d-6, &
            'ERROR(readevalqp): Inconsistent k-points!')
    end do ! ik

    close(unit)

end subroutine
