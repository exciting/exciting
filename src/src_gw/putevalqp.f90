
subroutine putevalqp(fname, kset, ib, nb, eks, efks, eqp, efqp)

    use mod_kpointset
    use mod_large_io, only: inquire_large, open_direct_unformatted_large
    use precision, only: i32, dp, long_int
    implicit none
    character(*), intent(in) :: fname
    type(k_set),  intent(in) :: kset
    integer(i32), intent(in) :: ib
    integer(i32), intent(in) :: nb
    real(dp),     intent(in) :: eks(ib:nb,kset%nkpt)
    real(dp),     intent(in) :: efks
    real(dp),     intent(in) :: eqp(ib:nb,kset%nkpt)
    real(dp),     intent(in) :: efqp
    ! local
    integer(long_int) :: recl
    integer(i32) :: ik
    integer(i32) :: fid

    call inquire_large( recl, [kset%nkpt, ib, nb], kset%vkl(:,1), &
                        eqp(ib:nb,1), eks(ib:nb,1), [efqp, efks] )

    call open_direct_unformatted_large( fid, trim(fname), "write", recl, "replace" )
    do ik = 1, kset%nkpt
        write(fid, Rec=ik) kset%nkpt, ib, nb, &
                           kset%vkl(:,ik), &
                           eqp(ib:nb,ik), &
                           eks(ib:nb,ik), &
                           efqp, efks
    end do ! ik
    close(fid)

    return
end subroutine
