
subroutine putevalqp(fname, kset, ib, nb, eks, efks, eqp, efqp)

    use mod_kpointset
    use m_getunit
    implicit none
    character(*), intent(in) :: fname
    type(k_set),  intent(in) :: kset
    integer(4),   intent(in) :: ib
    integer(4),   intent(in) :: nb
    real(8),      intent(in) :: eks(ib:nb,kset%nkpt)
    real(8),      intent(in) :: efks
    real(8),      intent(in) :: eqp(ib:nb,kset%nkpt)
    real(8),      intent(in) :: efqp
    ! local
    integer(4) :: recl
    integer(4) :: ik
    integer(4) :: fid

    inquire(IoLength=recl) kset%nkpt, ib, nb, &
                           kset%vkl(:,1), &
                           eqp(ib:nb,1), &
                           eks(ib:nb,1), &
                           efqp, efks

    call getunit(fid)
    open(fid, File=trim(fname), Action='WRITE', Form='UNFORMATTED', &
         Access='DIRECT', status='REPLACE', Recl=recl)
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
