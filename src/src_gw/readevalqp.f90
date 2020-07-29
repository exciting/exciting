
subroutine readevalqp(fname, kset, ib, nb, eks, efks, eqp, efqp)

    use mod_kpointset
    implicit none
    character(*), intent(in)  :: fname
    type(k_set),  intent(in)  :: kset
    integer(4),   intent(in)  :: ib
    integer(4),   intent(in)  :: nb
    real(8),      intent(out) :: eks(ib:nb,kset%nkpt)
    real(8),      intent(out) :: efks
    real(8),      intent(out) :: eqp(ib:nb,kset%nkpt)
    real(8),      intent(out) :: efqp
    ! local
    integer(4)    :: ik, nk0, ib0, nb0
    real(8)       :: vkl(3)
    integer(4)    :: recl
    logical       :: exist

    !-----------------------------------------------------------------------------
    ! Read the file
    !-----------------------------------------------------------------------------
    inquire( File=trim(fname), Exist=exist )
    if ( .not.exist ) then
        write(*,*)'ERROR(readevalqp): File ', trim(fname), ' does not exist!'
        stop
    end if

    inquire( IoLength=recl ) nk0, ib0, nb0
    open(70, File=trim(fname), Action='READ', Form='UNFORMATTED', &
         Access='DIRECT', Recl=recl)
    read(70, Rec=1) nk0, ib0, nb0
    close(70)

    ! Consistency check
    if ( nk0 /= kset%nkpt ) then
        write(*,*) 'ERROR(readevalqp): Inconsistent number of k-points!'
        write(*,*) 'nk0 = ', nk0, ' nkpt = ', kset%nkpt
        stop
    end if
    if ( ib0 /= ib ) then
        write(*,*) 'ERROR(readevalqp): Inconsistent number of GW states!'
        write(*,*) 'ib0 = ', ib0, ' ib = ', ib
        stop
    end if
    if ( nb0 /= nb0 ) then
        write(*,*) 'ERROR(readevalqp): Inconsistent number of GW states!'
        write(*,*) 'nb0 = ', nb0, ' nb = ', nb
        stop
    end if

    inquire( IoLength=recl) nk0, ib0, nb0, vkl, &
             eqp(ib:nb,1), eks(ib:nb,1), &
             efqp, efks

    open(70, File=trim(fname), Action='READ', Form='UNFORMATTED', &
         Access='DIRECT', Recl=recl)

    do ik = 1, kset%nkpt
        read(70, Rec=ik) nk0, ib0, nb0, vkl, &
             eqp(ib:nb,ik), eks(ib:nb,ik), &
             efqp, efks
        if ( abs(sum(vkl(:)-kset%vkl(:,ik))) > 1.d-6 ) then
            write(*,*) 'ERROR(readevalqp): Inconsistent k-points!'
            stop
        end if
    end do ! ik

    close(70)

end subroutine