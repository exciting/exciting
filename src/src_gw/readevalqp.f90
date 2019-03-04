
subroutine readevalqp(fname)

    use modgw, only : kset, ibgw, nbgw, evalqp, evalks, eferqp, eferks
    implicit none
    character(*), intent(in) :: fname
    integer(4)    :: ik, nk, ib, nb
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
      
    inquire( IoLength=recl ) nk, ib, nb
    open(70, File=trim(fname), Action='READ', Form='UNFORMATTED', &
         Access='DIRECT', Recl=recl)
    read(70, Rec=1) nk, ib, nb
    close(70)

    ! Consistency check
    if ( nk /= kset%nkpt ) then
        write(*,*) 'ERROR(readevalqp): Inconsistent number of k-points!'
        write(*,*) 'nk = ', nk, ' nkpt = ', kset%nkpt
        stop
    end if
    if ( ib /= ibgw ) then
        write(*,*) 'ERROR(readevalqp): Inconsistent number of GW states!'
        write(*,*) 'ib = ', ib, ' ibgw = ', ibgw
        stop
    end if
    if ( nb /= nbgw ) then
        write(*,*) 'ERROR(readevalqp): Inconsistent number of GW states!'
        write(*,*) 'nb = ', nb, ' nbgw = ', nbgw
        stop
    end if

    inquire( IoLength=recl) nk, ib, nb, vkl, &
             evalqp(ibgw:nbgw,1), evalks(ibgw:nbgw,1), &
             eferqp, eferks
  
    open(70, File=trim(fname), Action='READ', Form='UNFORMATTED', &
         Access='DIRECT', Recl=recl)
 
    do ik = 1, nk
        read(70, Rec=ik) nk, ib, nb, vkl, &
             evalqp(ibgw:nbgw,ik), evalks(ibgw:nbgw,ik), &
             eferqp, eferks
        if ( abs(sum(vkl(:)-kset%vkl(:,ik))) > 1.d-6 ) then
            write(*,*) 'ERROR(readevalqp): Inconsistent k-points!'
            stop
        end if
    end do ! ik
    
    close(70)

end subroutine