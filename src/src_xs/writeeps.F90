
module m_writeeps
  implicit none
contains

  subroutine writeeps(iq,w,eps,fn)
    use modxs
    use m_getunit
    use m_tdwriteh
    implicit none
    ! arguments
    integer, intent(in) :: iq
    real(8), intent(in) :: w(:)
    complex(8), intent(in) :: eps(:)
    character(*), intent(in) :: fn
    ! local variables
    character(*), parameter :: thisnam = 'writeeps'
    integer :: n1(1),n,iw
    real(8), allocatable :: imeps(:),kkeps(:)

    if (any(shape(w).ne.shape(eps))) then
       write(unitout,'(a)') 'Error('//thisnam//'): input arrays have &
            &diffenrent shape'
       call terminate()
    end if

    n1=shape(w)
    n=n1(1)

    allocate(imeps(n),kkeps(n))

    ! Kramers-Kronig transform imaginary part to obtain real part
    imeps(:)=aimag(eps(:))
!!    call kramkron(1,1,1.d-8,n,w,imeps,kkeps)

    call getunit(unit1)
    open(unit1,file=trim(fn),action='write')
    ! write parameters as header to file
    call tdwriteh(unit1,iq)
    ! write data to file
!!!    write(unit1,'(3g18.10)') (w(iw)*escale,eps(iw),iw=1,n)
    write(unit1,'(4g18.10)') (w(iw)*escale,eps(iw),kkeps(iw),iw=1,n)
    ! close files
    close(unit1)

    deallocate(imeps,kkeps)

  end subroutine writeeps

end module m_writeeps
