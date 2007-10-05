
module m_writeloss
  implicit none
contains

  subroutine writeloss(iq,w,loss,fn)
    use modmain
    use modtddft
    use m_getunit
    use m_tdwriteh
    implicit none
    ! arguments
    integer, intent(in) :: iq
    real(8), intent(in) :: w(:)
    real(8), intent(in) :: loss(:)
    character(*), intent(in) :: fn
    ! local variables
    character(*), parameter :: thisnam = 'writeloss'
    integer :: n1(1),n,iw

    if (any(shape(w).ne.shape(loss))) then
       write(unitout,'(a)') 'Error('//thisnam//'): input arrays have &
            &diffenrent shape'
       call terminate()
    end if

    n1=shape(w)
    n=n1(1)

    call getunit(unit1)
    open(unit1,file=trim(fn),action='write')
    ! write parameters as header to file
    call tdwriteh(unit1,iq)

!!$    ! write data to file
!!$    write(unit1,'(2g18.10)') (w(iw),loss(iw),iw=1,n)
    ! include dynamical structure factor
    !*** like in weissker2006
    write(unit1,'(3g18.10)') (w(iw)*escale,loss(iw),loss(iw)* &
         (gqc(1,iq)**2/(4.d0*pi**2*chgtot/omega)),iw=1,n)
!!$    ! *** like in Idoia thesis
!!$    write(unit1,'(3g18.10)') (w(iw),loss(iw),loss(iw)* &
!!$         gqc(1,iq)**2*omega/(twopi),iw=1,n)

    ! close files
    close(unit1)

  end subroutine writeloss

end module m_writeloss
