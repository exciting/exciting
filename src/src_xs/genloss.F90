
module m_genloss
  implicit none
contains

  subroutine genloss(eps,loss)
    use modxs
    implicit none
    ! arguments
    complex(8), intent(in) :: eps(:)
    real(8), intent(out) :: loss(:)
    ! local variables
    character(*), parameter :: thisnam = 'genloss'

    if (any(shape(eps).ne.shape(loss))) then
       write(unitout,'(a)') 'Error('//thisnam//'): input and output arrays &
            &have diffenrent shape'
       call terminate
    end if

    ! loss function
    loss(:) = - aimag(1/eps(:))

  end subroutine genloss

end module m_genloss
