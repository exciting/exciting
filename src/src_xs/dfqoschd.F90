
module m_dfqoschd
  implicit none
contains

  subroutine dfqoschd(pou,puo,you,yuo)
    use modmain
    use modxs
    implicit none
    ! arguments
    complex(8), intent(in) :: pou(3),puo(3)
    complex(8), intent(out) :: you,yuo
    ! local variables
    integer :: oc,i,j
    real(8) :: s(3,3)
    complex(8) :: zt

    ! symmetrization matrix

    ! old stuff <= version 0.9.74
!!$    s(:,:)=0.5d0*(symdfq0(:,:)+transpose(symdfq0))

    s(:,:)=symdfq0(:,:)
    you=(0.d0,0.d0)
    yuo=(0.d0,0.d0)
    do i=1,3
       do j=1,3
          you=you+s(i,j)*pou(i)*conjg(pou(j))
          yuo=yuo+s(i,j)*puo(i)*conjg(puo(j))
       end do
    end do

  end subroutine dfqoschd

end module m_dfqoschd
