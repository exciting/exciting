
module m_zaxpyc
  implicit none
contains

  subroutine zaxpyc(n,za,zx,incx,zy,incy)
    complex(8),intent(in) :: zx(:),za
    complex(8),intent(inout) :: zy(:)
    integer,intent(in) :: incx,incy,n
    ! automatic arrays
    complex(8) :: zt(size(zx))

    zt(:)=conjg(zx(:))
    call zaxpy(n,za,zt,incx,zy,incy)

  end subroutine zaxpyc


end module m_zaxpyc
