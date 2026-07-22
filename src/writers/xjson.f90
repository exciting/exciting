!> Convert numerical data to JSON complient string.
module xjson
  use precision, only: sp, dp

  implicit none
  private

  !> format for double precision numbers
  character(8), parameter :: FMT_REAL_DP = 'e26.16e3'
  !> format for single precision numbers
  character(8), parameter :: FMT_REAL_SP = 'e18.8e2'

  public :: to_json

  interface to_json
    procedure :: integer_tensor_to_json, real_dp_tensor_to_json, complex_dp_tensor_to_json
  end interface to_json

contains

  !> convert integer rank 0 tensor to JSON string
  function integer_scalar_to_json( a ) result( json )
    !> tensor
    integer, intent(in) :: a
    !> JSON string
    character(:), allocatable :: json

    character(32) :: buff
    
    write( buff, '(i32)' ) a
    json = trim( adjustl( buff ) )
  end function integer_scalar_to_json

  !> convert real(dp) rank 0 tensor to JSON string
  function real_dp_scalar_to_json( a, precision ) result( json )
    !> tensor
    real(dp), intent(in) :: a
    !> precision (default: `dp`)
    integer, optional, intent(in) :: precision
    !> JSON string
    character(:), allocatable :: json

    integer :: p = dp
    character(26) :: buff
    
    if (present( precision )) p = precision

    if (p == sp) then
      write( buff, '('//FMT_REAL_SP//')' ) a
    else
      write( buff, '('//FMT_REAL_DP//')' ) a
    end if
    json = trim( adjustl( buff ) )
  end function real_dp_scalar_to_json

  !> convert complex(dp) rank 0 tensor to JSON string
  function complex_dp_scalar_to_json( a, precision ) result( json )
    !> tensor
    complex(dp), intent(in) :: a
    !> precision (default: `dp`)
    integer, optional, intent(in) :: precision
    !> JSON string
    character(:), allocatable :: json

    json = '[' // real_dp_scalar_to_json( a%re, precision ) // ', ' // real_dp_scalar_to_json( a%im, precision ) // ']'
  end function complex_dp_scalar_to_json

  !> convert real(dp) tensor with rank <= 5 to JSON string
  recursive function integer_tensor_to_json( a ) result( json )
    !> tensor
    integer, intent(in) :: a(..)
    !> JSON string
    character(len=:), allocatable :: json
  
    integer :: i, n

    if (rank(a) == 0) then
      n = 1
    else
      n = size( a, dim=rank(a) )
    end if

    select rank (a)
      rank (0)
        json = integer_scalar_to_json( a )
      rank default
        json = '['
        select rank (a)
          rank (1)
            do i = 1, n-1
              json = json // integer_tensor_to_json( a(i) ) // ', '
            end do
            json = json // integer_tensor_to_json( a(n) )
          rank (2)
            do i = 1, n-1
              json = json // integer_tensor_to_json( a(:,i) ) // ', '
            end do
            json = json // integer_tensor_to_json( a(:,n) )
          rank (3)
            do i = 1, n-1
              json = json // integer_tensor_to_json( a(:,:,i) ) // ', '
            end do
            json = json // integer_tensor_to_json( a(:,:,n) )
          rank (4)
            do i = 1, n-1
              json = json // integer_tensor_to_json( a(:,:,:,i) ) // ', '
            end do
            json = json // integer_tensor_to_json( a(:,:,:,n) )
          rank (5)
            do i = 1, n-1
              json = json // integer_tensor_to_json( a(:,:,:,:,i) ) // ', '
            end do
            json = json // integer_tensor_to_json( a(:,:,:,:,n) )
        end select
        json = json // ']'
    end select
  end function integer_tensor_to_json

  !> convert real(dp) tensor with rank <= 5 to JSON string
  recursive function real_dp_tensor_to_json( a, precision ) result( json )
    !> tensor
    real(kind=dp), intent(in) :: a(..)
    !> precision (default: `dp`)
    integer, optional, intent(in) :: precision
    !> JSON string
    character(len=:), allocatable :: json
  
    integer :: i, n

    if (rank(a) == 0) then
      n = 1
    else
      n = size( a, dim=rank(a) )
    end if

    select rank (a)
      rank (0)
        json = real_dp_scalar_to_json( a, precision )
      rank default
        json = '['
        select rank (a)
          rank (1)
            do i = 1, n-1
              json = json // real_dp_tensor_to_json( a(i), precision ) // ', '
            end do
            json = json // real_dp_tensor_to_json( a(n), precision )
          rank (2)
            do i = 1, n-1
              json = json // real_dp_tensor_to_json( a(:,i), precision ) // ', '
            end do
            json = json // real_dp_tensor_to_json( a(:,n), precision )
          rank (3)
            do i = 1, n-1
              json = json // real_dp_tensor_to_json( a(:,:,i), precision ) // ', '
            end do
            json = json // real_dp_tensor_to_json( a(:,:,n), precision )
          rank (4)
            do i = 1, n-1
              json = json // real_dp_tensor_to_json( a(:,:,:,i), precision ) // ', '
            end do
            json = json // real_dp_tensor_to_json( a(:,:,:,n), precision )
          rank (5)
            do i = 1, n-1
              json = json // real_dp_tensor_to_json( a(:,:,:,:,i), precision ) // ', '
            end do
            json = json // real_dp_tensor_to_json( a(:,:,:,:,n), precision )
        end select
        json = json // ']'
    end select
  end function real_dp_tensor_to_json

  !> convert complex(dp) tensor with rank <= 5 to JSON string
  recursive function complex_dp_tensor_to_json( a, precision ) result( json )
    !> tensor
    complex(kind=dp), intent(in) :: a(..)
    !> precision (default: `dp`)
    integer, optional, intent(in) :: precision
    !> JSON string
    character(len=:), allocatable :: json
  
    integer :: i, n

    if (rank(a) == 0) then
      n = 1
    else
      n = size( a, dim=rank(a) )
    end if

    select rank (a)
      rank (0)
        json = complex_dp_scalar_to_json( a, precision )
      rank default
        json = '['
        select rank (a)
          rank (1)
            do i = 1, n-1
              json = json // complex_dp_tensor_to_json( a(i), precision ) // ', '
            end do
            json = json // complex_dp_tensor_to_json( a(n), precision )
          rank (2)
            do i = 1, n-1
              json = json // complex_dp_tensor_to_json( a(:,i), precision ) // ', '
            end do
            json = json // complex_dp_tensor_to_json( a(:,n), precision )
          rank (3)
            do i = 1, n-1
              json = json // complex_dp_tensor_to_json( a(:,:,i), precision ) // ', '
            end do
            json = json // complex_dp_tensor_to_json( a(:,:,n), precision )
          rank (4)
            do i = 1, n-1
              json = json // complex_dp_tensor_to_json( a(:,:,:,i), precision ) // ', '
            end do
            json = json // complex_dp_tensor_to_json( a(:,:,:,n), precision )
          rank (5)
            do i = 1, n-1
              json = json // complex_dp_tensor_to_json( a(:,:,:,:,i), precision ) // ', '
            end do
            json = json // complex_dp_tensor_to_json( a(:,:,:,:,n), precision )
        end select
        json = json // ']'
    end select
  end function complex_dp_tensor_to_json

end module xjson
