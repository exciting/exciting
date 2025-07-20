! Created by  on 24/03/2022.

!> Module for converting native fortran types to chars.

!TODO(Bene): Issue #128: implement matrix conversion for real/complex(sp), integer, and logical.
module to_char_conversion
  use precision, only: sp, dp, str_16, str_32, i32, long_int
  use unit_conversion, only: sp_to_char, dp_to_char

  implicit none

  private
  public :: to_char

  ! Default strings are chosen for yaml file generation, to be parsed in with python code.

  !> Default string for `.true.`
  character(*), parameter :: default_true = "TRUE"
  !> Default string for `.false.`
  character(*), parameter :: default_false = "FALSE"
  !> Default string for imaginary part indicator ('i.e.' 3.0+4.0i).
  character(*), parameter :: default_imag_indicator = "i"
  !> Default string for separating array elements
  character(*), parameter :: default_sep = ","
  !> Default strings for indicating the start of an array
  character(*), parameter :: default_arr_start = "["
  !> Default strings for indicating the end of an array
  character(*), parameter :: default_arr_end = "]"
  !> String to be returned in case of empty array (regardless of array's type)
  character(*), parameter :: empty_array_return = default_arr_start // default_sep // default_arr_end
  !> String to be returned in case of empty matrix (regardless of array's type)
  character(*), parameter :: empty_matrix_return = default_arr_start // empty_array_return // default_arr_end

  interface to_char
    procedure :: convert_logical, convert_integer_i32, convert_integer_long_int, &
            convert_real_sp, convert_real_dp, &
            convert_complex_sp, convert_complex_dp, &
            convert_vec_logical, convert_vec_integer_i32, convert_vec_integer_long_int, &
            convert_vec_real_sp, convert_vec_real_dp, convert_matr_real_dp, &
            convert_vec_complex_sp, convert_vec_complex_dp, convert_matr_complex_dp
  end interface to_char

  contains

  ! Basic type conversion
  !> Convert **logical** to **character**
  pure function convert_logical(bool) result(char)
    !> **logical** to convert
    logical, intent(in) :: bool

    character(:), allocatable :: char

    if (bool) then
      char = default_true
    else
      char = default_false
    end if
  end function convert_logical

  !> Convert **integer(i32)** to **character**
  pure function convert_integer_i32(i) result(char)
    !> **integer** to convert
    integer(i32), intent(in) :: i

    character(:), allocatable :: char

    character(str_16) :: char_

    write(char_, *) i
    char = trim(adjustl(char_))
  end function convert_integer_i32

  !> Convert **integer(long_int)** to **character**
  pure function convert_integer_long_int(i) result(char)
    !> **integer** to convert
    integer(long_int), intent(in) :: i

    character(:), allocatable :: char

    character(str_32) :: char_

    write(char_, *) i
    char = trim(adjustl(char_))
  end function convert_integer_long_int

  !> Convert **real(sp)** to **character**
  !> 
  !> The output format is given such that `13.352 --> "1.335200E+01"`
  pure function convert_real_sp(x) result(char)
    !> **real(sp)** to convert.
    real(sp), intent(in) :: x
    character(:), allocatable :: char

    character(str_16) :: char_

    write(char_, sp_to_char) x
    char = trim(adjustl(char_))
  end function convert_real_sp


  !> Convert **real(dp)** to **character**
  !>
  !> The output format is given such that `13.352 --> "1.33520000000000E+01"`.
  pure function convert_real_dp(d) result(char)
    !> **real(dp)** to convert.
    real(dp), intent(in) :: d
    character(:), allocatable :: char

    character(str_32) :: char_

    write(char_, dp_to_char) d
    char = trim(adjustl(char_))
  end function convert_real_dp


  !> Convert **complex(sp)** to **character**
  !>
  !> The output format is given such that `1.0 + 0.5 i --> "1.000000E+00+5.000000E-01i"`
  pure function convert_complex_sp(c) result(char)
    !> **complex(sp)** to convert.
    complex(sp), intent(in) :: c
    character(:), allocatable :: char

    character(:), allocatable :: rval, imval
    character(1) :: sign
    character(*), parameter :: imag_indicator = default_imag_indicator

    rval = convert_real_sp(real(c))
    imval = convert_real_sp(abs(imag(c)))
    sign = "+"
    if (imag(c) < 0) sign = "-"

    char = trim(adjustl(rval // sign // imval // imag_indicator))
  end function convert_complex_sp


  !> Convert **complex(dp)** to **character**
  !>
  !> The output format is given such that `1.0 + 0.5 i --> "1.00000000000000E+00+5.00000000000000E-01i"`
  pure function convert_complex_dp(z) result(char)
    !> **complex(dp)** to convert.
    complex(dp), intent(in) :: z

    character(:), allocatable :: char

    character(:), allocatable :: rval, imval
    character(1) :: sign
    character(*), parameter :: imag_indicator = default_imag_indicator

    rval = convert_real_dp(real(z))
    imval = convert_real_dp(abs(imag(z)))
    sign = "+"
    if (imag(z) < 0) sign = "-"

    char = trim(adjustl(rval // sign // imval // imag_indicator))
  end function convert_complex_dp


  ! Vector conversion
  !> Convert vector filled with **logical**s to **character**.
  pure function convert_vec_logical(vec) result(char)
    !> **logical** vector to convert.
    logical, intent(in) :: vec(:)
    character(:), allocatable :: char

    character(:), allocatable :: char_
    character(*), parameter :: sep = default_sep
    character(*), parameter :: arr_start = default_arr_start
    character(*), parameter :: arr_end = default_arr_end
    integer(i32) :: i, len_char, size_vec
    integer(i32), parameter :: max_len_var = max(len(default_true), len(default_false))

    size_vec = size(vec)
    if( size_vec == 0 ) then
      char = empty_array_return
    else
      len_char = len_vec_char(size_vec, max_len_var, len(sep), len(arr_start), len(arr_end))
      allocate(character(len_char) :: char_)
      char_ = adjustl(arr_start)
      do i = 1, size_vec-1
        char_ = adjustl(char_ // convert_logical(vec(i)) // sep)
      end do
      char = trim(adjustl(char_ // convert_logical(vec(size_vec)) // arr_end))
    end if
  end function convert_vec_logical


  !> Convert vector filled with **integer(i32)**s to **character**.
  pure function convert_vec_integer_i32(vec) result(char)
    !> **integer(i32)** vector to convert.
    integer(i32), intent(in) :: vec(:)
    character(:), allocatable :: char

    character(:), allocatable :: char_
    character(str_32) :: largest_int
    character(*), parameter :: sep = default_sep
    character(*), parameter :: arr_start = default_arr_start
    character(*), parameter :: arr_end = default_arr_end
    integer(i32) :: i, len_char, size_vec
    integer(i32) :: max_len_var

    size_vec = size(vec)
    if( size_vec == 0 ) then
      char = empty_array_return
    else
      write(largest_int,'(I0)') huge(0_i32)
      max_len_var = len_trim(largest_int)
      len_char = len_vec_char(size_vec, max_len_var, len(sep), len(arr_start), len(arr_end))
      allocate(character(len_char) :: char_)
      char_ = adjustl(arr_start)
      do i = 1, size_vec - 1
        char_ = adjustl(char_ // convert_integer_i32(vec(i)) // sep)
      end do
      char = trim(adjustl(char_ // convert_integer_i32(vec(size_vec)) // arr_end))
    end if
  end function convert_vec_integer_i32

  !> Convert vector filled with **integer(long_int)**s to **character**.
  pure function convert_vec_integer_long_int(vec) result(char)
    !> **integer(long_int)** vector to convert.
    integer(long_int), intent(in) :: vec(:)
    character(:), allocatable :: char

    character(:), allocatable :: char_
    character(str_32) :: largest_long_int
    character(*), parameter :: sep = default_sep
    character(*), parameter :: arr_start = default_arr_start
    character(*), parameter :: arr_end = default_arr_end
    integer(i32) :: i, len_char, size_vec
    integer(i32) :: max_len_var

    size_vec = size(vec)
    if( size_vec == 0 ) then
      char = empty_array_return
    else
      write(largest_long_int,'(I0)') huge(0_long_int)
      max_len_var = len_trim(largest_long_int)
      len_char = len_vec_char(size_vec, max_len_var, len(sep), len(arr_start), len(arr_end))
      allocate(character(len_char) :: char_)
      char_ = adjustl(arr_start)
      do i = 1, size_vec - 1
        char_ = adjustl(char_ // convert_integer_long_int(vec(i)) // sep)
      end do
      char = trim(adjustl(char_ // convert_integer_long_int(vec(size_vec)) // arr_end))
    end if
  end function convert_vec_integer_long_int


  !> Convert vector filled with **real(sp)**s to **character**
  pure function convert_vec_real_sp(vec) result(char)
    !> **real(sp)** vector to convert.
    real(sp), intent(in) :: vec(:)
    character(:), allocatable :: char

    character(:), allocatable :: char_
    character(*), parameter :: sep = default_sep
    character(*), parameter :: arr_start = default_arr_start
    character(*), parameter :: arr_end = default_arr_end
    integer(i32) :: i, len_char, size_vec
    integer(i32) :: max_len_var

    size_vec = size(vec)
    if( size_vec == 0 ) then
      char = empty_array_return
    else
      max_len_var = len_string_with_real_sp()
      len_char = len_vec_char(size_vec, max_len_var, len(sep), len(arr_start), len(arr_end))
      allocate(character(len_char) :: char_)
      char_ = adjustl(arr_start)
      do i = 1, size_vec - 1
        char_ = adjustl(char_ // convert_real_sp(vec(i)) // sep)
      end do
      char = trim(adjustl(char_ // convert_real_sp(vec(size_vec)) // arr_end))
    end if
  end function convert_vec_real_sp


  !> Convert vector filled with **real(dp)**s to **character**
  pure function convert_vec_real_dp(vec) result(char)
    !> **real(dp)** vector to convert.
    real(dp), intent(in) :: vec(:)
    character(:), allocatable :: char

    character(:), allocatable :: char_
    character(*), parameter :: sep = default_sep
    character(*), parameter :: arr_start = default_arr_start
    character(*), parameter :: arr_end = default_arr_end
    integer(i32) :: i, len_char, size_vec
    integer(i32) :: max_len_var

    size_vec = size(vec)
    if( size_vec == 0 ) then
      char = empty_array_return
    else
      max_len_var = len_string_with_real_dp()
      len_char = len_vec_char(size_vec, max_len_var, len(sep), len(arr_start), len(arr_end))
      allocate(character(len_char) :: char_)
      char_ = adjustl(arr_start)
      do i=1, size_vec - 1
        char_ = adjustl(char_ // convert_real_dp(vec(i)) // sep)
      end do
      char = trim(adjustl(char_ // convert_real_dp(vec(size_vec)) // arr_end))
    end if
  end function convert_vec_real_dp


  !> Convert vector filled with **complex(sp)**s to **character**
  pure function convert_vec_complex_sp(vec) result(char)
    !> **complex(sp)** vector to convert.
    complex(sp), intent(in) :: vec(:)
    character(:), allocatable :: char

    character(:), allocatable :: char_
    character(*), parameter :: sep = default_sep
    character(*), parameter :: arr_start = default_arr_start
    character(*), parameter :: arr_end = default_arr_end
    integer(i32) :: i, len_char, size_vec
    integer(i32) :: max_len_var

    size_vec = size(vec)
    if( size_vec == 0 ) then
      char = empty_array_return
    else
      max_len_var = 2*len_string_with_real_sp() + 2
      len_char = len_vec_char(size_vec, max_len_var, len(sep), len(arr_start), len(arr_end))
      allocate(character(len_char) :: char_)
      char_ = adjustl(arr_start)
      do i=1, size_vec - 1
        char_ = adjustl(char_ // convert_complex_sp(vec(i)) // sep)
      end do
      char = trim(adjustl(char_ // convert_complex_sp(vec(size_vec)) // arr_end))
    end if
  end function convert_vec_complex_sp


  !> Convert vector filled with **complex(dp)**s to **character**
  pure function convert_vec_complex_dp(vec) result(char)
    !> **complex(dp)** vector to convert.
    complex(dp), intent(in) :: vec(:)
    character(:), allocatable :: char

    character(:), allocatable :: char_
    character(*), parameter :: sep = default_sep
    character(*), parameter :: arr_start = default_arr_start
    character(*), parameter :: arr_end = default_arr_end
    integer(i32) :: i, len_char, size_vec
    integer(i32) :: max_len_var

    size_vec = size(vec)
    if( size_vec == 0 ) then
      char = empty_array_return
    else
      max_len_var = 2*len_string_with_real_dp() + 2
      len_char = len_vec_char(size_vec, max_len_var, len(sep), len(arr_start), len(arr_end))
      allocate(character(len_char) :: char_)
      char_ = adjustl(arr_start)
      do i=1, size_vec - 1
        char_ = adjustl(char_ // convert_complex_dp(vec(i)) // sep)
      end do
      char = trim(adjustl(char_ // convert_complex_dp(vec(size_vec)) // arr_end))
    end if
  end function convert_vec_complex_dp


  ! Matrix conversion
  !> Convert **real(dp)** matrix to character following column like order.
  pure function convert_matr_real_dp(matr) result(char)
    !> **real(dp)** matrix to convert.
    real(dp), intent(in) :: matr(:, :)
    character(:), allocatable :: char

    character(:), allocatable :: char_
    character(*), parameter :: sep = default_sep
    character(*), parameter :: arr_start = default_arr_start
    character(*), parameter :: arr_end = default_arr_end
    integer(i32) :: i, len_char, len_col_char, size_matr, shape_matr(2)
    integer(i32) :: max_len_var

    size_matr = size(matr)
    if( size_matr == 0 ) then
      char = empty_matrix_return
    else 
      max_len_var = len_string_with_real_dp()
      shape_matr = shape(matr)
      len_col_char = len_vec_char(shape_matr(2), max_len_var, len(sep), len(arr_start), len(arr_end))
      len_char = len_vec_char(shape_matr(1), len_col_char, len(sep), len(arr_start), len(arr_end))
      allocate(character(len_char) :: char_)
      char_ = adjustl(arr_start)
      do i = 1, shape_matr(2)-1
        char_ = adjustl(char_ // convert_vec_real_dp(matr(:, i)) // sep)
      end do
      char = trim(adjustl(char_ // convert_vec_real_dp(matr(:, shape_matr(2))) // arr_end))
    end if
  end function convert_matr_real_dp


  !> Convert **complex(dp)** matrix to character following column like order.
  pure function convert_matr_complex_dp(matr) result(char)
    !> **complex(dp)** matrix to convert.
    complex(dp), intent(in) :: matr(:, :)
    character(:), allocatable :: char

    character(:), allocatable :: char_
    character(*), parameter :: sep = default_sep
    character(*), parameter :: arr_start = default_arr_start
    character(*), parameter :: arr_end = default_arr_end
    integer(i32) :: i, len_char, len_col_char, size_matr, shape_matr(2)
    integer(i32) :: max_len_var

    size_matr = size(matr)
    if( size_matr == 0 ) then
      char = empty_matrix_return
    else 
      max_len_var = 2*len_string_with_real_dp() + 2
      shape_matr = shape(matr)
      len_col_char = len_vec_char(shape_matr(2), max_len_var, len(sep), len(arr_start), len(arr_end))
      len_char = len_vec_char(shape_matr(1), len_col_char, len(sep), len(arr_start), len(arr_end))
      allocate(character(len_char) :: char_)
      char_ = adjustl(arr_start)
      do i = 1, shape_matr(2)-1
        char_ = adjustl(char_ // convert_vec_complex_dp(matr(:, i)) // sep)
      end do
      char = trim(adjustl(char_ // convert_vec_complex_dp(matr(:, shape_matr(2))) // arr_end))
    end if
  end function convert_matr_complex_dp

  
  ! utils
  !> Calculate the number of single **character**s needed to convert a vector.
  pure integer(i32) function len_vec_char(size_vec, len_var, len_sep, len_arr_start, len_arr_end)
    !> Size of the vector
    integer(i32), intent(in) :: size_vec
    !> Number of **character**s needed to represent the **type** as string 
    integer(i32), intent(in) :: len_var
    !> Number of **character**s that represent the delimiter between the vectors elements (_i.e._ `"," -> 1`)
    integer(i32), intent(in) :: len_sep
    !> Number of **character**s that represent the start indicator of the vector (_i.e._ `"[" -> 1`)
    integer(i32), intent(in) :: len_arr_start
    !> Number of **character**s that represent the end indicator of the vector (_i.e._ `"]" -> 1`)
    integer(i32), intent(in) :: len_arr_end

    len_vec_char = len_arr_start + len_var * size_vec + len_sep * (size_vec - 1) + len_arr_end
  end function len_vec_char

  !> (private) Return the length required to store a real(dp) as string
  pure integer(i32) function len_string_with_real_dp()
    character(str_32) :: string_with_real_dp

    write(string_with_real_dp, dp_to_char) 0._dp
    len_string_with_real_dp = len_trim(string_with_real_dp)
  end function

  !> (private) Return the length required to store a real(sp) as string
  pure integer(i32) function len_string_with_real_sp()
    character(str_32) :: string_with_real_sp

    write(string_with_real_sp, sp_to_char) 0._sp
    len_string_with_real_sp = len_trim(string_with_real_sp)
  end function

end module to_char_conversion
