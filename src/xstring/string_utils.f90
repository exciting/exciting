!> Routines for workong with/on strings.
module string_utils

  implicit none

  private
  public :: split, validate_filename

  contains

  !> Split a string into two at the first occurence of a delimiter. The delimiter is not part of the results.
  pure subroutine split(string, delimiter, first, second)
    !> String to be split.
    character(*), intent(in) :: string
    !> Separator.
    character(1), intent(in) :: delimiter
    !> First part of the separated string.
    character(:), allocatable, intent(out) :: first
    !> Second part of the separated string.
    character(:), allocatable, intent(out) :: second

    integer :: first_end, second_start

    first_end = scan(string, delimiter) - 1
    second_start = scan(string, delimiter) + 1

    first = string(1 : first_end)
    second = string(second_start : len(string))
  end subroutine split  

  !> Validate a file name by checking if the body contains only allowed characters (all upper and lower case 
  !> alphabetical and `_`) and if the file name has the expected ending.
  pure logical function validate_filename(filename, expected_ending)

    use constants, only: lower_case_alphabet_set, upper_case_alphabet_set, digit_set

    !> File name to be checked
    character(*), intent(in) :: filename
    !> Expected file ending. If the string starts with a period (`.`), it will be removed.
    character(*), intent(in) :: expected_ending

    character(*), parameter :: allowed_special_characters = '_-', delimiter = '.'

    character(:), allocatable :: filename_local, allowed_characters, body, ending, expected_ending_local

    logical :: body_is_valid, ending_is_valid

    ! If expected_ending starts with a period, remove it.
    expected_ending_local = expected_ending
    if (expected_ending(1:1) == '.') expected_ending_local = expected_ending(2 :)

    filename_local = trim(adjustl(filename))
    allowed_characters = lower_case_alphabet_set // upper_case_alphabet_set // digit_set // allowed_special_characters

    call split(filename_local, delimiter, body, ending)

    body_is_valid = (verify(body, allowed_characters) == 0)
    ending_is_valid = (verify(ending, allowed_characters) == 0)

    ending_is_valid = (ending_is_valid .and. (expected_ending_local == ending))

    validate_filename = (body_is_valid .and. ending_is_valid)
  end function validate_filename



end module string_utils