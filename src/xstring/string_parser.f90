module string_parser


  contains 

  subroutine split_string(input_string, delimiter, words)
    implicit none
    character(len=*), intent(in)  :: input_string
    character(len=1), intent(in)  :: delimiter
    character(len=256), allocatable :: words(:)

    integer :: start, end_pos, count, str_len
    character(len=100) :: temp_words(100)  ! Temporary fixed-size array

    ! Initialize
    count = 0
    start = 1
    str_len = len_trim(input_string)

    ! Loop to extract words
    do
        end_pos = scan(input_string(start:), delimiter)
        if (end_pos == 0) then
            count = count + 1
            temp_words(count) = trim(input_string(start:))  ! Last word
            exit
        else
            count = count + 1
            temp_words(count) = trim(input_string(start:start + end_pos - 2))
            start = start + end_pos  ! Move past delimiter
        end if
    end do

    ! Allocate and copy words to final array
    allocate(words(count))
    words(:count) = temp_words(:count)

  end subroutine split_string 

end module string_parser 