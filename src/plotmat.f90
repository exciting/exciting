module m_plotmat
  implicit none
  
  contains
    subroutine plotmat( mat, matlab)
      complex(8), intent( in) :: mat(:,:)
      logical, intent( in), optional :: matlab

      integer :: sx, sy, x, y
      logical :: ml
      
      sx = size( mat, 1)
      sy = size( mat, 2)
      !write(*,*) sx, sy

      if( present( matlab)) then
        if( matlab) then
          ml = .true.
        else 
          ml = .false.
        end if
      else
        ml = .false.
      end if
      
      if( ml) then
        write( *, '("M = [")', advance='no')
        do x = 1, sx-1
          do y = 1, sy-1
            write( *, '(SP,E23.16,E23.16,"*1i, ")', advance='no') mat( x, y)
          end do
          write( *, '(SP,E23.16,E23.16,"*1i; ")', advance='no') mat( x, y)
        end do
        do y = 1, sy-1
          write( *, '(SP,E23.16,E23.16,"*1i, ")', advance='no') mat( x, y)
        end do
        write( *, '(SP,E23.16,E23.16,"*1i];")', advance='no') mat( x, y)
      else
        do x = 1, sx
          do y = 1, sy
            write( *, '(SP,F13.6,F13.6,"i",5x)', advance='no') mat( x, y)
          end do
          write(*,*)
        end do
      end if
    end subroutine plotmat
end module m_plotmat
