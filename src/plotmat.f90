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

    !subroutine writemat( mat, xdim, ydim, fname, r)
    !  integer, intent( in) :: xdim, ydim, r
    !  complex(8), intent( in) :: mat( xdim, ydim)
    !  character(*), intent( in) :: fname

    !  integer :: recln

    !  inquire( iolength=recln) mat
    !  open( 50, file=trim( fname)//".mat", action='WRITE', access='DIRECT', recl=recln, position='append')
    !  write( 50) mat
    !  close( 50)
    !  return
    !end subroutine writemat

    !subroutine readmat( mat, xdim, ydim, fname, r)
    !  integer, intent( in) :: xdim, ydim, r
    !  complex(8), intent( out) :: mat( xdim, ydim)
    !  character(*), intent( in) :: fname

    !  integer :: recln

    !  inquire( iolength=recln) mat
    !  open( 50, file=trim( fname)//".mat", action='READ', access='DIRECT', status='OLD', recl=recln)
    !  read( 50, rec=r) mat
    !  close( 50)
    !  return
    !end subroutine readmat
    subroutine writemat( mat, xdim, ydim, fname)
      integer, intent( in) :: xdim, ydim
      complex(8), intent( in) :: mat( xdim, ydim)
      character(*), intent( in) :: fname

      integer :: x, y

      open( 50, file=trim( fname)//".mat", action='WRITE', form='FORMATTED', position='APPEND')
      do x = 1, xdim
        do y = 1, ydim
          write( 50, '(SP,2F23.16)') mat( x, y)
        end do
      end do
      write( 50, *)
      close( 50)
      return
    end subroutine writemat

    subroutine readmat( mat, xdim, ydim, fname)
      integer, intent( in) :: xdim, ydim
      complex(8), intent( out) :: mat( xdim, ydim)
      character(*), intent( in) :: fname

      integer :: x, y
      !logical, save :: rwnd = .true.

      !if( rwnd) then
        open( 50, file=trim( fname)//".mat", action='READ', form='FORMATTED', status='OLD', position='REWIND')
      !  rwnd = .false.
      !else
      !  open( 50, file=trim( fname)//".mat", action='READ', form='FORMATTED', status='OLD', position='ASIS')
      !end if
      do x = 1, xdim
        do y = 1, ydim
          read( 50, '(SP,2F23.16)') mat( x, y)
        end do
      end do
      read( 50, *)
      close( 50)
      return
    end subroutine readmat

    subroutine writematlab( mat, fname)
      complex(8), intent( in) :: mat(:,:)
      character(*), intent( in) :: fname

      integer :: sx, sy, x, y
      logical :: ml
      
      sx = size( mat, 1)
      sy = size( mat, 2)

      open( 50, file=trim( fname)//".mat", action='WRITE', form='FORMATTED')
      do x = 1, sx
        do y = 1, sy
          write( 50, '(SP,E23.16,E23.16,"i",5x)', advance='no') mat( x, y)
        end do
        write( 50, *)
      end do
      close( 50)
      return
    end subroutine writematlab
end module m_plotmat
