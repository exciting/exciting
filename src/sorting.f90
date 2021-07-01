!> Module for sorting procedures
module sorting
  use precision, only: dp

  implicit none

  private

  interface sort_index_1d
    module procedure :: sort_index_1d_int, sort_index_1d_real
  end interface

  public :: sort_index_1d, sort_index_2d

  contains

    !> Returns the index map that sorts the elements of a 1d integer array
    !> in ascending order. Employs the Heapsort algorithm.
    function sort_index_1d_int( n, a, inc) result( idx)
      !> number of elements in the array
      integer, intent( in) :: n
      !> the array to be sorted
      integer, intent( in) :: a(*)
      !> increment (consider only every inc-th element)
      integer, optional, intent( in) :: inc
      !> index map that sorts the array
      integer :: idx(n)

      integer :: incr
      integer :: i, j, k, l, m

      incr = 1
      if( present( inc)) incr = inc

      do i = 1, n
        idx(i) = i
      end do
      if( n == 1) return
      l = n/2 + 1
      k = n
      do while( .true.)
        if( l > 1) then
          l = l - 1
          m = idx(l)
        else
          m = idx(k)
          idx(k) = idx(1)
          k = k - 1
          if( k == 1) then
            idx(1) = m
            exit
          end if
        end if
        i = l
        j = l + l
        do while( j <= k)
          if( j < k) then
            if( a( incr*(idx(j)-1)+1) < a( incr*(idx(j+1)-1)+1)) j = j + 1
          end if
          if( a( incr*(m-1)+1) < a( incr*(idx(j)-1)+1)) then
            idx(i) = idx(j)
            i = j
            j = j + j
          else
            j = k + 1
          end if
        end do
        idx(i) = m
      end do
    end function sort_index_1d_int

    !> Returns the index map that sorts the elements of a 1d real array
    !> in ascending order. Employs the Heapsort algorithm
    function sort_index_1d_real( n, a, inc) result( idx)
      !> number of elements in the array
      integer, intent( in) :: n
      !> the array to be sorted
      real(dp), intent( in) :: a(*)
      !> increment (consider only every inc-th element)
      integer, optional, intent( in) :: inc
      !> index map that sorts the array
      integer :: idx(n)

      integer :: incr
      integer :: i, j, k, l, m

      incr = 1
      if( present( inc)) incr = inc

      do i = 1, n
        idx(i) = i
      end do
      if( n == 1) return
      l = n/2 + 1
      k = n
      do while( .true.)
        if( l > 1) then
          l = l - 1
          m = idx(l)
        else
          m = idx(k)
          idx(k) = idx(1)
          k = k - 1
          if( k == 1) then
            idx(1) = m
            exit
          end if
        end if
        i = l
        j = l + l
        do while( j <= k)
          if( j < k) then
            if( a( incr*(idx(j)-1)+1) < a( incr*(idx(j+1)-1)+1)) j = j + 1
          end if
          if( a( incr*(m-1)+1) < a( incr*(idx(j)-1)+1)) then
            idx(i) = idx(j)
            i = j
            j = j + j
          else
            j = k + 1
          end if
        end do
        idx(i) = m
      end do
    end function sort_index_1d_real

    !> Returns the index map that sorts the columns of an 
    !> $n_r \times n_c$ integer array in ascending order.
    !> Thereby, each column represents a numbers and the
    !> $nr$ values in each row are treated as the digits of 
    !> this number. 
    recursive function sort_index_2d( nr, nc, a, lda) result( idx)
      !> number of rows in the array
      integer, intent( in) :: nr
      !> number of columns in the array
      integer, intent( in) :: nc
      !> leading dimension of the array as allocated in the calling scope
      integer, intent( in) :: lda
      !> the integer array
      integer, intent( in) :: a(lda,*)
      !> index map that sorts the columns of the array
      integer :: idx(nc)
    
      integer :: i, j, k, l
      integer, allocatable :: tmp(:,:), srt(:)
    
      idx = sort_index_1d_int( nc, a(1,1), inc=lda)
      if( nr .eq. 1) return
      l = 1
      do j = 1, nc
        if( j < nc) then
          if( a(1,idx(j)) == a(1,idx(j+1))) then
            l = l + 1
            cycle
          end if
        end if
        if( l > 1) then
          allocate( tmp( nr-1, l), srt(l))
          i = j - l + 1
          do k = 1, l
            tmp(:,k) = a(2:nr,idx(i+k-1))
          end do
          srt = sort_index_2d( nr-1, l, tmp, nr-1)
          idx(i:j) = idx(i+srt-1)
          l = 1
          deallocate( tmp, srt)
        end if
      end do    
      return
    end function sort_index_2d

end module sorting
