!> Module for sorting procedures.
!>
!> This module a) restores the original heapsort routine used throughout exciting prior to GITHASH
!> 57033f53f651c9bbdbbedc8d5512a69643d681ec
!> and b) retains cleaner (no GOTO) implementations for future use.
!>
!> # Legacy Routine `sortidx`
!>
!> The original implementation is retained because a tolerance is required to ensure a consistent
!> ordering for degenerate values. In the context of using sorting on G-vectors with degenerate magnitiudes,
!> numerical rounding errors/errors associated with the fast arithmetic in AVX, can lead to a
!> different ordering of G-vectors each time the routine is called, if a tolerance is not employed.
!> This behaviour is particularly difficult to catch with testing.
!>
!> One may be able to introduce a tolerance into new routines and utilise these where `sortidx` is already present,
!> however this would require testing with -O3 on various HPC architectures/build stacks, to confirm there
!> are no problems.
!>
module sorting
  use precision, only: dp

  implicit none

  private

  interface sort_index_1d
     module procedure :: sort_index_1d_int, sort_index_1d_real
  end interface sort_index_1d

  public :: sortidx, sort_index_1d, sort_index_2d

contains

    !>  Finds the permutation index `idx` which sorts the real array `a` in ascending order.
    !>  Uses the heapsort algorthim. Array `a` is preserved.
    !>
    !>   Created October 2002 (JKD)
    !>   Included tolerance eps, April 2006 (JKD)
  Subroutine sortidx (n, a, idx)
    Implicit None
    !>  Number of elements in array
    Integer, Intent (In) :: n
    !> Real array to sort
    Real (8), Intent (In) :: a (n)
    !>  Permutation index
    Integer, Intent (Out) :: idx (n)

    Integer :: i, j, k, l, m
    !> Relative tolerance factor for deciding if one number is smaller than another
    Real (8), Parameter :: eps = 1.d0 + 1.d-12
    If (n .Le. 0) Then
       Write (*,*)
       Write (*, '("Error(sortidx): n <= 0 : ", I8)') n
       Write (*,*)
       Stop
    End If
    Do i = 1, n
       idx (i) = i
    End Do
    If (n .Eq. 1) Return
    l = n / 2 + 1
    k = n
10  Continue
    If (l .Gt. 1) Then
       l = l - 1
       m = idx (l)
    Else
       m = idx (k)
       idx (k) = idx (1)
       k = k - 1
       If (k .Eq. 1) Then
          idx (1) = m
          Return
       End If
    End If
    i = l
    j = l + l
20  Continue
    If (j .Le. k) Then
       If (j .Lt. k) Then
          If (a(idx(j)) .Lt. eps * a(idx(j+1))) j = j + 1
       End If
       If (a(m) .Lt. eps * a(idx(j))) Then
          idx (i) = idx (j)
          i = j
          j = j + j
       Else
          j = k + 1
       End If
       Go To 20
    End If
    idx (i) = m
    Go To 10
  End Subroutine sortidx


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
    if( nr == 1) return
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
