recursive subroutine sortidxColumn( nr, nc, a, idx)
  ! Returns the index map that sorts the columns of an 
  ! nr times nc integer array in ascending order.
  ! Thereby, each column represents a numbers and the
  ! nr values in each row are treated as the digits of 
  ! this number.
  integer, intent( in)  :: nr, nc
  integer, intent( in)  :: a(nr,*)
  integer, intent( out) :: idx(nc)

  integer :: i, j, k, l, srt( nc)
  integer :: tmp( nr-1, nc)

  call sortidx( nc, dble( a(1,1:nc)), idx)
  if( nr .eq. 1) return
  j = 1
  do i = 1, nc-1
    if( (a( 1, idx( i)) .ne. a( 1, idx( i+1))) .or. (i+1 .eq. nc)) then
      k = i
      if( (a( 1, idx( i)) .eq. a( 1, idx( i+1))) .and. (i+1 .eq. nc)) k = nc
      do l = j, k
        tmp(:,l) = a( 2:nr, idx(j))
      end do
      call sortidxColumn( nr-1, k-j+1, tmp, srt( j:k))
      idx( j:k) = idx( j-1+srt( j:k))
      j = i+1
    end if
  end do    
  return
end subroutine sortidxColumn

