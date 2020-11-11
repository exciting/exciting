recursive subroutine sortidxColumn( nr, nc, a, lda, idx)
  ! Returns the index map that sorts the columns of an 
  ! nr times nc integer array in ascending order.
  ! Thereby, each column represents a numbers and the
  ! nr values in each row are treated as the digits of 
  ! this number.
  integer, intent( in)  :: nr, nc, lda
  integer, intent( in)  :: a(lda,*)
  integer, intent( out) :: idx(nc)

  integer :: i, j, k, l, srt(nc)
  integer :: tmp( nr-1, nc)

  call sortidxi( nc, a(1,1), lda, idx)
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
      i = j - l + 1
      do k = 1, l
        tmp(:,k) = a(2:nr,idx(i+k-1))
      end do
      call sortidxColumn( nr-1, l, tmp, nr-1, srt(1:l))
      idx(i:j) = idx(i+srt(1:l)-1)
      l = 1
    end if
  end do    
  return
end subroutine sortidxColumn

