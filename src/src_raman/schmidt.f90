subroutine schmidt(mat,n,m,nvo)
! orthogonalizes the m vectors in matrix mat(n,m)
! using a simple Schmidt orthogonalization
! the resulting nvo orthonormal vectors are returned
! again in mat(n,1:nvo)
implicit none
integer, intent(in) :: n,m
real(8), intent(inout) :: mat(n,m)
integer, intent(out) :: nvo
real(8) :: tempmat(n,m)
integer :: i,j
!
tempmat = 0.d0
nvo = 1
do i = 1,m
   tempmat(:,i) = mat(:,i)
   do j = 1,i-1
      tempmat(:,i) = tempmat(:,i) - dot_product(mat(:,i),tempmat(:,j))*tempmat(:,j)
   enddo
!  write(*,*) tempmat(:,i)
   mat(:,i) = norm(tempmat(:,i),n)*tempmat(:,i)
   if (any(abs(mat(:,i)) .gt. 1.d-8)) then
      do j = i-1,1,-1
         if(any(abs(mat(:,i)-mat(:,j)) .gt. 1.d-3)) then
            nvo = nvo + 1
            exit
         endif
      enddo
   endif
enddo
!
contains

function norm(vec,n)
! computes the norm or a vector
implicit none
integer, intent(in) :: n
real(8), intent(in) :: vec(n)
real(8) :: norm
integer :: i
!
norm = 0.d0
do i = 1,n
   norm = norm + vec(i)*vec(i)
enddo
!write(*,*) 'norm ',norm
if (norm .gt. 1.d-8) norm = 1.d0 / sqrt(norm)
return
end function norm

end subroutine schmidt
