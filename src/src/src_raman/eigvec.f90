! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created July 2014 by Stefan Kontur
!
module m_eigvec
implicit none

contains

   subroutine eigvec(rotmat,vec,eigval)
      ! computes the eigenvector of a rotation matrix
      use raman_params, only : eps
      implicit none
      real(8), dimension(3,3), intent(in) :: rotmat
      real(8), dimension(3), intent(out) :: vec
      integer, intent(in) :: eigval
      ! local variables
      integer, dimension(3) :: ipiv
      real(8), dimension(3,3) :: mat
      logical, dimension(3) :: free
      integer :: i,j
      !
      ! form R - \lambda E
      mat = rotmat
      do i = 1,3
        mat(i,i) = mat(i,i) - dble(eigval)
      enddo
      !
      ! form row echelon form
      do i = 1,3
         call pivot(mat,ipiv)
         if (ipiv(i) .le. 3) then
            if (abs(mat(i,ipiv(i))) .gt. eps) then
               mat(i,:) = mat(i,:) / mat(i,ipiv(i))
            endif
         else
           cycle
         endif
         do j = i+1,3
            mat(j,:) = mat(j,:) - mat(j,ipiv(i))*mat(i,:)
         enddo
      enddo
      !
      ! determine free variables
      call pivot(mat,ipiv)
      free = .true.
      do i = 1,3
         if(ipiv(i) .le. 3) free(ipiv(i)) = .false.
      enddo
      !
      ! determine variables by backsubstitution
      vec = 0.
      do i = 3,1,-1
         if (free(i)) then
            vec(i) = 1.
            cycle
         else
            do j = i+1,3
               vec(i) = vec(i) - vec(j)*mat(i,j)
            enddo
         endif
      enddo
      !
      ! orientation conventions for reference
      if (abs(vec(3)) .lt. eps) then
         if (abs(vec(2)) .lt. eps) then
            if (vec(1) .lt. 0.d0) vec = -vec
         else
            if (vec(2) .lt. 0.d0) vec = -vec
         endif
      else
         if (vec(3) .lt. 0.d0) vec = -vec
      endif
   end subroutine eigvec
!
   subroutine pivot(mat,ipiv)
      use raman_params, only : eps
      implicit none
      real(8), dimension(3,3), intent(inout) :: mat
      integer, dimension(3), intent(out) :: ipiv
      real(8), dimension(3) :: temp
      integer :: i,j,k
      !
      ipiv = 4
      do i = 1,3
         do j = 1,3
            if (abs(mat(i,j)) .gt. eps) then
               ipiv(i) = j
               exit
            endif
         enddo
      enddo
      do i = 1,3
         do j = i+1,3
            if (ipiv(i) .gt. ipiv(j)) then
               temp(:) = mat(j,:)
               mat(j,:) = mat(i,:)
               mat(i,:) = temp(:)
               k = ipiv(j)
               ipiv(j) = ipiv(i)
               ipiv(i) = k
            endif
         enddo
      enddo
   end subroutine pivot

end module m_eigvec
