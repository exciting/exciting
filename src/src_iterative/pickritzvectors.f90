! Copyright (C) 2015-2023 exciting team (Berlin and Riga)
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine pickritzvectors(n,segmentsize,evals,firstvec)
!> Selects a given number of the lowest eigenpairs that are larger than a predefined threshold value.
!> 
      use modmpi, only: terminate_mpi_env, mpiglobal
      use mod_APW_LO
      implicit none
      integer, intent(in) :: n             ! size of the subspace 
      integer, intent(in) :: segmentsize   ! number of eigenvalues to keep
      real(8), intent(in) :: evals(n)      ! list of eigenvalues
      integer, intent(out) :: firstvec     ! first eigenvalue that matches the requirements

      integer :: i
      real(8) :: ta,tb

      call timesec(ta)      
      i=1      
      do while ((i.le.n-segmentsize).and.(evals(i).lt.mine0))
       i=i+1
      enddo
      if (i.le.n-segmentsize+1) then
        firstvec=i
      else
        write(*,*) 'Eigenvalue filtering:'
        write(*,*) 'Too many eigenvalues under the permitted bound of ',mine0
        call terminate_mpi_env(mpiglobal) 
      endif 
      call timesec(tb)
!write(*,*) 'pickritzvectors', tb-ta

end subroutine pickritzvectors
