module m_sqrtmat

  use modmpi
  use m_hesolver
!  use m_writecmplxparts

  implicit none

  contains

    subroutine sqrtmat_hepd(hepdmat)

      complex(8), intent(inout) :: hepdmat(:,:)

      integer(4) :: m, n, i, j
      complex(8), allocatable :: evecs(:,:)
      real(8), allocatable :: evals(:)
      complex(8), parameter :: zone = (1.0d0, 0.0d0)
      complex(8), parameter :: zzero = (0.0d0, 0.0d0)

      m = size(hepdmat,1)
      n = size(hepdmat,2)

      if(m /= n) then 
        write(*,*) "Error(sqrtmat_hepd): m /= n"
        call terminate
      end if

      allocate(evecs(m,m))
      allocate(evals(m))

      ! Diagonalize hermitian matrix
      call hesolver(hepdmat, evecs, evals)

      ! Take square root of eigenvalues
      if(any(evals < 0.0d0)) then 
        write(*,*) "Error(sqrtmat_hepd): Matrix is not positive definit"
        write(*,'(E10.3)') evals
        call terminate
      else
        ! Take D^{1/4} for: A^1/2 = Q D^1/2 Q^H = (Q D^1/4)*(Q D^1/4)^H
        evals = sqrt(sqrt(evals))
      end if

      ! Make Q D^1/4
      do j = 1, m
        do i = 1, m
          evecs(i,j) = evecs(i,j)*evals(j)
        end do
      end do
      
      ! Construct square root matrix
      call zgemm('N','C', m, m, m, zone, evecs, m,&
        & evecs, m, zzero, hepdmat, m)

      deallocate(evecs)
      deallocate(evals)

    end subroutine sqrtmat_hepd

end module m_sqrtmat
