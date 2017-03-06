module m_sqrtzmat

  use modmpi
  use modscl
  use m_dhesolver
  use m_dzgemm

use m_writecmplxparts

  implicit none

  private

  complex(8), parameter :: zone = (1.0d0, 0.0d0)
  complex(8), parameter :: zzero = (0.0d0, 0.0d0)

  public :: sqrtzmat_hepd, sqrtdzmat_hepd

  contains

    subroutine sqrtzmat_hepd(hepdmat)

      complex(8), intent(inout) :: hepdmat(:,:)

      integer(4) :: m, n, i, j
      complex(8), allocatable :: evecs(:,:)
      real(8), allocatable :: evals(:)

      m = size(hepdmat,1)
      n = size(hepdmat,2)

      if(m /= n) then 
        write(*,*) "Error(sqrtmat_hepd): m /= n"
        call terminate
      end if

      allocate(evecs(m,m))
      allocate(evals(m))
      evals = 0.0d0

      ! Diagonalize hermitian matrix
      call hesolver(hepdmat, evals, evec=evecs)

      !write(*,*) "Writing amp evals"
      !call writecmplxparts("nd_amb_evals", revec=evals, veclen=size(evals))

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

    end subroutine sqrtzmat_hepd

    subroutine sqrtdzmat_hepd(hepdmat, binfo, eecs)

      type(dzmat), intent(inout) :: hepdmat
      type(blacsinfo), intent(in) :: binfo
      integer(4), intent(in), optional :: eecs

      integer(4) :: m, n, i, j, jg, clustersize
      type(dzmat) :: evecs
      real(8), allocatable :: evals(:)
      
      if(present(eecs)) then
        clustersize = eecs
      else
        clustersize = 3
      end if

      if(hepdmat%isdistributed == .false.) then
        call sqrtzmat_hepd(hepdmat%za)
        return
      end if

      if(hepdmat%context /= binfo%context) then
        write(*,*) "Error(sqrtdzmat_hepd): Wrong BLACS context."
        write(*,*) "  mat%context =", hepdmat%context
        write(*,*) "  binfo%context =", binfo%context
        call terminate
      end if

      m = hepdmat%nrows
      n = hepdmat%ncols

      if(m /= n) then 
        write(*,*) "Error(sqrtdzmat_hepd): m /= n"
        call terminate
      end if

      allocate(evals(m))
      evals = 0.0d0

      call new_dzmat(evecs, m, n, binfo, hepdmat%mblck, hepdmat%nblck)

      ! Diagonalize hermitian matrix
      call dhesolver(hepdmat, evals, binfo, evecs, eecs=clustersize)

      !if(mpiglobal%rank == 0) then 
      !  call writecmplxparts('amb_evals', revec=evals, veclen=size(evals))
      !end if

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
      do j = 1, evecs%ncols_loc
        jg = evecs%c2g(j)
        do i = 1, evecs%nrows_loc
          evecs%za(i,j) = evecs%za(i,j)*evals(jg)
        end do
      end do
      
      ! Construct square root matrix
      call dzgemm(evecs, evecs, hepdmat, transb='C')

      deallocate(evals)
      call del_dzmat(evecs)

    end subroutine sqrtdzmat_hepd

end module m_sqrtzmat
