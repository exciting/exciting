!     BOP
!     
!     !ROUTINE: sqrtmatrix
!     
!     !INTERFACE:
      subroutine sqrtmatrix(matsiz,matrix,sqrmat,nsol)

!     !DESCRIPTION:
!     
!     This subroutine calculates the square root of a matrix.
!     nsol must be between 1 and 2**matsiz
!     For a given nsol, the solution is constructed in the following way:
!     - construct the binary representation of nsol: the maximum number (2**matsiz) would
!       be represented by 111...111 (matsiz digits).
!     - assign to each eigenvalue index the corresponding digit of the previously generated
!       binary number.
!     - if the digit corresponding to an eigenvalue is zero, then it chooses the square root
!       with positive real part, otherwise it chooses the sqare root with negative real part.
!     For example, if msiz = 3, the for nsol = 3 (011 in binary written with three digits)
!     square root of the matrix will be constructed with -sqrt(L1), -sqrt(L2), +sqrt(L3).
!     As the fortran sqrt function returns the solution with positive real part, then the
!     -sqrt(x) will have negative real part.
!     with a set of square roots of eigenvalues
! !REVISION HISTORY:
!   Created October 2015 (SR)
!EOP
!BOC
!     !USES:

      use modmain
      use modxs
      Use invert
!     
!     !INPUT/OUTPUT PARAMETERS:

      implicit none

      integer :: matsiz
      complex(8) :: matrix(matsiz,matsiz)
      complex(8) :: sqrmat(matsiz,matsiz)
      integer :: nsol
!     
!     !LOCAL VARIABLES:

      logical :: mtest
      integer(4) :: info
      integer(4) :: i, j, lwork,rwsize

      !real(8) :: tstart, tend

!     for diagonalization
      complex(8), allocatable :: tvec(:,:),itvec(:,:),itvec2(:,:)
      complex(8), allocatable :: tval(:),tmat(:,:)
      real(8), allocatable :: rwork(:)
      complex(8), allocatable :: work(:)


      mtest = .false.
!     
!     !REVISION HISTORY:
!     
!     Created April 2013 (SR)
!     
!     EOP
!     BOC      
      !call cpu_time(tstart)
      !if (tstart.lt.0.0d0) write(6,*)'warning, tstart < 0'
      
!----------------------------------
!     Calculate sqrtr of a matrix
!----------------------------------
!     
!     Set up the workspace for the diagonalization subroutine
!     
      allocate(tmat(matsiz,matsiz),itvec2(matsiz,matsiz))
      allocate(tvec(matsiz,matsiz),itvec(matsiz,matsiz))
      allocate(tval(matsiz))

      lwork=2*matsiz
      rwsize=2*matsiz
      allocate(work(lwork),rwork(rwsize))
!     
!     Diagonalize the matrix
!     
      do i = 1, matsiz
        do j = 1, matsiz
           tmat(i,j)=matrix(i,j)
        end do
      end do
      !tmat(:,:)=matrix(:,:)
!     query optimal lwork from zgeev
      call zgeev('n', 'v', matsiz, tmat, matsiz, tval, tvec, &
      & matsiz, tvec, matsiz, work, -1, rwork, info)
      lwork = work(1)
      deallocate(work)
      allocate(work(lwork))
!     diagonalize tmat
      call zgeev('n', 'v', matsiz, tmat, matsiz, tval, tvec, &
      & matsiz, tvec, matsiz, work, lwork, rwork, info)
      call errmsg(info.ne.0,'SQRTMATRIX',&
      & "Fail to diagonalize the matrix by ZGEEV !!!")

!     invert tvec
      Call zinvert_lapack (tvec, itvec)

      
!     
!     calculate the square root of the eigenvalues
!     
      do i = 1, matsiz
         tval(i)=(1-2*ibits(nsol-1,i-1,1))*sqrt(tval(i))
      enddo  

!     calculate tval*itvec

      do i = 1, matsiz
         itvec2(i,:)=tval(i)*itvec(i,:)
      end do
      call zgemm('n','n',matsiz,matsiz,matsiz,zone,tvec,matsiz,itvec2,&
      & matsiz,zzero,sqrmat,matsiz)
!     
!     deallocate local arrays
!     

      if(mtest) then
      do i = 1, matsiz
        do j = 1, matsiz
           write(*,*) "A(",i,",",j,") = ",matrix(i,j)
        end do
      end do
      write(*,*) "matsiz = ",matsiz,", lwork = ",lwork
      call zgemm('n','n',matsiz,matsiz,matsiz,zone,sqrmat,matsiz,sqrmat,&
      & matsiz,zzero,itvec2,matsiz)
      do i = 1, matsiz
        do j = 1, matsiz
            if (abs(itvec2(i,j)-matrix(i,j)) .gt. 1.d-8) then
                write(*,*) "|sqrt(A)**2-A| = ",abs(itvec2(i,j)-matrix(i,j)), "; for i,j = ",i,",",j
            end if
        end do
      end do
      end if

      deallocate(work,rwork)
      deallocate(tval,tmat,tvec,itvec,itvec2)

      !call cpu_time(tend)
      !if (tend.lt.0.0d0) write(6,*) 'warning, tend < 0'
      !call write_cputime(6,tend-tstart,'SQRTMATRIX')
      
      end subroutine sqrtmatrix
!     EOC 
