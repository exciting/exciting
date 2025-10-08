! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Module for setting up the eigensystem
! it is designed in a way that all other subroutines
! dealing with setting up and solving the system can acsess the
! data transparently allowing to choose from different datatypes
! more easily
module  modfvsystem
  use precision, only: i32, dp, sp, str_256
  use constants, only: zzero, zone
#include "offload.fpp"
  implicit none

  type HermitianMatrix
      integer(i32) :: rank
      logical :: packed, ludecomposed, sp
      integer(i32), pointer :: ipiv (:)
      complex(dp), pointer, contiguous :: za (:, :)
      complex(dp), pointer :: zap (:)
      complex(dp), pointer :: eigenvalues(:)
      complex(sp), pointer :: ca (:, :), cap (:)
      complex(sp), pointer :: eigenvaluessp(:)

  end type HermitianMatrix
  !
  type evsystem
      type (HermitianMatrix) :: hamilton, overlap
      complex(dp), allocatable :: apwi(:,:,:)
  end type evsystem

  integer(i32), parameter :: LUdecomp = 1
  integer(i32), parameter :: LDLdecomp = 2
  integer(i32), parameter :: LLdecomp = 3 ! Cholesky decomposition
  integer(i32), parameter :: Diagdecomp = 4 ! decomposition through diagonalization
  integer(i32), parameter :: InverseDiag = 5 ! inverse is approximated by a diagonal matrix 
  integer(i32), parameter :: InvertOnce = 6 ! matrix inversion is done just once, result is reused 

contains
  !
  !
  subroutine newmatrix (self, packed, rank)
      implicit none
      type (HermitianMatrix), intent (inout) :: self
      logical, intent (in) :: packed
      integer(i32), intent (in) :: rank
      self%rank = rank
      self%packed = packed
      self%ludecomposed = .false.
      nullify(self%ca)
      nullify(self%cap)
      nullify(self%za)
      nullify(self%zap)
      nullify(self%ipiv)
      if (self%sp) then
        if (packed) then
          allocate (self%cap(rank*(rank+1)/2))
          self%zap = 0.0
        else
          allocate (self%ca(rank, rank))
          self%ca = 0.0
        endif
      else
        if (packed) then
          allocate (self%zap(rank*(rank+1)/2))
          self%zap = 0.0
        else
          allocate (self%za(rank, rank))
          self%za = 0.0
        end if
      endif
  end subroutine newmatrix
  !
  !
  subroutine deletematrix (self)
      implicit none
      type (HermitianMatrix), intent (inout) :: self
      if (associated(self%zap)) deallocate (self%zap)
      if (associated(self%za))  deallocate (self%za)
      if (associated(self%cap)) deallocate (self%cap)
      if (associated(self%ca))  deallocate (self%ca)  
      if (associated(self%ipiv)) deallocate (self%ipiv)
  end subroutine deletematrix
  !
  !
  subroutine newsystem (self, packed, rank)
      implicit none
      type (evsystem), intent (out) :: self
      logical, intent (in) :: packed
      integer(i32), intent (in) :: rank
      self%hamilton%sp=.false.
      self%overlap%sp=.false.
      call newmatrix (self%hamilton, packed, rank)
      call newmatrix (self%overlap, packed, rank)
  end subroutine newsystem
  !
  !
  subroutine deletesystem (self)
      implicit none
      type (evsystem), intent (inout) :: self
      call deletematrix (self%hamilton)
      call deletematrix (self%overlap)
  end subroutine deletesystem
  !
  !
  subroutine Hermitianmatrix_rank2update (self, n, alpha, x, y) 
      implicit none
      type (HermitianMatrix), intent (inout) :: self
      integer(i32), intent (in) :: n
      complex(dp), intent (in) :: alpha, x (:), y (:)
      !
      if (self%packed) then
          call ZHPR2 ('U', n, alpha, x, 1, y, 1, self%zap)
      else
          call ZHER2 ('U', n, alpha, x, 1, y, 1, self%za, self%rank)
      end if
  end subroutine Hermitianmatrix_rank2update
  !
  !
  subroutine Hermitianmatrix_indexedupdate (self, i, j, z)
      implicit none
      type (HermitianMatrix), intent (inout) :: self
      integer(i32) :: i, j
      complex(dp) :: z
      integer(i32) :: ipx
      if (self%packed .eqv. .true.) then
          ipx = ((i-1)*i) / 2 + j
          self%zap (ipx) = self%zap(ipx) + z
      else
          if (j .le. i) then
              self%za (j, i) = self%za(j, i) + z
          else
              write (*,*) "warning lower part of hamilton updated"
          end if
      end if
      return
  end subroutine Hermitianmatrix_indexedupdate
  !
  !
  subroutine Hermitianmatrixvector (self, alpha, vin, beta, vout)
#ifdef USEOMP
      use omp_lib
#endif
      implicit none
      type (HermitianMatrix), intent (inout) :: self
      complex(dp), intent (in) :: alpha, beta
      complex(dp), intent (inout) :: vin (:)
      complex(dp), intent (inout) :: vout (:)
      integer(i32) ::  nthreads, whichthread, nst, nfin,bandsize, i
      complex(dp), allocatable :: outcome(:)
      !
      if (self%packed .eqv. .true.) then
          call zhpmv ("U", self%rank, alpha, self%zap, vin, 1, beta, &
          vout, 1)
      else
#ifdef USEOMP
        vout=beta*vout           
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(nthreads,whichthread,bandsize,nst,nfin,outcome) SHARED(self,alpha,vin,beta,vout)
        allocate(outcome(self%rank))
        nthreads=omp_get_num_threads()
        whichthread=omp_get_thread_num()
        bandsize=self%rank/nthreads
        nst=1+whichthread*bandsize
        if (whichthread+1.eq.nthreads) then
          nfin=self%rank
        else
          nfin=nst+bandsize-1
        end if

        call zgemv ("N", self%rank, nfin-nst+1,alpha,self%za(1,nst),self%rank,vin(nst),1,zzero,outcome,1)
        do i= 0, nthreads - 1
          if (i == whichthread) vout=vout+outcome
!$OMP BARRIER
        end do
        deallocate(outcome)
!$OMP END PARALLEL 
#else
        call zhemv ("U", self%rank, alpha, self%za, self%rank, vin, 1, beta, vout, 1)
#endif
      end if
  end subroutine Hermitianmatrixvector
  !
  !
  function ispacked (self)
      implicit none
      logical :: ispacked
      type (HermitianMatrix) :: self
      ispacked = self%packed
  end function ispacked
  !
  !
  function getrank (self)
      implicit none
      integer(i32) :: getrank
      type (HermitianMatrix) :: self
      getrank = self%rank
  end function getrank
  !
  !
  subroutine HermitianMatrixMatrix(self,zm1,zm2,ldzm,naa,ngp)
    implicit none
    type (HermitianMatrix), intent(inout) :: self
    complex(dp),intent(in) :: zm1(:,:),zm2(:,:)
    integer(i32),intent(in) :: ldzm,ngp,naa

    ! ZGEMM  performs one of the matrix-matrix operations
    !        C := alpha*op( A )*op( B ) + beta*C,
    call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
               'N', &           ! TRANSB = 'N'  op( B ) = B.
               ngp, &           ! M ... rows of op( A ) = rows of C
               ngp, &           ! N ... cols of op( B ) = cols of C
               naa, &           ! K ... cols of op( A ) = rows of op( B )
               zone, &          ! alpha
               zm1, &           ! A
               ldzm,&           ! LDA ... leading dimension of A
               zm2, &           ! B
               ldzm, &          ! LDB ... leading dimension of B
               zone, &          ! beta
               self%za(1,1), &  ! C
               self%rank &      ! LDC ... leading dimension of C
              )

  end subroutine HermitianMatrixMatrix
  !
  !
  subroutine HermitianmatrixInvert (self,Method)
      implicit none
      type (HermitianMatrix) :: self
      integer(i32), intent(inout) :: Method
      integer(i32) :: lwork, info
      complex(dp), allocatable :: work(:)
      complex(dp) :: worksize
      integer(i32) ::  i,j
     
      if ( ispacked(self)) then
        write(*,*) 'Packed matrices are not supported'
        stop
      endif
      if ( self%ludecomposed) then
        if (self%sp) then
          write(*,*) 'single precision not implemented yet (HermitianmatrixInvert)'
          stop
        else
          select case (Method)
          case (LUdecomp)
            call zgetri (self%rank, &            ! matrix size
                         self%za, &              ! matrix itself
                         self%rank, &            ! leading dimension
                         self%ipiv, &            ! factorization pivots
                         worksize, &                 ! workspace
                         -1, &                ! size of the workspace
                         info &                  ! error message
                        )
            lwork=int(worksize)
            allocate(work(lwork))
            call zgetri (self%rank, &            ! matrix size
                         self%za, &              ! matrix itself
                         self%rank, &            ! leading dimension
                         self%ipiv, &            ! factorization pivots
                         work, &                 ! workspace
                         lwork, &                ! size of the workspace
                         info &                  ! error message
                        )  
            deallocate(work)
          case (LDLdecomp)
            allocate(work(self%rank))
            call zhetri ('U', &                  ! upper or lower
                         self%rank, &            ! matrix size
                         self%za, &              ! matrix itself
                         self%rank, &            ! leading dimension
                         self%ipiv, &            ! factorization pivots
                         work, &                 ! workspace
                         info &                  ! error message
                        )
            deallocate(work)
          case (LLdecomp)
            call zpotri ('U', &                  ! upper or lower
                         self%rank, &            ! matrix size
                         self%za, &              ! matrix itself
                         self%rank, &            ! leading dimension
                         info &                  ! error message
                        )
          end select
        endif
      endif
      if (info /= 0) then
        write(*,*) 'INFO (HermitianmatrixInvert) =', info
        stop
      endif
      if ((Method == LLDecomp).or.(Method == LDLDecomp)) then
! fill the lower triangle too
        do i=2,self%rank
          do j=1,i-1
            self%za(i,j)=conjg(self%za(j,i))
          enddo
        enddo
      endif
  end subroutine HermitianmatrixInvert
  !
  !
  subroutine HermitianmatrixFactorize (self,Method)
      implicit none
      type (HermitianMatrix) :: self
      integer(i32), intent(inout) :: Method
      if ( ispacked(self)) then
        write(*,*) 'Packed matrices are not supported'
        stop
      endif
      if ( .not. self%ludecomposed) then
        select case (Method)
        case (LUdecomp)
          call HermitianmatrixLU (self)
        case (LDLdecomp)
          call HermitianmatrixLDL (self)
        case (LLdecomp)
          call HermitianmatrixLL (self)
          if (associated(self%ipiv)) Method=LDLdecomp
        end select
      self%ludecomposed = .true.
      endif
  end subroutine HermitianmatrixFactorize
  !
  !
  subroutine HermitianmatrixLU (self)
      implicit none
      type (HermitianMatrix) :: self
      integer(i32) :: info
      allocate (self%ipiv(self%rank))
      if (self%sp) then
        call CGETRF (self%rank, self%rank, self%ca, self%rank, self%ipiv, info)
      else
        call ZGETRF (self%rank, self%rank, self%za, self%rank, self%ipiv, info)
      endif
      if (info .ne. 0) then
        write (*,*) "error in iterativearpacksecequn  HermitianmatrixLU "                , info
        stop
      end if
  end subroutine HermitianmatrixLU
  !
  !
  subroutine HermitianmatrixLDL (self)
      implicit none
      type (HermitianMatrix) :: self
      integer(i32) :: info
      complex(dp), allocatable :: zwork(:)
      complex(sp), allocatable :: cwork(:)
      allocate (self%ipiv(self%rank))
      if (self%sp) then
        allocate(cwork(64*self%rank))
        call chetrf('U', &                      ! upper or lower part
                     self%rank, &               ! size of matrix
                     self%ca, &                 ! matrix
                     self%rank, &               ! leading dimension
                     self%ipiv,&                ! pivot indices
                     cwork,&                     ! work
                     64*self%rank, &            ! work size
                     info &                     ! error message
                    )
        deallocate(cwork)
      else
        allocate(zwork(64*self%rank))
        call zhetrf('U', &                      ! upper or lower part
                     self%rank, &               ! size of matrix
                     self%za, &                 ! matrix
                     self%rank, &               ! leading dimension
                     self%ipiv,&                ! pivot indices
                     zwork,&                     ! work
                     64*self%rank, &            ! work size
                     info &                     ! error message
                    )
        deallocate(zwork)
      endif
      if (info /= 0) then
        write (*,*) "error in iterativearpacksecequn  HermitianmatrixLDL "                , info
        stop
      end if
  end subroutine HermitianmatrixLDL 
  !
  !
  subroutine HermitianmatrixLL (self) !Cholesky decomposition
      implicit none
      type (HermitianMatrix) :: self
      integer(i32) :: info
      complex(dp), allocatable :: zwork(:,:)
      complex(sp), allocatable :: cwork(:,:)
      if (self%sp) then
        allocate(cwork(self%rank,self%rank))
        cwork=self%ca
        call cpotrf('U',&                       ! upper or lower part
                     self%rank, &               ! size of matrix
                     self%ca, &                 ! matrix
                     self%rank, &               ! leading dimension
                     info &                     ! error message
                    )
      else
        allocate(zwork(self%rank,self%rank))
        zwork=self%za
        call zpotrf('U',&                       ! upper or lower part
                     self%rank, &               ! size of matrix
                     self%za, &                 ! matrix
                     self%rank, &               ! leading dimension
                     info &                     ! error message
                    )
      endif
      if (info /= 0) then
        if (self%sp) self%ca=cwork
        if (.not. self%sp) self%za=zwork
        call HermitianmatrixLDL (self)
      endif
      if (self%sp) then
        deallocate(cwork)
      else
        deallocate(zwork)
      endif
  end subroutine HermitianmatrixLL
  !
  !
  subroutine Hermitianmatrixlinsolve (self, b, method)
      implicit none
      type (HermitianMatrix) :: self
      complex(dp), intent (inout) :: b (:)
      integer(i32), intent (in) :: method
      integer(i32) :: info
      complex(dp), allocatable :: outcome(:)

      if ( ispacked(self)) then
        write(*,*) 'Packed matrices are not supported'
        stop
      endif        
      if (self%ludecomposed) then
        select case (method)
        case (LUdecomp)
          call ZGETRS ('N', self%rank, 1, self%za, self%rank, self%ipiv, b, self%rank, info)
        case (LDLdecomp)
          call zhetrs('U', &                      ! upper or lower part
                       self%rank, &                       ! size
                       1, &                       ! number of right-hand sides
                       self%za, &      ! factorized matrix
                       self%rank, &                       ! leading dimension
                       self%ipiv, &    ! pivoting indices
                       b, &                    ! right-hand side / solution
                       self%rank, &                       ! leading dimension
                       info &                     ! error message
                     )
        case (LLdecomp)
          call zpotrs('U', &                      ! upper or lower part
                       self%rank,  &                      ! size
                       1,  &                      ! number of right-hand sides
                       self%za, &      ! factorized matrix
                       self%rank, &                       ! leading dimension
                       b, &                    ! right-hand side / solution
                       self%rank, &                       ! leading dimension
                       info &                     ! error message
                     )
        case (InvertOnce) 
          info=0
          allocate(outcome(self%rank))
          call Hermitianmatrixvector (self, zone, b, zzero, outcome)
          b=outcome
          deallocate(outcome)
        end select
          if (info .ne. 0) then
              write (*,*) "error in iterativearpacksecequn Hermitianmatrixlinsolve "                , info
              stop
          end if
      end if
  end subroutine Hermitianmatrixlinsolve
  !
  !
  subroutine HermitianMatrixAXPY (alpha, x, y)
      implicit none
      complex(dp) :: alpha
      type (HermitianMatrix), intent(in)    :: x
      type (HermitianMatrix), intent(inout) :: y
      integer(i32) :: mysize
      if (ispacked(x)) then
          mysize = (x%rank*(x%rank+1)) / 2
          call zaxpy (mysize, alpha, x%zap, 1, y%zap, 1)
      else
          mysize = x%rank * (x%rank)
          call zaxpy (mysize, alpha, x%za, 1, y%za, 1)
      end if
  end subroutine HermitianMatrixAXPY
  !
  !
  subroutine HermitianMatrixcopy (x, y) 
      implicit none
      complex(dp) :: alpha
      type (HermitianMatrix) :: x, y
      integer(i32) :: mysize
      if (ispacked(x)) then
          mysize = (x%rank*(x%rank+1)) / 2
          call zcopy (mysize, x%zap, 1, y%zap, 1)
      else
          mysize = x%rank * (x%rank)
          call zcopy (mysize, x%za, 1, y%za, 1)
      end if
  end subroutine HermitianMatrixcopy
  !
  !
  subroutine HermitianMatrixToFiles (self, prefix)
      implicit none
      type (HermitianMatrix), intent (in) :: self
      character (str_256), intent (in) :: prefix
      character (str_256) :: filename
      if (ispacked(self)) then
          filename = trim (prefix) // ".packed.real.OUT"
          open (888, file=filename)
          write (888,*) dble (self%zap)
      else
          filename = trim (prefix) // ".real.OUT"
          open (888, file=filename)
          write (888,*) dble (self%za)
      end if
      close (888)
      !
      if (ispacked(self)) then
          filename = trim (prefix) // ".packed.imag.OUT"
          open (888, file=filename)
          write (888,*) aimag (self%zap)
      else
          filename = trim (prefix) // ".imag.OUT"
          open (888, file=filename)
          write (888,*) aimag (self%za)
      end if
      close (888)
  end subroutine HermitianMatrixToFiles
  !
  !
  subroutine HermitianMatrixTruncate (self, threshold)
      implicit none
      type (HermitianMatrix), intent (inout) :: self
      real(dp), intent (in) :: threshold
      integer(i32) :: n, i, j
      n = self%rank
      if (ispacked(self)) then
          do i = 1, n * (n+1) / 2
              if (abs(dble(self%zap(i))) .lt. threshold) self%zap(i) = &
              self%zap(i) - dcmplx (dble(self%zap(i)), 0)
              if (abs(aimag(self%zap(i))) .lt. threshold) self%zap(i) &
              = self%zap(i) - dcmplx (0, aimag(self%zap(i)))
          end do
      else
          do j = 1, n
              do i = 1, n
                  if (abs(dble(self%za(i, j))) .lt. threshold) &
                  self%za(i, j) = self%za(i, j) - dcmplx &
                  (dble(self%za(i, j)), 0)
                  if (abs(aimag(self%za(i, j))) .lt. threshold) &
                  self%za(i, j) = self%za(i, j) - dcmplx (0, &
                  aimag(self%za(i, j)))
              end do
          end do
      end if
  end subroutine
  !
  !
  subroutine HermitianMatrixdiagonal (self, d)
      implicit none
      type (HermitianMatrix), intent (in) :: self
      complex(dp), intent (out) :: d (self%rank)
      integer(i32) :: i
      if (ispacked(self)) then
          do i = 1, self%rank
              d (i) = self%zap((i*(i+1))/2)
          end do
      else
          do i = 1, self%rank
              d (i) = self%za(i, i)
          end do
      end if
  end subroutine
  !
  subroutine solvewithlapack(system,nstfv,evecfv,evalfv)
      use mod_timing
      use modinput
      use mod_Gvector, only : ngrid,ngrtot,igfft
      use mod_gkvector, only : ngk
      use m_zfftifc, only: zfftifc
      use iso_c_binding,          only : c_ptr, c_f_pointer
      use mod_device_offload,     only : device_world
       use device_linalg_common_interface, only : zhegvx_gpu
      use m_memory_device,        only : allocate_device_memory, deallocate_device_memory, &
                                   bytes_double_complex, bytes_int, bytes_double_real, get_device_pointer
      implicit none

      type(evsystem), intent(inout) :: system
      integer(i32), intent(in)      :: nstfv
      real(dp), intent(out)         :: evalfv(nstfv)
      complex(dp), intent(out)      :: evecfv(:, :)
      !local
      integer(i32) :: is, ia, i, m, np, info, nmatp, nmatmax, lwork
      real(dp) :: vl, vu
      real(dp) :: ts0, ts1
       ! allocatable arrays
      integer(i32), allocatable :: iwork (:)
      integer(i32), allocatable :: ifail (:)
      real(dp), pointer, contiguous :: w (:)
      real(dp), allocatable :: rwork (:)
      complex(dp), pointer, contiguous :: work (:)
      complex(dp), allocatable :: zfft (:)
      type(c_ptr) :: work_cptr, rwork_cptr, iwork_cptr, ifail_cptr, w_cptr
      integer(i32) :: my_device

      call timesec (ts0)
      nmatmax = size( evecfv, dim=1 )
      if (system%hamilton%packed) then
        vl = 0.0_dp
        vu = 0.0_dp
        ! lapack 3.0 call
        !nmatmax
        nmatp=system%hamilton%rank
        allocate (iwork(5*nmatp))
        allocate (ifail(nmatp))
        allocate (w(nmatp))
        allocate (rwork(7*nmatp))
        allocate (work(2*nmatp))
        call zhpgvx (1, 'V', 'I', 'U', nmatp, system%hamilton%zap, &
                    system%overlap%zap, vl, vu, 1, nstfv, &
                    input%groundstate%solver%evaltol, m, w, evecfv, nmatmax, work, &
                    rwork, iwork, ifail, info)
        evalfv (1:nstfv) = w (1:nstfv)

        if (info /= 0) then
            write (*,*)
            write (*, '("Error(seceqnfv): diagonalisation failed")')
            write (*, '(" ZHPGVX returned INFO = ", I8)') info
            if (info > nmatp) then
                i = info - nmatp
                write (*, '(" The leading minor of the overlap matrix of or&der ", I8)') i
                write (*, '("  is not positive definite")')
                write (*, '(" Order of overlap matrix : ", I8)') nmatp
                write (*,*)
            end if
            stop
        end if
        call timesec (ts1)
        timefv = timefv + ts1 - ts0
        call deletesystem (system)
        deallocate (iwork, ifail, w, rwork, work)
      else
        ! lapack
        ! nmatp=system%hamilton%rank
        ! allocate (iwork(5*nmatp))
        ! allocate (ifail(nmatp))
        ! allocate (w(nmatp))
        ! allocate (rwork(7*nmatp))
        ! allocate (v(1))
        ! allocate (work(2*nmatp))

        !! This segment tests linear dependence of basis and plots the most singular component.
        !! It is meant for educational purposes.
        !
        !write(*,*) ngk(1, 1),nmatp
        !call zheev('V','U',nmatp,system%overlap%za, nmatp,w,work,2*nmatp,rwork,info)
        !write(*,*) w(1)
        !write(*,*)
        !allocate(zfft(ngrtot))
        !
        !zfft(:)=0d0
        !do i = 1, ngk(1, 1)
        !  zfft(igfft(i))=system%overlap%za(i,1)
        !end do
        !
        !call zfftifc (3, ngrid, 1, zfft)
        !write(*,*) sum(zfft)
        !write(*,*)
        !do i=1,48
        !  write(*,*) dble(zfft(i)),dimag(zfft(i))
        !enddo
        !stop

        my_device = device_world%get_device()
        nmatp = system%hamilton%rank
        vl = 0.0_dp
        vu = 0.0_dp

        OMP_OFFLOAD target data map(to:   system%hamilton%za, system%overlap%za)
        OMP_OFFLOAD target data map(from: evecfv)

        ! The workspace is allocated only on the device for GPU-aware compilation
        ! Notice that depending on the vendor some of these arrays will remain unused
        call allocate_device_memory(work_cptr,  bytes_double_complex, my_device)
        call allocate_device_memory(rwork_cptr, 7 * nmatp * bytes_double_real, my_device)
        call allocate_device_memory(iwork_cptr, 5 * nmatp * bytes_int, my_device)
        call allocate_device_memory(ifail_cptr, nmatp * bytes_int, my_device)
        ! The eigenvalues are also stored in the device
        call allocate_device_memory(w_cptr, nmatp * bytes_double_real, my_device) 

        call zhegvx_gpu(1, 'V', 'I', 'U', nmatp, get_device_pointer(system%hamilton%za, my_device), nmatp, &
                        get_device_pointer(system%overlap%za, my_device), nmatp, &
                        vl, vu, 1, nstfv, input%groundstate%solver%evaltol, &
                        m, w_cptr, get_device_pointer(evecfv, my_device), &
                        nmatmax, work_cptr, -1, rwork_cptr, iwork_cptr, ifail_cptr, info, device_world)

        call device_world%synchronize()

        ! Retrieve the proper workspace
        ! For GPU compilation that information lives in the device
        ! and we need to retrieve it from the work c-pointer
        call c_f_pointer(work_cptr, work, [1])
        OMP_OFFLOAD target map(from: lwork) has_device_addr(work)
        lwork = work(1)
        OMP_OFFLOAD end target
        nullify(work)

        ! Resize the workspace pointer
        call deallocate_device_memory(work_cptr, my_device)
        call allocate_device_memory(work_cptr, lwork * bytes_double_complex, my_device)

        call zhegvx_gpu(1, 'V', 'I', 'U', nmatp, get_device_pointer(system%hamilton%za, my_device), nmatp, &
                        get_device_pointer(system%overlap%za, my_device), nmatp, &
                        vl, vu, 1, nstfv, input%groundstate%solver%evaltol, &
                        m, w_cptr, get_device_pointer(evecfv, my_device), &
                        nmatmax, work_cptr, lwork, rwork_cptr, iwork_cptr, ifail_cptr, info, device_world)
        call device_world%synchronize()

        OMP_OFFLOAD end target data
        OMP_OFFLOAD end target data

        call c_f_pointer(w_cptr, w, [nmatp])
        OMP_OFFLOAD target map(from: evalfv) has_device_addr(w)
        evalfv(1:nstfv) = w(1:nstfv)
        OMP_OFFLOAD end target

        call deallocate_device_memory(w_cptr,  my_device)
        call deallocate_device_memory(work_cptr,  my_device)
        call deallocate_device_memory(rwork_cptr, my_device)
        call deallocate_device_memory(iwork_cptr, my_device)
        call deallocate_device_memory(ifail_cptr, my_device)

        if (info /= 0) then
            write (*,*)
            write (*, '("Error(seceqnfv): diagonalization failed")')
            write (*, '(" ZHGVX returned INFO = ", I8)') info
            if (info > nmatp) then
                i = info - nmatp
                write (*, '(" The leading minor of the overlap matrix of order ", I8)') i
                write (*, '("  is not positive definite")')
                write (*, '(" Order of overlap matrix : ", I8)') nmatp
                write (*,*)
            end if
            stop
        end if
        call timesec (ts1)
        timefv = timefv + ts1 - ts0
        call deletesystem (system)
      endif
  end subroutine
end module  modfvsystem
