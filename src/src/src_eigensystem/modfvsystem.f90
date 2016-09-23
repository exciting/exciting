
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Module for setting up the eigensystem
! it is designed in a way that all other subroutines
! dealing with setting up and solving the system can acsess the
! data transparently allowing to choose from different datatypes
! more easily
Module modfvsystem
   use mod_constants
    Implicit None
    !
    Type HermitianMatrix
        !
        Integer :: rank
        Logical :: packed, ludecomposed, sp
        Integer, Pointer :: ipiv (:)
        Complex(8), Pointer :: za (:, :), zap (:)
        complex(8), pointer :: eigenvalues(:)
        Complex(4), Pointer :: ca (:, :), cap (:)
        complex(4), pointer :: eigenvaluessp(:)

    End Type HermitianMatrix
    !
    Type evsystem
        Type (HermitianMatrix) :: hamilton, overlap
    End Type evsystem
    !
!    Type MatrixArray
!        Type (HermitianMatrix),pointer :: matrix(:)
!    end type

    integer, parameter :: LUdecomp = 1
    integer, parameter :: LDLdecomp = 2
    integer, parameter :: LLdecomp = 3 ! Cholesky decomposition
    integer, parameter :: Diagdecomp = 4 ! decomposition through diagonalization
    integer, parameter :: InverseDiag = 5 ! inverse is approximated by a diagonal matrix 
    integer, parameter :: InvertOnce = 6 ! matrix inversion is done just once, result is reused 
!
! (H-sigma*S)^-1 needed in ARPACK if we choose to store it
    Type (HermitianMatrix), pointer :: arpackinverse(:)
!
Contains
    !
    !
    Subroutine newmatrix (self, packed, rank)
        implicit none
        type (HermitianMatrix), Intent (Inout) :: self
        Logical, Intent (In) :: packed
        Integer, Intent (In) :: rank
        self%rank = rank
        self%packed = packed
        self%ludecomposed = .False.
        nullify(self%ca)
        nullify(self%cap)
        nullify(self%za)
        nullify(self%zap)
        nullify(self%ipiv)
        if (self%sp) then
          If (packed) Then
            Allocate (self%cap(rank*(rank+1)/2))
            self%zap = 0.0
          Else
            Allocate (self%ca(rank, rank))
            self%ca = 0.0
          endif
        else
          If (packed) Then
            Allocate (self%zap(rank*(rank+1)/2))
            self%zap = 0.0
          Else
            Allocate (self%za(rank, rank))
            self%za = 0.0
          End If
        endif
    End Subroutine newmatrix
    !
    !
    Subroutine deletematrix (self)
        implicit none
        Type (HermitianMatrix), Intent (Inout) :: self
        if (associated(self%zap)) Deallocate (self%zap)
        if (associated(self%za))  Deallocate (self%za)
        if (associated(self%cap)) Deallocate (self%cap)
        if (associated(self%ca))  Deallocate (self%ca)  
        If (associated(self%ipiv)) Deallocate (self%ipiv)
    End Subroutine deletematrix
    !
    !
    Subroutine newsystem (self, packed, rank)
        implicit none
        Type (evsystem), Intent (Out) :: self
        Logical, Intent (In) :: packed
        Integer, Intent (In) :: rank
        self%hamilton%sp=.false.
        self%overlap%sp=.false.
        Call newmatrix (self%hamilton, packed, rank)
        Call newmatrix (self%overlap, packed, rank)
    End Subroutine newsystem
    !
    !
    Subroutine deletesystem (self)
        implicit none
        Type (evsystem), Intent (Inout) :: self
        Call deletematrix (self%hamilton)
        Call deletematrix (self%overlap)
    End Subroutine deletesystem
    !
    !
    Subroutine Hermitianmatrix_rank2update (self, n, alpha, x, y) 
        implicit none
        Type (HermitianMatrix), Intent (Inout) :: self
        Integer, Intent (In) :: n
        Complex (8), Intent (In) :: alpha, x (:), y (:)
        !
        If (self%packed) Then
            Call ZHPR2 ('U', n, alpha, x, 1, y, 1, self%zap)
        Else
            Call ZHER2 ('U', n, alpha, x, 1, y, 1, self%za, self%rank)
        End If
    End Subroutine Hermitianmatrix_rank2update
    !
    !
    Subroutine Hermitianmatrix_indexedupdate (self, i, j, z)
        implicit none
        Type (HermitianMatrix), Intent (Inout) :: self
        Integer :: i, j
        Complex (8) :: z
        Integer :: ipx
        If (self%packed .Eqv. .True.) Then
            ipx = ((i-1)*i) / 2 + j
            self%zap (ipx) = self%zap(ipx) + z
        Else
            If (j .Le. i) Then
                self%za (j, i) = self%za(j, i) + z
            Else
                Write (*,*) "warning lower part of hamilton updated"
            End If
        End If
        Return
    End Subroutine Hermitianmatrix_indexedupdate
    !
    !
    Subroutine Hermitianmatrixvector (self, alpha, vin, beta, vout)
#ifdef USEOMP
        use omp_lib
#endif
        Implicit None
        Type (HermitianMatrix), Intent (Inout) :: self
        Complex (8), Intent (In) :: alpha, beta
        Complex (8), Intent (Inout) :: vin (:)
        Complex (8), Intent (Inout) :: vout (:)
        integer nthreads,whichthread,nst,nfin,bandsize,i
        Complex (8), allocatable :: outcome(:)
        !
        If (self%packed .Eqv. .True.) Then
            Call zhpmv ("U", self%rank, alpha, self%zap, vin, 1, beta, &
            vout, 1)
        Else
!         call zgemv ("N", self%rank,self%rank,alpha,self%za,self%rank,vin,1,beta,vout,1)
!         call zgemv ("N", self%rank, 5,alpha,self%za(1:self%rank,1:5),self%rank,vin(1:5),1,beta,vout(1:5),1)

          
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
          endif

          call zgemv ("N", self%rank, nfin-nst+1,alpha,self%za(1,nst),self%rank,vin(nst),1,zzero,outcome,1)
          do i=0,nthreads-1
           if (i.eq.whichthread) vout=vout+outcome
!$OMP BARRIER
          enddo
          deallocate(outcome)
!$OMP END PARALLEL 
#else
           Call zhemv ("U", self%rank, alpha, self%za, self%rank, vin, 1, beta, vout, 1)
#endif
        End If
    End Subroutine Hermitianmatrixvector
    !
    !
    Function ispacked (self)
        implicit none
        Logical :: ispacked
        Type (HermitianMatrix) :: self
        ispacked = self%packed
    End Function ispacked
    !
    !
    Function getrank (self)
        implicit none
        Integer :: getrank
        Type (HermitianMatrix) :: self
        getrank = self%rank
    End Function getrank
    !
    !
    subroutine HermitianMatrixMatrix(self,zm1,zm2,ldzm,naa,ngp)
      Implicit None
      Type (HermitianMatrix),intent(inout) :: self
      Complex(8),intent(in)::zm1(:,:),zm2(:,:)
      integer,intent(in)::ldzm,ngp,naa

      complex(8)::zone=(1.0,0.0)

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
    Subroutine HermitianmatrixInvert (self,Method)
        implicit none
        Type (HermitianMatrix) :: self
        Integer, intent(inout) :: Method
        integer :: lwork,info
        complex(8), allocatable :: work(:)
        complex(8) :: worksize
        integer i,j
       
        if ( ispacked(self)) Then
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
        if (info.ne.0) then
          write(*,*) 'INFO (HermitianmatrixInvert) =', info
          stop
        endif
        if ((Method.eq.LLDecomp).or.(Method.eq.LDLDecomp)) then
! fill the lower triangle too
          do i=2,self%rank
            do j=1,i-1
              self%za(i,j)=conjg(self%za(j,i))
            enddo
          enddo
        endif
    End Subroutine HermitianmatrixInvert
    !
    !
    Subroutine HermitianmatrixFactorize (self,Method)
        implicit none
        Type (HermitianMatrix) :: self
        Integer, intent(inout) :: Method
        if ( ispacked(self)) Then
          write(*,*) 'Packed matrices are not supported'
          stop
        endif
        if ( .Not. self%ludecomposed) then
          select case (Method)
          case (LUdecomp)
            call HermitianmatrixLU (self)
          case (LDLdecomp)
            call HermitianmatrixLDL (self)
          case (LLdecomp)
            call HermitianmatrixLL (self)
            if (associated(self%ipiv)) Method=LDLdecomp
          end select
        self%ludecomposed = .True.
        endif
    End Subroutine HermitianmatrixFactorize
    !
    !
    Subroutine HermitianmatrixLU (self)
        implicit none
        Type (HermitianMatrix) :: self
        Integer :: info
        allocate (self%ipiv(self%rank))
        if (self%sp) then
          Call CGETRF (self%rank, self%rank, self%ca, self%rank, self%ipiv, info)
        else
          Call ZGETRF (self%rank, self%rank, self%za, self%rank, self%ipiv, info)
        endif
        If (info .Ne. 0) Then
          Write (*,*) "error in iterativearpacksecequn  HermitianmatrixLU "                , info
          Stop
        End If
    End Subroutine HermitianmatrixLU
    !
    !
    Subroutine HermitianmatrixLDL (self)
        implicit none
        Type (HermitianMatrix) :: self
        Integer :: info
        complex(8), allocatable :: zwork(:)
        complex(4), allocatable :: cwork(:)
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
        If (info .Ne. 0) Then
          Write (*,*) "error in iterativearpacksecequn  HermitianmatrixLDL "                , info
          Stop
        End If
    End Subroutine HermitianmatrixLDL 
    !
    !
    Subroutine HermitianmatrixLL (self) !Cholesky decomposition
        implicit none
        Type (HermitianMatrix) :: self
        Integer :: info
        complex(8), allocatable :: zwork(:,:)
        complex(4), allocatable :: cwork(:,:)
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
        If (info .Ne. 0) Then
          if (self%sp) self%ca=cwork
          if (.not.self%sp) self%za=zwork
          call HermitianmatrixLDL (self)
        endif
        if (self%sp) then
          deallocate(cwork)
        else
          deallocate(zwork)
        endif
    End Subroutine HermitianmatrixLL
    !
    !
    Subroutine Hermitianmatrixlinsolve (self, b, method)
        implicit none
        Type (HermitianMatrix) :: self
        Complex (8), Intent (Inout) :: b (:)
        integer, Intent (In) :: method
        Integer :: info
        Complex (8), allocatable :: outcome(:)
        complex(8)::zone=(1.0,0.0)
        if ( ispacked(self)) Then
          write(*,*) 'Packed matrices are not supported'
          stop
        endif        
        If (self%ludecomposed) Then
          Select case (method)
          case (LUdecomp)
            Call ZGETRS ('N', self%rank, 1, self%za, self%rank, self%ipiv, b, self%rank, info)
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
            If (info .Ne. 0) Then
                Write (*,*) "error in iterativearpacksecequn Hermitianmatrixlinsolve "                , info
                Stop
            End If
        End If
    End Subroutine Hermitianmatrixlinsolve
    !
    !
    Subroutine HermitianMatrixAXPY (alpha, x, y)
        implicit none
        Complex (8) :: alpha
        Type (HermitianMatrix) :: x, y
        Integer :: mysize
        If (ispacked(x)) Then
            mysize = (x%rank*(x%rank+1)) / 2
            Call zaxpy (mysize, alpha, x%zap, 1, y%zap, 1)
        Else
            mysize = x%rank * (x%rank)
            Call zaxpy (mysize, alpha, x%za, 1, y%za, 1)
        End If
    End Subroutine HermitianMatrixAXPY
    !
    !
    Subroutine HermitianMatrixcopy (x, y) 
        implicit none
        Complex (8) :: alpha
        Type (HermitianMatrix) :: x, y
        Integer :: mysize
        If (ispacked(x)) Then
            mysize = (x%rank*(x%rank+1)) / 2
            Call zcopy (mysize, x%zap, 1, y%zap, 1)
        Else
            mysize = x%rank * (x%rank)
            Call zcopy (mysize, x%za, 1, y%za, 1)
        End If
    End Subroutine HermitianMatrixcopy
    !
    !
    Subroutine HermitianMatrixToFiles (self, prefix)
        Implicit None
        Type (HermitianMatrix), Intent (In) :: self
        Character (256), Intent (In) :: prefix
        Character (256) :: filename
        If (ispacked(self)) Then
            filename = trim (prefix) // ".packed.real.OUT"
            Open (888, File=filename)
            Write (888,*) dble (self%zap)
        Else
            filename = trim (prefix) // ".real.OUT"
            Open (888, File=filename)
            Write (888,*) dble (self%za)
        End If
        Close (888)
        !
        If (ispacked(self)) Then
            filename = trim (prefix) // ".packed.imag.OUT"
            Open (888, File=filename)
            Write (888,*) aimag (self%zap)
        Else
            filename = trim (prefix) // ".imag.OUT"
            Open (888, File=filename)
            Write (888,*) aimag (self%za)
        End If
        Close (888)
    End Subroutine HermitianMatrixToFiles
    !
    !
    Subroutine HermitianMatrixTruncate (self, threshold)
        Implicit None
        Type (HermitianMatrix), Intent (Inout) :: self
        Real (8), Intent (In) :: threshold
        Integer :: n, i, j
        n = self%rank
        If (ispacked(self)) Then
            Do i = 1, n * (n+1) / 2
                If (Abs(dble(self%zap(i))) .Lt. threshold) self%zap(i) = &
                self%zap(i) - dcmplx (dble(self%zap(i)), 0)
                If (Abs(aimag(self%zap(i))) .Lt. threshold) self%zap(i) &
                = self%zap(i) - dcmplx (0, aimag(self%zap(i)))
            End Do
        Else
            Do j = 1, n
                Do i = 1, n
                    If (Abs(dble(self%za(i, j))) .Lt. threshold) &
                    self%za(i, j) = self%za(i, j) - dcmplx &
                    (dble(self%za(i, j)), 0)
                    If (Abs(aimag(self%za(i, j))) .Lt. threshold) &
                    self%za(i, j) = self%za(i, j) - dcmplx (0, &
                    aimag(self%za(i, j)))
                End Do
            End Do
        End If
    End Subroutine
    !
    !
    Subroutine HermitianMatrixdiagonal (self, d)
        Implicit None
        Type (HermitianMatrix), Intent (In) :: self
        Complex (8), Intent (Out) :: d (self%rank)
        Integer :: i
        If (ispacked(self)) Then
            Do i = 1, self%rank
                d (i) = self%zap((i*(i+1))/2)
            End Do
        Else
            Do i = 1, self%rank
                d (i) = self%za(i, i)
            End Do
        End If
    End Subroutine
    !
    Subroutine solvewithlapack(system,nstfv,evecfv,evalfv)
        use mod_timing
        use modinput
        use mod_eigensystem, only:nmatmax
        Use mod_Gvector, only : ngrid,ngrtot,igfft
        use mod_gkvector, only : ngk
        implicit none
        type(evsystem)::system
        integer::nstfv

        Real (8), Intent (Out) :: evalfv (nstfv)
        Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv)
        !local
        Integer :: is, ia, i, m, np, info ,nmatp
        Real (8) :: vl, vu
        Real (8) :: ts0, ts1
         ! allocatable arrays
        Integer, Allocatable :: iwork (:)
        Integer, Allocatable :: ifail (:)
        Real (8), Allocatable :: w (:)
        Real (8), Allocatable :: rwork (:)
        Complex (8), Allocatable :: v (:)
        Complex (8), Allocatable :: work (:)
        Complex (8), Allocatable :: zfft (:)
        Call timesec (ts0)
        if (system%hamilton%packed) then
        vl = 0.d0
        vu = 0.d0
        ! LAPACK 3.0 call
        !nmatmax
        nmatp=system%hamilton%rank
        Allocate (iwork(5*nmatp))
        Allocate (ifail(nmatp))
        Allocate (w(nmatp))
        Allocate (rwork(7*nmatp))
        Allocate (v(1))
        Allocate (work(2*nmatp))
        Call zhpgvx (1, 'V', 'I', 'U', nmatp, system%hamilton%zap, &
        system%overlap%zap, vl, vu, 1, nstfv, &
        input%groundstate%solver%evaltol, m, w, evecfv, nmatmax, work, &
        rwork, iwork, ifail, info)
        evalfv (1:nstfv) = w (1:nstfv)
        !
        !
        !
        If (info .Ne. 0) Then
            Write (*,*)
            Write (*, '("Error(seceqnfv): diagonalisation failed")')
            Write (*, '(" ZHPGVX returned INFO = ", I8)') info
            If (info .Gt. nmatp) Then
                i = info - nmatp
                Write (*, '(" The leading minor of the overlap matrix of or&
           &der ", I8)'                ) i
                Write (*, '("  is not positive definite")')
                Write (*, '(" Order of overlap matrix : ", I8)') nmatp
                Write (*,*)
            End If
            Stop
        End If
        Call timesec (ts1)
        !$OMP CRITICAL
!        timefv = timefv + ts1 - ts0
        !$OMP END CRITICAL
        timefv = timefv + ts1 - ts0
        Call deletesystem (system)
        Deallocate (iwork, ifail, w, rwork, v, work)
        else
        vl = 0.d0
        vu = 0.d0
        ! LAPACK 3.0 call
        !nmatmax
        nmatp=system%hamilton%rank
        Allocate (iwork(5*nmatp))
        Allocate (ifail(nmatp))
        Allocate (w(nmatp))
        Allocate (rwork(7*nmatp))
        Allocate (v(1))
        Allocate (work(2*nmatp))

! This segment tests linear dependence of basis and plots the most singular component.
! It is meant for educational purposes. 
if (.false.) then
        write(*,*) ngk(1, 1),nmatp
        call zheev('V','U',nmatp,system%overlap%za, nmatp,w,work,2*nmatp,rwork,info)
        write(*,*) w(1)
        write(*,*)
        allocate(zfft(ngrtot))

        zfft(:)=0d0
        Do i = 1, ngk(1, 1)
          zfft(igfft(i))=system%overlap%za(i,1)
        End Do

        Call zfftifc (3, ngrid, 1, zfft)        
        write(*,*) sum(zfft)
        write(*,*)
        do i=1,48
          write(*,*) dble(zfft(i)),dimag(zfft(i))
        enddo
        stop 
endif

        !Call zhpgvx (1, 'V', 'I', 'U', nmatp, system%hamilton%zap, &
        !system%overlap%zap, vl, vu, 1, nstfv, &
        !input%groundstate%solver%evaltol, m, w, evecfv, nmatmax, work, &
        !rwork, iwork, ifail, info)
        call ZHEGVX(1, 'V', 'I', 'U', nmatp, system%hamilton%za, nmatp, system%overlap%za, nmatp,&
         vl, vu, 1, nstfv, input%groundstate%solver%evaltol, &
         m, w, evecfv, nmatmax, work, 2*nmatp, rwork, iwork, ifail, info )
        evalfv (1:nstfv) = w (1:nstfv)
        !
        !
        !
        If (info .Ne. 0) Then
            Write (*,*)
            Write (*, '("Error(seceqnfv): diagonalisation failed")')
            Write (*, '(" ZHPGVX returned INFO = ", I8)') info
            If (info .Gt. nmatp) Then
                i = info - nmatp
                Write (*, '(" The leading minor of the overlap matrix of or&
           &der ", I8)'                ) i
                Write (*, '("  is not positive definite")')
                Write (*, '(" Order of overlap matrix : ", I8)') nmatp
                Write (*,*)
            End If
            Stop
        End If
        Call timesec (ts1)
        !$OMP CRITICAL
!        timefv = timefv + ts1 - ts0
        !$OMP END CRITICAL
        timefv = timefv + ts1 - ts0
        Call deletesystem (system)
        Deallocate (iwork, ifail, w, rwork, v, work)
        endif
    end subroutine
End Module modfvsystem
