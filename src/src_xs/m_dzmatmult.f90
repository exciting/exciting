module m_dzmatmult

  implicit none

  contains

    !BOP
    ! !ROUTINE: dzmatmult
    ! !INTERFACE:
    subroutine dzmatmult(zma, zmb, zmc,&
       & m, n, k, ia, ja, ib, jb, ic, jc,&
       & alpha, beta, transa, transb)
    ! !USES:
      use modmpi 
      use modscl
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    !   type(dzmat) :: zma    ! A distributed complex matrix
    !   type(dzmat) :: zmb    ! B distributed complex matrix
    !   integer(4), optional :: m, n, k ! Dimensions for matrix matrix multiplication
    !   integer(4), optional :: ix, jx  ! Matrix subselection for global arrays (x=a,b,c)
    !                                   ! Indicates upper left corner of sub-matrix in 
    !                                   ! gobal matrix indices.
    !   complex(8), optional :: alpha, beta      ! Scaling factors for alpha*A*B+beta*C=C
    !   character(1), optional :: transa, transb ! Op(A), Op(B), Op='N','T','C'
    ! IN/OUT:
    !   type(dzmat) :: zmc  ! C distributed complex matrix
    !
    ! !DESCRIPTION:
    !   Wrapper for BLAS's complex matrix-matrix multiply routine {\tt ZGEEM}
    !   and its PBLAS's counterpart {\tt PZGEEM}. The $\alpha$ and $\beta$ parameters
    !   default to $1$ and $0$. The {\tt TRANSA TRANSB} parameter default no ``N''.
    !   The ix and jx defalut to 1 and 1. m, n, k defalut values are set
    !   to the respective global dimensions of OP(M), where M is the respective matrix.
    ! 
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      implicit none

      ! Arguments
      type(dzmat), intent(in) :: zma, zmb
      type(dzmat), intent(inout) :: zmc
      integer(4), intent(in), optional :: m, n, k
      integer(4), intent(in), optional :: ia,ja,ib,jb,ic,jc
      complex(8), intent(in), optional :: alpha, beta
      character(1), intent(in), optional :: transa, transb

      ! Local
      character(*), parameter :: thisname = "dzmatmult"
      complex(8) :: a, b
      character(1) :: ta, tb
      integer(4) :: ropa, copa, ropb, copb, rmatc, cmatc
#ifdef SCAL
      integer(4) :: ra,ca,rb,cb,rc,cc
#endif
      logical :: verbose

      verbose = .false.

      ! Setting defaults

      ! Scaling factors
      if(present(alpha)) then
        a = alpha
      else
        a = (1.0d0,0.0d0)
      end if
      if(present(beta)) then
        b = beta
      else
        b = (0.0d0,0.0d0)
      end if

      ! Op A B
      if(present(transa)) then
        ta = transa
      else
        ta = 'N'
      end if
      if(present(transb)) then
        tb = transb
      else
        tb = 'N'
      end if

#ifdef SCAL
      ! Check contexts
      if(zma%context /= zmb%context .or. zma%context /= zmc%context &
        & .or. zmb%context /= zmc%context) then
        write(*, '("Error(",a,"):&
          & PBALS does not perform inter-context operations.")') trim(thisname)
        call terminate
      end if
      ! Subselections of global arrays
      if(present(ia)) then
        ra = ia
      else
        ra = 1
      end if
      if(present(ja)) then
        ca = ja
      else
        ca = 1
      end if
      if(present(ib)) then
        rb = ib
      else
        rb = 1
      end if
      if(present(jb)) then
        cb = jb
      else
        cb = 1
      end if
      if(present(ic)) then
        rc = ic
      else
        rc = 1
      end if
      if(present(jc)) then
        cc = jc
      else
        cc = 1
      end if
#endif

      ! Default values for m,n,k
      ! C = alpha* Op(A)*Op(B) + beta*C
      if(present(m)) then
        ropa = m
      else
        select case(ta)
          case('N','n')
            ropa = zma%nrows
          case('T','t','C','c')
            ropa = zma%ncols
          case default
            write(*, '("Error(",a,"): TRANSA invalid")') trim(thisname)
            call terminate
        end select
      end if
      if(present(n)) then 
        copb = n
      else
        select case(tb)
          case('N','n')
            copb = zmb%ncols
          case('T','t','C','c')
            copb = zmb%nrows
          case default
            write(*, '("Error(",a,"): TRANSB invalid")') trim(thisname)
            call terminate
        end select
      end if
      if(present(k)) then 
        copa = k
      else
        select case(ta)
          case('N','n')
            copa = zma%ncols
          case('T','t','C','c')
            copa = zma%nrows
          case default
            write(*, '("Error(",a,"): TRANSA invalid")') trim(thisname)
            call terminate
        end select
      end if

      ! Checks of dimensions
      select case(tb)
        case('N','n')
          ropb = zmb%nrows
        case('T','t','C','c')
          ropb = zmb%ncols
      end select
      rmatc = zmc%nrows
      cmatc = zmc%ncols
      if(copa > ropb) then 
        write(*,'("dzmatmult@rank",i2,":(ERROR) copa > robp",2i4)') copa, ropb 
        call terminate
      else if(copa /= ropb) then
        write(*,'("dzmatmult@rank",i2,":(Warning) copa /= robp",2i4)') copa, ropb 
      end if
      if(ropa > rmatc) then 
        write(*,'("dzmatmult@rank",i2,":(ERROR) ropa > rmatc",2i4)') ropa, rmatc 
        call terminate
      else if(ropa /= rmatc) then 
        write(*,'("dzmatmult@rank",i2,":(Warning) ropa /= rmatc",2i4)') ropa, rmatc 
      end if
      if(copb > cmatc) then 
        write(*,'("dzmatmult@rank",i2,":(ERROR) copb > cmatc",2i4)') copb, cmatc 
        call terminate
      else if(copb > cmatc) then 
        write(*,'("dzmatmult@rank",i2,":(Warning) copb /= cmatc",2i4)') copb, cmatc 
      end if

#ifdef SCAL
      if(.not. zma%isdistributed) then
        call zgemm(ta, tb, ropa, copb, copa, a, zma%za, zma%nrows,&
          & zmb%za, zmb%nrows, b, zmc%za, zmc%nrows)
      else
        call pzgemm(ta, tb, ropa, copb, copa, a, zma%za, ra, ca, zma%desc,&
          & zmb%za, rb, cb, zmb%desc, b, zmc%za, rc, cc, zmc%desc)
      end if
#else
      call zgemm(ta, tb, ropa, copb, copa, a, zma%za, zma%nrows,&
        & zmb%za, zmb%nrows, b, zmc%za, zmc%nrows)
#endif

    end subroutine dzmatmult
    !EOC
    
end module m_dzmatmult
