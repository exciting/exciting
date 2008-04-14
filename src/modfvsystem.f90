! Module for setting up the eigensystem
! it is designed in a way that all other subroutines
! dealing with setting up and solving the system can acsess the
! data transparently allowing to choose from different datatypes
! more easily
module modfvsystem
  implicit none

  type HermiteanMatrix

     integer:: rank
     logical:: packed,ludecomposed
     integer,pointer::ipiv(:)
     complex(8), pointer:: za(:,:),zap(:)
  end type HermiteanMatrix

  type evsystem
     type (HermiteanMatrix) ::hamilton, overlap
  end type evsystem

contains
  subroutine newmatrix(self,packed,rank)
    type (HermiteanMatrix),intent(inout)::self
    logical,intent(in)::packed
    integer,intent(in)::rank
    self%rank=rank
    self%packed=packed
    self%ludecomposed=.false.
    if(packed.eqv..true.) then
       allocate(self%zap(rank*(rank+1)/2))
       self%zap=0.0
    else
       allocate(self%za(rank,rank))
       self%za=0.0
    endif
  end subroutine newmatrix

  subroutine deletematrix(self)
    type (HermiteanMatrix),intent(inout)::self
    if(self%packed.eqv..true.) then
       deallocate(self%zap)
    else
       deallocate(self%za)
    endif
    if(self%ludecomposed) deallocate(self%ipiv)
  end subroutine deletematrix




  subroutine newsystem(self,packed,rank)
    type (evsystem),intent(out)::self
    logical,intent(in)::packed
    integer,intent(in)::rank
    call newmatrix(self%hamilton,packed,rank)
    call newmatrix(self%overlap,packed,rank)
  end subroutine newsystem

  subroutine deleteystem(self)
    type(evsystem),intent(inout)::self
    call deletematrix(self%hamilton)
    call deletematrix(self%overlap)
  end subroutine deleteystem

  subroutine Hermiteanmatrix_rank2update(self,n,alpha,x,y)
    type (HermiteanMatrix),intent(inout)::self
    integer,intent(in)::n
    complex(8),intent(in)::alpha,x(:),y(:)

    if(self%packed) then
       call ZHPR2 ( 'U', n, alpha, x, 1, y, 1, self%zap )
    else
       call ZHER2 ( 'U', n, alpha, x, 1, y, 1, self%za,self%rank)
    endif
  end subroutine Hermiteanmatrix_rank2update


  subroutine Hermiteanmatrix_indexedupdate(self,i,j,z)
    type (HermiteanMatrix),intent(inout)::self
    integer::i,j
    complex(8)::z
    integer ipx
    if(self%packed.eqv..true.)then
       ipx=((i-1)*i)/2 + j
       self%zap(ipx)=self%zap(ipx)+z
    else
       if(j.le.i)then
          self%za(j,i)=self%za(j,i)+z
       else
          write(*,*)"warning lower part of hamilton updated"
       endif
    endif
    return
  end subroutine Hermiteanmatrix_indexedupdate

  subroutine Hermiteanmatrixvector(self,alpha,vin,beta,vout)
    implicit none
    type (HermiteanMatrix),intent(inout)::self
    complex(8),intent(in)::alpha,beta
    complex(8),intent(inout)::vin(:)
    complex(8),intent(inout)::vout(:)

    if(self%packed.eqv..true.)then
       call zhpmv("U",self%rank,alpha,self%zap,vin, 1,beta,vout, 1)
    else
       call zhemv("U",self%rank,alpha,self%za,self%rank,vin, 1,beta,vout, 1)
    endif
  end subroutine Hermiteanmatrixvector

  function ispacked(self)
    logical::ispacked
    type(HermiteanMatrix)::self
    ispacked=self%packed
  end function ispacked
  function getrank(self)
    integer:: getrank
    type(HermiteanMatrix)::self
    getrank=self%rank
  end function getrank

  subroutine HermiteanmatrixLU(self)
    type(HermiteanMatrix)::self
    integer info
    if(.not.self%ludecomposed) allocate(self%ipiv(self%rank))

    if(.not.self%ludecomposed)then
       if(.not.ispacked(self))then
          call ZGETRF( self%rank,self%rank, self%za, self%rank, self%IPIV, INFO )
       else
          call ZHPTRF('U',self%rank, self%zap, self%IPIV, INFO )
       endif
       if (info.ne.0)then
          write(*,*)"error in iterativearpacksecequn  HermiteanmatrixLU ",info
          stop
       endif
       self%ludecomposed=.true.
    endif
  end subroutine HermiteanmatrixLU

  subroutine Hermiteanmatrixlinsolve(self,b)
    type(HermiteanMatrix)::self
    complex(8),intent(inout)::b(:)
    integer info
    if(self%ludecomposed) then
       if(.not.ispacked(self))then
          call ZGETRS( 'N', self%rank, 1, self%za, self%rank, self%IPIV, &
               b , self%rank, INFO)
       else
          call ZHPTRS( 'U', self%rank, 1, self%zap, self%IPIV, b, self%rank, INFO )
       endif
       if (info.ne.0)then
          write(*,*)"error in iterativearpacksecequn Hermiteanmatrixlinsolve ",info
          stop
       endif
    endif
  end subroutine Hermiteanmatrixlinsolve

  subroutine HermiteanMatrixAXPY(alpha,x,y)
    complex(8)::alpha
    type(HermiteanMatrix)::x,y
    integer:: mysize
    if (ispacked(x)) then 
       mysize=(x%rank*(x%rank+1))/2
       call zaxpy(mysize,alpha,x%zap,1,y%zap,1)
    else
       mysize=x%rank*(x%rank)
       call zaxpy(mysize,alpha,x%za,1,y%za,1)
    endif
  end subroutine HermiteanMatrixAXPY

  subroutine HermiteanMatrixcopy(x,y)
    complex(8)::alpha
    type(HermiteanMatrix)::x,y
    integer:: mysize
    if (ispacked(x)) then 
       mysize=(x%rank*(x%rank+1))/2
       call zcopy(mysize,x%zap,1,y%zap,1)
    else
       mysize=x%rank*(x%rank)
       call zcopy(mysize,x%za,1,y%za,1)
    endif
  end subroutine HermiteanMatrixcopy

  subroutine HermiteanMatrixToFiles(self,prefix)
    implicit none
    type(HermiteanMatrix),intent(in)::self
    character(256),intent(in)::prefix
    character(256)::filename
    if(ispacked(self)) then
       filename=trim(prefix)//".packed.real.OUT"
       open(888,file=filename)
       write(888,*)dble(self%zap)
    else
       filename=trim(prefix)//".real.OUT"
       open(888,file=filename)
       write(888,*)dble(self%za)
    endif
    close (888)

    if(ispacked(self)) then
       filename=trim(prefix)//".packed.imag.OUT"
       open(888,file=filename)
       write(888,*)aimag(self%zap)
    else
       filename=trim(prefix)//".imag.OUT"
       open(888,file=filename)
       write(888,*)aimag(self%za)
    endif
    close (888)
  end subroutine HermiteanMatrixToFiles

  subroutine HermiteanMatrixTruncate(self,threshold)  
    implicit none
    type(HermiteanMatrix),intent(inout)::self
    real(8),intent(in)::threshold
    integer ::n,i,j
    n= self%rank
    if(ispacked(self)) then 
       do i=1,n*(n+1)/2
          if(abs(dble(self%zap(i))).lt.threshold) self%zap(i)=self%zap(i)-dcmplx(dble(self%zap(i)),0)
          if(abs(aimag(self%zap(i))).lt.threshold) self%zap(i)=self%zap(i)-dcmplx(0,aimag(self%zap(i)))
       end do
    else 
       do j=1,n
          do i=1,n
             if(abs(dble(self%za(i,j))).lt.threshold) self%za(i,j)=self%za(i,j)-dcmplx(dble(self%za(i,j)),0)
             if(abs(aimag(self%za(i,j))).lt.threshold) self%za(i,j)=self%za(i,j)-dcmplx(0,aimag(self%za(i,j)))
          end do
       end do
	endif
  end subroutine  
end module modfvsystem
