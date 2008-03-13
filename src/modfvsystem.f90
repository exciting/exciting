module modfvsystem
  implicit none
  type HermiteanMatrix
     private
     integer:: rank
     logical:: packed
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
    complex(8),intent(in)::alpha,vin(:),beta
    complex(8),intent(inout)::vout(:)
    if(self%packed.eqv..true.)then
       call zhpmv("U",self%rank,alpha,self%zap(1),vin(1), 1,beta,vout(1), 1)
    else
       call zhemv("U",self%rank,alpha,self%zap(1),self%rank,vin(1), 1,beta,vout(1), 1)
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

  function getpackedpointer(self)
    complex(8),pointer::getpackedpointer(:)
    type(HermiteanMatrix)::self
    if(ispacked(self))then
       getpackedpointer=>self%zap
    else
       write(*,*)"error in etpackedpointer"
    endif
  end function getpackedpointer

  function get2dpointer(self)
    complex(8),pointer::get2dpointer(:,:)
    type(HermiteanMatrix)::self
    if(.not.ispacked(self))then
       get2dpointer=>self%za
    else
       write(*,*)"error in get2dpointer"
    endif
  end function get2dpointer

end module modfvsystem
