module modfvsystem
  implicit none
type HermiteanMatrix

integer:: rank
logical:: packed
complex(8), pointer:: za(:,:),zap(:)
end type

type evsystem

type (HermiteanMatrix) ::hamilton, overlap
end type 

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
end subroutine
subroutine deletematrix(self)
type (HermiteanMatrix),intent(inout)::self
if(self%packed.eqv..true.) then
deallocate(self%zap)
else
deallocate(self%za)
endif
end subroutine




subroutine newsystem(self,packed,rank)
type (evsystem),intent(out)::self
logical,intent(in)::packed
integer,intent(in)::rank
call newmatrix(self%hamilton,packed,rank)
call newmatrix(self%overlap,packed,rank)
end subroutine

subroutine deleteystem(self)
type(evsystem),intent(inout)::self
call deletematrix(self%hamilton)
call deletematrix(self%overlap)
end subroutine

subroutine Hermiteanmatrix_rank2update(self,n,alpha,x,y)
type (HermiteanMatrix),intent(inout)::self
integer,intent(in)::n
complex(8),intent(in)::alpha,x(:),y(:)

if(self%packed) then
 call ZHPR2 ( 'U', n, alpha, x, 1, y, 1, self%zap )
    else
 call ZHER2 ( 'U', n, alpha, x, 1, y, 1, self%za,self%rank)
endif
end subroutine


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
end subroutine
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
end subroutine
end module modfvsystem
