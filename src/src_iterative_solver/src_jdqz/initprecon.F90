subroutine initprecon(n,sigma)
use modmain,only:zzero,zone
use modfvsystem
use jacobidavidsoncommon,only: system,p,pdiag
implicit none
integer,intent(in)::n
complex(8),intent(in)::sigma

call newmatrix(p,ispacked(system%hamilton),getrank(system%hamilton))
allocate(pdiag(getrank(system%hamilton)))
 
  call HermiteanMatrixcopy(system%hamilton,p)
  call HermiteanMatrixAXPY(-sigma,system%overlap,p)
  call  HermiteanMatrixdiagonal(p,pdiag)
 ! call HermiteanmatrixLU(p)
 
end subroutine

subroutine deallocprecon()
use jacobidavidsoncommon
use modfvsystem
implicit none
call deletematrix(p)
deallocate(pdiag)
end subroutine
