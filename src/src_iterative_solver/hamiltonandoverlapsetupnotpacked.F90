
subroutine hamiltonandoverlapsetupnotpacked(n,ngp,apwalm,igpig,vgpc,hamilton,overlap)
  use modmain, only:ngkmax,apwordmax,lmmaxapw,natmtot
  use diisinterfaces
  implicit none
  integer, intent(in)::n,ngp
  complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  integer, intent(in) :: igpig(ngkmax)
  real(8), intent(in) :: vgpc(3,ngkmax)
  complex(8),intent(inout)::hamilton(n,n),overlap(n,n)
  complex(8)::hp(n*(n+1)/2),op(n*(n+1)/2)

  integer ::packedi, icolumn
  hamilton=0
  overlap=0
  call hamiltonandoverlapsetup((n*n+1)/2,ngp,apwalm,igpig,vgpc,hp,op)
  packedi=1
  Do icolumn=1,n
     call zcopy(icolumn,hp(packedi),1,hamilton(1,icolumn),1)
     call zcopy(icolumn,op(packedi),1,overlap(1,icolumn),1)
     packedi=packedi+icolumn
  end do


 
end subroutine hamiltonandoverlapsetupnotpacked
