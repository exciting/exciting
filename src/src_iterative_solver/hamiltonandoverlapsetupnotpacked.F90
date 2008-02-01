
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

integer ::ip, icolumn,irow
hamilton=0
overlap=0
call hamiltonandoverlapsetup((n*(n+1))/2,ngp,apwalm,igpig,vgpc,hp,op)

ip=1
do icolumn=1,n
if(icolumn.gt.1)then
   do irow=1 ,icolumn-1
      hamilton(irow,icolumn)=hp(ip)
      hamilton(icolumn,irow)=conjg(hp(ip))
      ip=ip+1
   end do
  endif
     hamilton(icolumn,icolumn)= dcmplx(dble(hp(ip)))
     ip=ip+1
end do
ip=1
do icolumn=1,n
if(icolumn.gt.1)then
   do irow=1,icolumn-1
      overlap(irow,icolumn)=op(ip)
      overlap(icolumn,irow)=conjg(op(ip))
      ip=ip+1 
   end do
   endif
      overlap(icolumn,icolumn)=dcmplx(dble(op(ip)))
   ip=ip+1
end do


end subroutine 