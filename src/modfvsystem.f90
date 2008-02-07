module modfvsystem
implicit none
complex(8),allocatable::hamiltonp(:),overlapp(:),hamilton(:,:),overlap(:,:)
logical::packed
integer ::ohrank
interface
subroutine hmlaan(is,ia,ngp,apwalm)
use modmain
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)

end subroutine
end interface

interface
subroutine hmlalon(is,ia,ngp,apwalm)
use modmain

integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)

end subroutine
end interface

interface
subroutine hmlistln(ngp,igpig,vgpc)
use modmain

integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)

end subroutine
end interface
interface
subroutine hmllolon(is,ia,ngp)
use modmain
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp


end subroutine
end interface

interface
subroutine olpaan(is,ia,ngp,apwalm)
use modmain
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
end subroutine
end interface
interface
subroutine olpalon(is,ia,ngp,apwalm)
use modmain
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
end subroutine
end interface

interface
subroutine olpistln(ngp,igpig)
use modmain
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
end subroutine
end interface
interface
subroutine olplolon(is,ia,ngp)
use modmain
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
end subroutine
end interface
interface
subroutine hamiltonandoverlapsetup(ngp,apwalm,igpig,vgpc)
use modmain
implicit none
integer, intent(in)::ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
end subroutine
end interface
end module