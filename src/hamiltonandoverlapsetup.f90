subroutine hamiltonandoverlapsetup(ngp,apwalm,igpig,vgpc)
use modfvsystem
use modmain
implicit none
integer, intent(in)::ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
integer ::n

!local variables
integer::i,is,ia
complex(8) v(1)
real:: cpu0,cpu1
!----------------------------------------!
!     Hamiltonian and overlap set up     !
!----------------------------------------!


call cpu_time(cpu0)
! set the matrices to zero

! muffin-tin contributions
do is=1,nspecies
  do ia=1,natoms(is)
    call hmlaan(is,ia,ngp,apwalm)
    call hmlalon(is,ia,ngp,apwalm)
    call hmllolon(is,ia,ngp)
    call olpaan(is,ia,ngp,apwalm)
    call olpalon(is,ia,ngp,apwalm)
    call olplolon(is,ia,ngp)
  end do
end do
	do is=1,ohrank
   		write(888,*)hamiltonp(is),overlapp(is)
 	end do
 	write(*,*)"wrote",ohrank
	stop
! interstitial contributions
call hmlistln(ngp,igpig,vgpc)
call olpistln(ngp,igpig)

if(.not. packed)then
 	call hamiltonoverlapocopy_UL
else
 	do is=1,ohrank
   		write(888,*)hamiltonp(is),overlapp(is)
 	end do
 	write(*,*)"wrote",ohrank
	stop
endif
call cpu_time(cpu1)
 timemat= timemat+cpu1-cpu0
 !do is=1,np
 !write(888,*)h(is),o(is)
 !end do 
 !stop
end subroutine
