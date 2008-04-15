subroutine hamiltonandoverlapsetup(system,ngp,apwalm,igpig,vgpc)
use modfvsystem
use modmain
implicit none
type(evsystem)::system
integer, intent(in)::ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
integer ::n
character(256)::prefix
!local variables
integer::i,is,ia
complex(8) v(1)
real:: cpu0,cpu1
real(8)::threshold
!----------------------------------------!
!     Hamiltonian and overlap set up     !
!----------------------------------------!


call cpu_time(cpu0)
! set the matrices to zero

! muffin-tin contributions
do is=1,nspecies
  do ia=1,natoms(is)
    call hmlaan(system%hamilton,is,ia,ngp,apwalm)
    call hmlalon(system%hamilton,is,ia,ngp,apwalm)
    call hmllolon(system%hamilton,is,ia,ngp)
    call olpaan(system%overlap,is,ia,ngp,apwalm)
    call olpalon(system%overlap,is,ia,ngp,apwalm)
    call olplolon(system%overlap,is,ia,ngp)
  end do
end do
! interstitial contributions
call hmlistln(system%hamilton,ngp,igpig,vgpc)
call olpistln(system%overlap,ngp,igpig)
threshold=1e-16
!call HermiteanMatrixTruncate(system%hamilton,threshold)
!call HermiteanMatrixTruncate(system%overlap,threshold)

!

if(.not.ispacked(system%hamilton))then
 	call hamiltonoverlapocopy_UL(system)
endif
#ifdef DEBUGHO
write(*,*)"apwalm", apwalm
prefix="H"
 call HermiteanMatrixToFiles(system%hamilton,prefix)
prefix="O"
 call HermiteanMatrixToFiles(system%overlap,prefix)		
 	write(*,*)"wrote" 
	stop
#endif 

call cpu_time(cpu1)
 timemat= timemat+cpu1-cpu0

end subroutine
