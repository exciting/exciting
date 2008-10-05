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
integer,save::ikc
real(8),save :: cputot
real(8):: cpuaa,cpualo,cpulolo,cpui,cpu00,cpu01
integer::i,is,ia
complex(8) v(1)
real(8):: cpu0,cpu1
real(8)::threshold
!----------------------------------------!
!     Hamiltonian and overlap set up     !
!----------------------------------------!


call timesec(cpu0)
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

call timesec(cpu1)
 timemat= timemat+cpu1-cpu0


!write(60,*)
!write(60,'("Muffin-tin Hamiltonian setup; Timings (CPU seconds) :")')
!write(60,'(" k-point",T40,": ",I8)') ikc
!write(60,'(" APW-APW",T40,": ",F12.2)') cpuaa
!write(60,'(" APW- lo",T40,": ",F12.2)') cpualo
!write(60,'(" lo - lo",T40,": ",F12.2)') cpulolo
!write(60,'(" interstitial",T40,": ",F12.2)') cpui
!write(60,'(" total",T40,": ",F12.2)') cpuaa+cpualo+cpulolo+cpui
!cputot=cputot+cpuaa+cpualo+cpulolo+cpui
!write(60,'(" cumulative total",T40,": ",F12.2)') cputot
!write(60,*)



end subroutine
