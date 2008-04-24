subroutine hamiltonandoverlapsetup(np,ngp,apwalm,igpig,vgpc,h,o)
use modmain
implicit none
integer, intent(in)::np,ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8),intent(inout)::h(np),o(np)

!local variables
integer::i,is,ia
complex(8) v(1)
real:: cpu0,cpu1,cpuaa,cpualo,cpulolo,cpui,cpu00,cpu01
real,save :: cputot
data cputot /0.d0/
integer,save::ikc
data ikc /0/

ikc=ikc+1

!----------------------------------------!
!     Hamiltonian and overlap set up     !
!----------------------------------------!

call cpu_time(cpu0)
! set the matrices to zero
do i=1,np
 h(i)=0.0
 o(i)=0.0
end do

cpuaa=0.d0
cpualo=0.d0
cpulolo=0.d0

! muffin-tin contributions
do is=1,nspecies
  do ia=1,natoms(is)

    call cpu_time(cpu00)
    call hmlaa(.false.,is,ia,ngp,apwalm,v,h)
    call cpu_time(cpu01)
    cpuaa=cpuaa+cpu01-cpu00

    call hmlalo(.false.,is,ia,ngp,apwalm,v,h)
    call cpu_time(cpu00)
    cpualo=cpualo+cpu00-cpu01

    call hmllolo(.false.,is,ia,ngp,v,h)
    call cpu_time(cpu01)
    cpulolo=cpulolo+cpu01-cpu00

    call olpaa(.false.,is,ia,ngp,apwalm,v,o)
    call olpalo(.false.,is,ia,ngp,apwalm,v,o)
    call olplolo(.false.,is,ia,ngp,v,o)

  end do
end do

! interstitial contributions
call cpu_time(cpu00)
call hmlistl(.false.,ngp,igpig,vgpc,v,h)
call cpu_time(cpu01)
cpui=cpu01-cpu00

call olpistl(.false.,ngp,igpig,v,o)
call cpu_time(cpu1)
 timemat= timemat+cpu1-cpu0


write(60,*)
write(60,'("Muffin-tin Hamiltonian setup; Timings (CPU seconds) :")')
write(60,'(" k-point",T40,": ",I8)') ikc
write(60,'(" APW-APW",T40,": ",F12.2)') cpuaa
write(60,'(" APW- lo",T40,": ",F12.2)') cpualo
write(60,'(" lo - lo",T40,": ",F12.2)') cpulolo
write(60,'(" interstitial",T40,": ",F12.2)') cpui
write(60,'(" total",T40,": ",F12.2)') cpuaa+cpualo+cpulolo+cpui
cputot=cputot+cpuaa+cpualo+cpulolo+cpui
write(60,'(" cumulative total",T40,": ",F12.2)') cputot
write(60,*)



end subroutine
