subroutine hamiltonandoverlapsetup(np,ik,ispn,apwalm,h,o)
use modmain
complex(8),intent(inout)::h(np),o(np)
integer, intent(in)::np,ik,ispn
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
!local variables
integer::i,is,ia
complex(8) v(1)
!----------------------------------------!
!     Hamiltonian and overlap set up     !
!----------------------------------------!

call cpu_time(cpu0)
! set the matrices to zero
do i=1,np
 h(i)=0.0
 o(i)=0.0
end do
! muffin-tin contributions
do is=1,nspecies
  do ia=1,natoms(is)
    call hmlaa(.false.,is,ia,ngk(ik,ispn),apwalm,v,h)
    call hmlalo(.false.,is,ia,ngk(ik,ispn),apwalm,v,h)
    call hmllolo(.false.,is,ia,ngk(ik,ispn),v,h)
    call olpaa(.false.,is,ia,ngk(ik,ispn),apwalm,v,o)
    call olpalo(.false.,is,ia,ngk(ik,ispn),apwalm,v,o)
    call olplolo(.false.,is,ia,ngk(ik,ispn),v,o)
  end do
end do
! interstitial contributions
call hmlistl(.false.,ngk(ik,ispn),igkig(1,ik,ispn),vgkc(1,1,ik,ispn),v,h)
call olpistl(.false.,ngk(ik,ispn),igkig(1,ik,ispn),v,o)
call cpu_time(cpu1)
end subroutine