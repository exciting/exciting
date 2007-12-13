module diisinterfaces
  implicit none

  interface
     subroutine  hamiltonandoverlapsetupnotpacked(n,ngp,apwalm,igpig,vgpc,hamilton,overlap)
       use modmain, only:ngkmax,apwordmax,lmmaxapw,natmtot
       integer, intent(in)::n,ngp
       complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
       integer, intent(in) :: igpig(ngkmax)
       real(8), intent(in) :: vgpc(3,ngkmax)
       complex(8),intent(inout)::hamilton(n,n),overlap(n,n)
     end subroutine hamiltonandoverlapsetupnotpacked
  end interface

  interface
     subroutine prerotate_preconditioner(n,h,o,evecfv,X,w)
       use modmain ,only:nstfv,nmatmax
       integer, intent(in) :: n
       complex(8), intent(in)::h(n,n),o(n,n),evecfv(nstfv)
       complex(8),intent(inout)::X(nmatmax,nmatmax)
       real(8),intent(inout)::w(n)
     end subroutine prerotate_preconditioner
  end interface

  interface
     subroutine       readprecond(ik,n,X,w)
       use modmain
       use modmpi	
       integer, intent(in)::n,ik
       complex(8), intent(out)::X(nmatmax,nmatmax)
       real(8),intent(out)::w(nmatmax)
     end subroutine readprecond
  end interface
 interface
    subroutine setuphsvect(n,hamilton,overlap,evecfv,h,s)
      use modmain, only : nmatmax,nstfv
      integer ,intent(in):: n
      complex(8), intent(in):: hamilton(n,n),overlap(n,n),evecfv(nmatmax,nstfv)
      complex(8), intent(out)::h(n,nstfv),s(n,nstfv)
    end subroutine setuphsvect
 end interface
interface
subroutine rayleighqotient(n,m,evecfv, h,s,evalfv)
  use modmain, only: nstfv
  implicit none
  integer, intent(in)::n,m
  complex(8) ,intent(in)::h(n),s(n),evecfv(n,nstfv)
  real(8) ,intent(out)::evalfv(nstfv)
    end subroutine 
 end interface
interface
subroutine residualvectors(n,iunconverged,h,s,evalfv,r,rnorms)
use modmain, only: nmatmax,nstfv
implicit none
integer , intent (in)::n,iunconverged
!packed ut
complex(8),intent(in)::h(n,nstfv),s(n,nstfv) 
complex(8),intent(out)::r(n,nstfv)
real(8),intent(in)::evalfv(nstfv)
real(8),intent(out)::rnorms(nstfv)
   end subroutine 
 end interface

end module diisinterfaces
