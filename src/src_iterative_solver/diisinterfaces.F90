module diisinterfaces

  implicit none
  complex(8) zdotc
  real(8) dlamch
  external zdotc,dlamch
  
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
     subroutine seceqfvprecond  (n,h,o,X,w,evalfv,evecfv)
       use modmain, only: nmatmax,nstfv
    
       implicit none
       integer, intent(in)::n
       complex(8),intent(in)::h(n,n),o(n,n)
       complex(8),intent(out)::evecfv(nmatmax,nstfv)
       real(8),intent(out)::evalfv(nstfv), w(n)
       complex(8),intent (OUT)::X(nmatmax,nmatmax)

     end subroutine seceqfvprecond
  end interface

  interface
     subroutine prerotate_preconditioner(n,m,h,evecfv,P)
       use modmain ,only:nstfv,nmatmax
       implicit none
       integer, intent(in) :: n,m
       complex(8), intent(in)::h(n,n),evecfv(nstfv)
       complex(8),intent(inout)::P(nmatmax,nmatmax)

     end subroutine prerotate_preconditioner
  end interface

  interface
     subroutine precondspectrumupdate(n,m,hamilton,overlap,P,w)
       use modmain ,only:nstfv,nmatmax
       implicit none
       integer, intent(in) :: n,m
       complex(8), intent(in)::hamilton(n,n),overlap(n,n)
       complex(8),intent(in)::P(nmatmax,nmatmax)
       real(8),intent(inout)::w(nmatmax) 
     end subroutine precondspectrumupdate
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
     subroutine  writeprecond(ik,n,X,w)
       use modmain
       use modmpi
       integer, intent(in)::n,ik
       complex(8), intent(in)::X(nmatmax,nmatmax)
       real(8),intent(in)::w(nmatmax)
     end subroutine writeprecond
  end interface
  interface
     subroutine setuphsvect(n,m,hamilton,overlap,evecfv,h,s)
       use modmain, only : nmatmax,nstfv
       implicit none
       integer ,intent(in):: n,m
       complex(8), intent(in):: hamilton(n,n),overlap(n,n),evecfv(nmatmax,m)
       complex(8), intent(out)::h(n,m),s(n,m)

     end subroutine setuphsvect
  end interface
  interface
     subroutine rayleighqotient(n,m,evecfv, h,s,evalfv)
       use modmain, only: nstfv
       implicit none
       integer, intent(in)::n,m
       complex(8) ,intent(in)::h(n),s(n),evecfv(n,nstfv)
       real(8) ,intent(out)::evalfv(nstfv)
     end subroutine rayleighqotient
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
     end subroutine residualvectors
  end interface
  interface
     subroutine calcupdatevectors(n,iunconverged,P,w,r,evalfv,phi) 
       use modmain, only:nstfv,nmatmax
       integer ,intent (in)::n , iunconverged
       complex(8),intent(in)::P(nmatmax,nmatmax),r(n,nstfv)
       complex(8),intent(out)::phi(n,nstfv)
       real(8), intent(in)::w(nmatmax),evalfv(nstfv)
     end subroutine calcupdatevectors
  end interface
  interface
     subroutine   diisupdate(idiis,iunconverged,n,h,s,trialvec,evalfv ,evecfv)
       use modmain,only: nstfv,nmatmax
       implicit none
       integer ,intent(in)::idiis,iunconverged,n
       complex(8),intent(in)::h(n,nstfv,idiis),s(n,nstfv,idiis),trialvec(n,nstfv,idiis)
       real(8), intent(in):: evalfv(nstfv)
       complex(8),intent(out)::evecfv(nmatmax,nstfv)
     end subroutine diisupdate
  end interface

  interface
     subroutine solvediis(m,Pmatrix,Qmatrix,c)
       implicit none
       integer, intent(in)::m
       complex(8), intent(in)::Pmatrix(m,m),Qmatrix(m,m)
       complex(8), intent(out)::c(m)
     end subroutine solvediis
  end interface
end module diisinterfaces

