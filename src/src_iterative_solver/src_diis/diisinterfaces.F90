


module diisinterfaces

  implicit none
  complex(8) zdotc
  real(8)::dlamch
  external zdotc, dlamch
   real(8), external::dznrm2
#ifdef USEDIISIF
  interface


subroutine  DIISseceqnfv(ik, ispn, apwalm, vgpc, evalfv, evecfv)
   use modmain, only: nstfv, vkl, ngk, igkig, nmat, vgkl, timemat, npmat&
       , apwordmax, lmmaxapw, natmtot, nkpt, nmatmax, nspnfv, timefv, ngkmax
  integer,	intent(in)		:: ik
  integer,	intent(in)		:: ispn
  real(8),    intent(in)    :: vgpc(3, ngkmax)
  complex(8), intent(in)	:: apwalm(ngkmax, apwordmax, lmmaxapw, natmtot)
  real(8),	intent(inout)	:: evalfv(nstfv, nspnfv)
  complex(8), intent(inout) :: evecfv(nmatmax, nstfv, nspnfv)
       end subroutine
  end interface

  interface


subroutine seceqfvprecond  (n, system, X, w, evalfv, evecfv)
       use modmain, only: nmatmax, nstfv
      use  modfvsystem
	     implicit none
  type(evsystem)::system

       integer, intent(in)::n

       complex(8), intent(out)::evecfv(nmatmax, nstfv)
       real(8), intent(out)::evalfv(nstfv), w(nmatmax)
       complex(8), intent (OUT)::X(nmatmax, nmatmax)

     end subroutine seceqfvprecond
  end interface

  interface


subroutine prerotate_preconditioner(n, m, h, P)
       use modmain , only:nstfv, nmatmax
       implicit none
       integer, intent(in) :: n, m
       complex(8), intent(in)::h(n, n)
       complex(8), intent(inout)::P(nmatmax, nmatmax)

     end subroutine prerotate_preconditioner
  end interface

  interface


subroutine precondspectrumupdate(n, m, hamilton, overlap, P, w)
       use modmain , only:nstfv, nmatmax
       implicit none
       integer, intent(in) :: n, m
       complex(8), intent(in)::hamilton(n, n), overlap(n, n)
       complex(8), intent(in)::P(nmatmax, nmatmax)
       real(8), intent(inout)::w(nmatmax)
     end subroutine precondspectrumupdate
  end interface
  interface


subroutine	 readprecond(ik, n, X, w)
       use modmain
       use modmpi
       integer, intent(in)::n, ik
       complex(8), intent(out)::X(nmatmax, nmatmax)
       real(8), intent(out)::w(nmatmax)
     end subroutine readprecond
  end interface
  interface


subroutine  writeprecond(ik, n, X, w)
       use modmain
       use modmpi
       integer, intent(in)::n, ik
       complex(8), intent(in)::X(nmatmax, nmatmax)
       real(8), intent(in)::w(nmatmax)
     end subroutine writeprecond
  end interface
  interface


subroutine setuphsvect(n, m, system, evecfv, ldv, h, s)

  use modmain, only : nstfv, zone, zzero
   use modfvsystem
  implicit none
  integer , intent(in):: n, m, ldv
  type(evsystem)::system
  complex(8), intent(in):: evecfv(ldv, m)
  complex(8), intent(out)::h(n, m), s(n, m)

     end subroutine setuphsvect
  end interface
  interface


subroutine rayleighqotient(n, m, evecfv, h, s, evalfv)

  implicit none
  integer, intent(in)::n, m
  complex(8) , intent(in)::h(n, m), s(n, m), evecfv(n, m)
  real(8) , intent(out)::evalfv(m)
     end subroutine rayleighqotient
  end interface
  interface


subroutine residualvectors(n, iunconverged, h, s, evalfv, r, rnorms)
       use modmain, only: nmatmax, nstfv
       implicit none
       integer , intent (in)::n, iunconverged
       !packed ut
       complex(8), intent(in)::h(n, nstfv), s(n, nstfv)
       complex(8), intent(out)::r(n, nstfv)
       real(8), intent(in)::evalfv(nstfv)
       real(8), intent(out)::rnorms(nstfv)
     end subroutine residualvectors
  end interface
  interface


subroutine calcupdatevectors(n, iunconverged, P, w, r, evalfv, evecfv, phi)
  use modmain, only:nstfv, zzero, zone, nmatmax

  implicit none
  integer , intent (in)::n , iunconverged
  complex(8), intent(in)::P(nmatmax, nmatmax)
  complex(8), intent(in)::r(n, nstfv), evecfv(n, nstfv)
  complex(8), intent(out)::phi(n, nstfv)
  real(8), intent(in)::w(nmatmax), evalfv(nstfv)
  end subroutine calcupdatevectors

  end interface
  interface


subroutine    diisupdate(idiis, icurrent, iunconverged, n, h, s&
     , trialvec, evalfv , evecfv, infodiisupdate)

  use modmain, only: nstfv, zone, zzero
  use sclcontroll, only:diismax
  implicit none
  integer , intent(in)::idiis, icurrent, iunconverged, n
  complex(8), intent(in)::h(n, nstfv, diismax)
  complex(8), intent(in):: s(n, nstfv, diismax), trialvec(n, nstfv, diismax)
  real(8), intent(in):: evalfv(nstfv, diismax)
  complex(8), intent(out)::evecfv(n, nstfv)
  integer, intent(out)::infodiisupdate

     end subroutine diisupdate
  end interface
    interface


subroutine  normalize(n, m, overlap, evecfv, ldv)
use modmain, only:zone, zzero

implicit none
integer, intent (in)::n, m, ldv
complex(8), intent(in)::overlap(n, n)
complex(8), intent(out)::evecfv(ldv, m)
end subroutine
 end interface
  interface


subroutine solvediis(m, Pmatrix, Qmatrix, c)
       implicit none
       integer, intent(in)::m
     real(8), intent(in)::Pmatrix(m, m), Qmatrix(m, m)
       real(8), intent(out)::c(m)
     end subroutine solvediis
  end interface
  interface


subroutine exactupdatevectors(n, iunconverged, hamilton, overlap, r, rhizvalue, eigenvector, trialvecs)
!calculate update equation with linsolver

!solvefor dA:  dA=(H-e*S)\R

! dA 	Update step to zero residual
! H 	Hamilton
! S 	Overlap
! e 	Rhitz Value
! R 	Residual

!trialvecs=eigenvector+dA
use modfvsystem
use modmain, only:zone
integer, intent(in):: n, iunconverged
type(HermiteanMatrix), intent(in)::hamilton, overlap
complex(8), intent(in)::r(n, iunconverged), eigenvector(n, iunconverged)
real(8), intent(in)::rhizvalue(iunconverged)
complex(8), intent(out)::trialvecs(n, iunconverged)

end subroutine
end interface
#endif
end module diisinterfaces
