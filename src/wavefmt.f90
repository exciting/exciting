
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: wavefmt
! !INTERFACE:
subroutine wavefmt(lrstp,lmax,is,ia,ngp,apwalm,evecfv,ld,wfmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp  : radial step length (in,integer)
!   lmax   : maximum angular momentum required (in,integer)
!   is     : species number (in,integer)
!   ia     : atom number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evecfv : first-variational eigenvector (in,complex(nmatmax))
!   ld     : leading dimension (in,integer)
!   wfmt   : muffin-tin wavefunction (out,complex(ld,*))
! !DESCRIPTION:
!   Calculates the first-variational wavefunction in the muffin-tin in terms of
!   a spherical harmonic expansion. For atom $\alpha$ and a particular
!   $p$-point, the $r$-dependent $(l,m)$-coefficients of the wavefunction for
!   the $i$th state are given by
!   $$ \Psi^{i{\bf p}\alpha}_{lm}(r)=\sum_{\bf G}\Phi^{i{\bf p}}_{\bf G}
!    \sum_{j=1}^{M^{\alpha}_l}A^{\alpha}_{jlm}({\bf G+p})u^{\alpha}_{jl}(r)
!    +\sum_{j=1}^{N^{\alpha}}\Phi^{i{\bf p}}_{(\alpha,j,m)}v^{\alpha}_j(r)
!    \delta_{l,l_j}, $$
!   where $\Phi^{i{\bf p}}$ is the $i$th eigenvector returned from routine
!   {\tt seceqn}; $A^{\alpha}_{jlm}({\bf G+p})$ is the matching coefficient;
!   $M^{\alpha}_l$ is the order of the APW; $u^{\alpha}_{jl}$ is the APW radial
!   function; $N^{\alpha}$ is the number of local-orbitals; $v^{\alpha}_j$ is
!   the $j$th local-orbital radial function; and $(\alpha,j,m)$ is a compound
!   index for the location of the local-orbital in the eigenvector. See routines
!   {\tt genapwfr}, {\tt genlofr}, {\tt match} and {\tt seceqn}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!   Fixed description, October 2004 (C. Brouder)
!   Removed argument ist, November 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lrstp
integer, intent(in) :: lmax
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: evecfv(nmatmax)
integer, intent(in) :: ld
complex(8), intent(out) :: wfmt(ld,*)
! local variables
integer ias,l,m,lm,i
integer ir,nr,io,ilo
real(8) a,b
complex(8) zt1
! external functions
complex(8) zdotu
external zdotu
if (lmax.gt.lmaxapw) then
  write(*,*)
  write(*,'("Error(wavefmt): lmax > lmaxapw : ",I8)') lmax
  write(*,*)
  stop
end if
ias=idxas(ia,is)
! zero the wavefunction
nr=0
do ir=1,nrmt(is),lrstp
  nr=nr+1
  wfmt(:,nr)=0.d0
end do
! APW functions
do l=0,lmax
  do m=-l,l
    lm=idxlm(l,m)
    do io=1,apword(l,is)
      zt1=zdotu(ngp,evecfv,1,apwalm(1,io,lm,ias),1)
      a=dble(zt1)
      b=aimag(zt1)
      call wavefmt_add(nr,ld,wfmt(lm,1),a,b,lrstp,apwfr(1,1,io,l,ias))
    end do
  end do
end do
! local-orbital functions
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  if (l.le.lmax) then
    do m=-l,l
      lm=idxlm(l,m)
      i=ngp+idxlo(lm,ilo,ias)
      a=dble(evecfv(i))
      b=aimag(evecfv(i))
      call wavefmt_add(nr,ld,wfmt(lm,1),a,b,lrstp,lofr(1,1,ilo,ias))
    end do
  end if
end do
return
end subroutine
!EOC

!BOP
! !ROUTINE: wavefmt_add
! !INTERFACE:
subroutine wavefmt_add(nr,ld,wfmt,a,b,lrstp,fr)
! !INPUT/OUTPUT PARAMETERS:
!   nr     : number of radial mesh points (in,integer)
!   ld     : leading dimension (in,integer)
!   wfmt   : complex muffin-tin wavefunction passed in as a real array
!            (inout,real(2*ld,*))
!   a      : real part of complex constant (in,real)
!   b      : imaginary part of complex constant (in,real)
!   lrstp  : radial step length (in,integer)
!   fr     : real radial function (in,real(lrstp,*))
! !DESCRIPTION:
!   Adds a real function times a complex constant to a complex muffin-tin
!   wavefunction as efficiently as possible. See routine {\tt wavefmt}.
!
! !REVISION HISTORY:
!   Created December 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr
integer, intent(in) :: ld
real(8), intent(inout) :: wfmt(2*ld,*)
real(8), intent(in) :: a
real(8), intent(in) :: b
integer, intent(in) :: lrstp
real(8), intent(in) :: fr(lrstp,*)
! local variables
integer ir
! values smaller than eps are taken to be zero
real(8), parameter :: eps=1.d-14
if (abs(b).lt.eps) then
! zero constant
  if (abs(a).lt.eps) return
! pure real constant
  do ir=1,nr
    wfmt(1,ir)=wfmt(1,ir)+a*fr(1,ir)
  end do
else if (abs(a).lt.eps) then
! pure imaginary constant
  do ir=1,nr
    wfmt(2,ir)=wfmt(2,ir)+b*fr(1,ir)
  end do
else
! general complex constant
  do ir=1,nr
    wfmt(1,ir)=wfmt(1,ir)+a*fr(1,ir)
    wfmt(2,ir)=wfmt(2,ir)+b*fr(1,ir)
  end do
end if
return
end subroutine
!EOC

