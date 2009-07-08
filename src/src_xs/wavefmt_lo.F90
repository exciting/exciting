

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: wavefmt_lo
! !INTERFACE:


subroutine wavefmt_lo(lrstp, lmax, is, ia, ngp, evecfv, ld, wfmt)
! !USES:
use modinput
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lrstp  : radial step length (in,integer)
!   lmax   : maximum angular momentum required (in,integer)
!   is     : species number (in,integer)
!   ia     : atom number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   evecfv : first-variational eigenvector (in,complex(nmatmax))
!   ld     : leading dimension (in,integer)
!   wfmt   : muffin-tin wavefunction (out,complex(ld,*))
! !DESCRIPTION:
!   Muffin-tin wavefunction built up by local orbital contribution only.
!   Based upon the routine {\tt wavefmt}.
!
! !REVISION HISTORY:
!   Created May 2008 (Sagmeister)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lrstp
integer, intent(in) :: lmax
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp
complex(8), intent(in) :: evecfv(nmatmax)
integer, intent(in) :: ld
complex(8), intent(out) :: wfmt(ld, *)
! local variables
integer::ias, l, m, lm, i
integer::ir, nr, ilo
real(8)::a, b
! external functions
complex(8) zdotu
external zdotu
if (lmax.gt.input%groundstate%lmaxapw) then
  write(*, *)
  write(*, '("Error(wavefmt): lmax > lmaxapw : ", I8)') lmax
  write(*, *)
  stop
end if
ias=idxas(ia, is)
! zero the wavefunction
nr=0
do ir=1, nrmt(is), lrstp
  nr=nr+1
  wfmt(:, nr)=0.d0
end do
!!$! APW functions
!!$do l=0,lmax
!!$  do m=-l,l
!!$    lm=idxlm(l,m)
!!$    do io=1,apword(l,is)
!!$      zt1=zdotu(ngp,evecfv,1,apwalm(1,io,lm,ias),1)
!!$      a=dble(zt1)
!!$      b=aimag(zt1)
!!$      call wavefmt_add(nr,ld,wfmt(lm,1),a,b,lrstp,apwfr(1,1,io,l,ias))
!!$    end do
!!$  end do
!!$end do
! local-orbital functions
do ilo=1, nlorb(is)
  l=lorbl(ilo, is)
  if (l.le.lmax) then
    do m=-l, l
      lm=idxlm(l, m)
      i=ngp+idxlo(lm, ilo, ias)
      a=dble(evecfv(i))
      b=aimag(evecfv(i))
      call wavefmt_add(nr, ld, wfmt(lm, 1), a, b, lrstp, lofr(1, 1, ilo, ias))
    end do
  end if
end do
return
end subroutine
!EOC
