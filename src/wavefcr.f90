
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine wavefcr(lrstp,is,ia,ist,m,ld,wfcr)
use modmain
implicit none
! arguments
integer, intent(in) :: lrstp
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ist
integer, intent(in) :: m
integer, intent(in) :: ld
complex(8), intent(out) :: wfcr(lmmaxvr,ld,2)
! local variables
integer ias,k,l,lm,ir,irc
real(8) cg1,cg2,t1,t2
l=spl(ist,is)
k=spk(ist,is)
if ((m.lt.-k).or.(m.gt.k-1)) goto 10
ias=idxas(ia,is)
! zero the wavefunction
irc=0
do ir=1,nrmt(is),lrstp
  irc=irc+1
  wfcr(:,irc,:)=0.d0
end do
! calculate the Clebsch-Gordon coefficients
if (k.eq.l+1) then
  cg1=sqrt(dble(l+m+1)/dble(2*l+1))
  cg2=-sqrt(dble(l-m)/dble(2*l+1))
else if (k.eq.l) then
  cg1=sqrt(dble(l-m)/dble(2*l+1))
  cg2=sqrt(dble(l+m+1)/dble(2*l+1))
else
  goto 10
end if
! detemine the two-component spinor
irc=0
do ir=1,nrmt(is),lrstp
  irc=irc+1
! major component of radial wavefunction
  t1=rwfcr(ir,1,ist,ias)/spr(ir,is)
  if (abs(m).le.l) then
    lm=idxlm(l,m)
    t2=t1*cg1
    wfcr(1:lmmaxvr,irc,1)=t2*zbshtvr(1:lmmaxvr,lm)
  end if
  if (abs(m+1).le.l) then
    lm=idxlm(l,m+1)
    t2=t1*cg2
    wfcr(1:lmmaxvr,irc,2)=t2*zbshtvr(1:lmmaxvr,lm)
  end if
end do
return
10 continue
write(*,*)
write(*,'("Error(wavefcr): mismatched l, k or m : ",3I4)') l,k,m
write(*,'(" for species ",I4)') is
write(*,'(" atom ",I4)') ia
write(*,'(" and state ",I6)') ist
write(*,*)
stop
end subroutine

