
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genapwfr
! !INTERFACE:
subroutine genapwfr
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the APW radial functions. This is done by integrating the scalar
!   relativistic Schr\"{o}dinger equation (or its energy deriatives) at the
!   current linearisation energies using the spherical part of the effective
!   potential. The number of radial functions at each $l$-value is given by the
!   variable {\tt apword} (at the muffin-tin boundary, the APW functions have
!   continuous derivatives up to order ${\tt apword}-1$). Within each $l$, these
!   functions are orthonormalised with the Gram-Schmidt method. The radial
!   Hamiltonian is applied to the orthonormalised functions and the results are
!   stored in the global array {\tt apwfr}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,nr,ir
integer nn,l,io1,io2
real(8) t1
! automatic arrays
real(8) vr(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
real(8) p0(nrmtmax,apwordmax),p1(nrmtmax),p1s(apwordmax)
real(8) q0(nrmtmax,apwordmax),q1(nrmtmax,apwordmax)
real(8) hp0(nrmtmax)
do is=1,nspecies
  nr=nrmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    vr(1:nr)=veffmt(1,1:nr,ias)*y00
    do l=0,lmaxapw
      do io1=1,apword(l,is)
! integrate the radial Schrodinger equation
        call rschroddme(apwdm(io1,l,is),l,0,apwe(io1,l,ias),nprad,nr, &
         spr(1,is),vr,nn,p0(1,io1),p1,q0(1,io1),q1(1,io1))
! normalise radial functions
        do ir=1,nr
          fr(ir)=p0(ir,io1)**2
        end do
        call fderiv(-1,nr,spr(1,is),fr,gr,cf)
        t1=1.d0/sqrt(abs(gr(nr)))
        p0(1:nr,io1)=t1*p0(1:nr,io1)
        p1s(io1)=t1*p1(nr)
        q0(1:nr,io1)=t1*q0(1:nr,io1)
        q1(1:nr,io1)=t1*q1(1:nr,io1)
! subtract linear combination of previous vectors
        do io2=1,io1-1
          do ir=1,nr
            fr(ir)=p0(ir,io1)*p0(ir,io2)
          end do
          call fderiv(-1,nr,spr(1,is),fr,gr,cf)
          t1=gr(nr)
          p0(1:nr,io1)=p0(1:nr,io1)-t1*p0(1:nr,io2)
          p1s(io1)=p1s(io1)-t1*p1s(io2)
          q0(1:nr,io1)=q0(1:nr,io1)-t1*q0(1:nr,io2)
          q1(1:nr,io1)=q1(1:nr,io1)-t1*q1(1:nr,io2)
        end do
! normalise radial functions
        do ir=1,nr
          fr(ir)=p0(ir,io1)**2
        end do
        call fderiv(-1,nr,spr(1,is),fr,gr,cf)
        t1=abs(gr(nr))
        if (t1.lt.1.d-20) then
          write(*,*)
          write(*,'("Error(genapwfr): degenerate APW radial functions")')
          write(*,'(" for species ",I4)') is
          write(*,'(" atom ",I4)') ia
          write(*,'(" angular momentum ",I4)') l
          write(*,'(" and order ",I4)') io1
          write(*,*)
          stop
        end if
        t1=1.d0/sqrt(t1)
        p0(1:nr,io1)=t1*p0(1:nr,io1)
        p1s(io1)=t1*p1s(io1)
        q0(1:nr,io1)=t1*q0(1:nr,io1)
        q1(1:nr,io1)=t1*q1(1:nr,io1)
! apply the Hamiltonian
        call rschrodapp(l,nr,spr(1,is),vr,p0(1,io1),q0(1,io1),q1(1,io1),hp0)
! divide by r and store in global array
        do ir=1,nr
          t1=1.d0/spr(ir,is)
          apwfr(ir,1,io1,l,ias)=t1*p0(ir,io1)
          apwfr(ir,2,io1,l,ias)=t1*hp0(ir)
        end do
! derivative at the muffin-tin surface
        apwdfr(io1,l,ias)=(p1s(io1)-p0(nr,io1)*t1)*t1
      end do
    end do
  end do
end do
return
end subroutine
!EOC
