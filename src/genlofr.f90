
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genlofr
! !INTERFACE:
subroutine genlofr
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the local-orbital radial functions. This is done by integrating
!   the scalar relativistic Schr\"{o}dinger equation (or its energy deriatives)
!   at the current linearisation energies using the spherical part of the
!   effective potential. For each local-orbital, a linear combination of
!   {\tt lorbord} radial functions is constructed such that its radial
!   derivatives up to order ${\tt lorbord}-1$ are zero at the muffin-tin radius.
!   This function is normalised and the radial Hamiltonian applied to it. The
!   results are stored in the global array {\tt lofr}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer np,is,ia,ias,nr,ir
integer ilo,io1,io2
integer j,l,nn,info
real(8) t1
! automatic arrays
real(8) vr(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
real(8) p0(nrmtmax,maxlorbord),p1(nrmtmax)
real(8) q0(nrmtmax,maxlorbord),q1(nrmtmax,maxlorbord)
real(8) p0s(nrmtmax),q0s(nrmtmax),q1s(nrmtmax)
real(8) hp0(nrmtmax)
! allocatable arrays
integer, allocatable :: ipiv(:)
real(8), allocatable :: xa(:),ya(:)
real(8), allocatable :: a(:,:),b(:),c(:)
! external functions
real(8) polynom
external polynom
! polynomial order
np=max(maxlorbord+1,4)
allocate(ipiv(np))
allocate(xa(np),ya(np),c(np))
allocate(a(np,np),b(np))
do is=1,nspecies
  nr=nrmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    vr(1:nr)=veffmt(1,1:nr,ias)*y00
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      do io2=1,lorbord(ilo,is)
! integrate the radial Schrodinger equation
        call rschroddme(lorbdm(io2,ilo,is),l,0,lorbe(io2,ilo,ias),nprad,nr, &
         spr(1,is),vr,nn,p0(1,io2),p1,q0(1,io2),q1(1,io2))
! normalise radial functions
        do ir=1,nr
          fr(ir)=p0(ir,io2)**2
        end do
        call fderiv(-1,nr,spr(1,is),fr,gr,cf)
        t1=1.d0/sqrt(abs(gr(nr)))
        p0(1:nr,io2)=t1*p0(1:nr,io2)
        q0(1:nr,io2)=t1*q0(1:nr,io2)
        q1(1:nr,io2)=t1*q1(1:nr,io2)
! set up the matrix of radial derivatives
        do j=1,np
          ir=nr-np+j
          xa(j)=spr(ir,is)
          ya(j)=p0(ir,io2)/spr(ir,is)
        end do
        do io1=1,lorbord(ilo,is)
          a(io1,io2)=polynom(io1-1,np,xa,ya,c,rmt(is))
        end do
      end do
! set up the target vector
      b(:)=0.d0
      b(lorbord(ilo,is))=1.d0
      call dgesv(lorbord(ilo,is),1,a,np,ipiv,b,np,info)
      if (info.ne.0) then
        write(*,*)
        write(*,'("Error(genlofr): degenerate local-orbital radial functions")')
        write(*,'(" for species ",I4)') is
        write(*,'(" atom ",I4)') ia
        write(*,'(" and local-orbital ",I4)') ilo
        write(*,'(" ZGESV returned INFO = ",I8)') info
        write(*,*)
        stop
      end if
! generate linear superposition of radial functions
      p0s(:)=0.d0
      q0s(:)=0.d0
      q1s(:)=0.d0
      do io1=1,lorbord(ilo,is)
        t1=b(io1)
        p0s(1:nr)=p0s(1:nr)+t1*p0(1:nr,io1)
        q0s(1:nr)=q0s(1:nr)+t1*q0(1:nr,io1)
        q1s(1:nr)=q1s(1:nr)+t1*q1(1:nr,io1)
      end do
! normalise radial functions
      do ir=1,nr
        fr(ir)=p0s(ir)**2
      end do
      call fderiv(-1,nr,spr(1,is),fr,gr,cf)
      t1=1.d0/sqrt(abs(gr(nr)))
      p0s(1:nr)=t1*p0s(1:nr)
      q0s(1:nr)=t1*q0s(1:nr)
      q1s(1:nr)=t1*q1s(1:nr)
! apply the scalar relativistic Hamiltonian
      call rschrodapp(l,nr,spr(1,is),vr,p0s,q0s,q1s,hp0)
! divide by r and store in global array
      do ir=1,nr
        t1=1.d0/spr(ir,is)
        lofr(ir,1,ilo,ias)=t1*p0s(ir)
        lofr(ir,2,ilo,ias)=t1*hp0(ir)
      end do
    end do
  end do
end do
deallocate(ipiv,xa,ya,a,b,c)
return
end subroutine
!EOC
