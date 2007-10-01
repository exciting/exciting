
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: mossbauer
! !INTERFACE:
subroutine mossbauer
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the contact charge density and contact magnetic hyperfine field for
!   each atom and outputs the data to the file {\tt MOSSBAUER.OUT}. The nuclear
!   radius used for the contact quantities is approximated by the empirical
!   formula $R_{\rm N}=1.2 Z^{1/3}$ fm, where $Z$ is the atomic number.
!
! !REVISION HISTORY:
!   Created May 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ir,nr
real(8) t1,rn,vn,rho0,b
! allocatable arrays
real(8), allocatable :: fr(:)
real(8), allocatable :: gr(:)
real(8), allocatable :: cf(:,:)
! initialise universal variables
call init0
! read density and potentials from file
call readstate
! allocate local arrays
allocate(fr(nrmtmax))
allocate(gr(nrmtmax))
allocate(cf(3,nrmtmax))
open(50,file='MOSSBAUER.OUT',action='WRITE',form='FORMATTED')
do is=1,nspecies
!--------------------------------!
!     contact charge density     !
!--------------------------------!
! determine approximate nuclear radius : 1.2*A^(1/3) fm
  t1=1.2d-15/0.5291772108d-10
  rn=t1*abs(spzn(is))**(1.d0/3.d0)
  do ir=1,nrmt(is)
    if (spr(ir,is).gt.rn) goto 10
  end do
  write(*,*)
  write(*,'("Error(mossbauer): nuclear radius too large : ",G18.10)') rn
  write(*,'(" for species ",I4)') is
  write(*,*)
  stop
10 continue
  nr=ir
  rn=spr(nr,is)
! nuclear volume
  vn=(4.d0/3.d0)*pi*rn**3
  do ia=1,natoms(is)
    ias=idxas(ia,is)
!--------------------------------!
!     contact charge density     !
!--------------------------------!
    fr(1:nr)=rhomt(1,1:nr,ias)*y00
    do ir=1,nr
      fr(ir)=(fourpi*spr(ir,is)**2)*fr(ir)
    end do
    call fderiv(-1,nr,spr(1,is),fr,gr,cf)
    rho0=gr(nr)/vn
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
    write(50,'(" approximate nuclear radius              : ",G18.10)') rn
    write(50,'(" number of mesh points to nuclear radius : ",I6)') nr
    write(50,'(" contact charge density                  : ",G18.10)') rho0
!------------------------------------------!
!     contact magnetic hyperfine field     !
!------------------------------------------!
    if (spinpol) then
      do ir=1,nr
        if (ndmag.eq.3) then
          t1=sqrt(magmt(1,ir,ias,1)**2+magmt(1,ir,ias,2)**2 &
           +magmt(1,ir,ias,3)**2)
        else
          t1=magmt(1,ir,ias,1)
        end if
        fr(ir)=t1*y00*fourpi*spr(ir,is)**2
      end do
      call fderiv(-1,nr,spr(1,is),fr,gr,cf)
      b=gr(nr)/vn
      write(50,'(" contact magnetic hyperfine field (mu_B) : ",G18.10)') b
    end if
  end do
end do
close(50)
write(*,*)
write(*,'("Info(mossbauer):")')
write(*,'(" Mossbauer parameters written to MOSSBAUER.OUT")')
write(*,*)
deallocate(fr,gr,cf)
return
end subroutine
!EOC
