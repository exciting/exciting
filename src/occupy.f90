
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: occupy
! !INTERFACE:
subroutine occupy
! !USES:
use modmain
#ifdef TETRA
!<rga>
use modtetra
!</rga>
#endif
! !DESCRIPTION:
!   Finds the Fermi energy and sets the occupation numbers for the
!   second-variational states using the routine {\tt fermi}.
!
! !REVISION HISTORY:
!   Created February 2004 (JKD)
!   Modifiactions for tetrahedron method, November 2007 (RGA 
!     alias Ricardo Gomez-Abal)
!EOP
!BOC
implicit none
! local variables
integer, parameter :: maxit=1000
integer ik,ist,it
real(8) e0,e1,occ,chg,x,t1,fact
! external functions
real(8) sdelta,stheta
external sdelta,stheta
! find minimum and maximum eigenvalues
e0=evalsv(1,1)
e1=e0
do ik=1,nkpt
  do ist=1,nstsv
    e0=min(e0,evalsv(ist,ik))
    e1=max(e1,evalsv(ist,ik))
  end do
end do
if (e0.lt.evalmin) then
  write(*,*)
  write(*,'("Warning(occupy): valence eigenvalues below evalmin for s.c. &
   &loop ",I5)') iscl
end if
if (spinpol) then
  occ=1.d0
else
! double occupancy for spin unpolarised systems
  occ=2.d0
end if
#ifdef TETRA
!<rga>
if (.not.tetra) then
!</rga>
#endif
  t1=1.d0/swidth
  ! determine the Fermi energy using the bisection method
  do it=1,maxit
    efermi=0.5d0*(e0+e1)
    chg=0.d0
    do ik=1,nkpt
      do ist=1,nstsv
        if (evalsv(ist,ik).gt.evalmin) then
          x=(efermi-evalsv(ist,ik))*t1
          occsv(ist,ik)=occ*stheta(stype,x)
          chg=chg+wkpt(ik)*occsv(ist,ik)
        else
          occsv(ist,ik)=0.d0
        end if
      end do
    end do
    if (chg.lt.chgval) then
      e0=efermi
    else
      e1=efermi
    end if
    if ((e1-e0).lt.epsocc) goto 10
  end do
  write(*,*)
  write(*,'("Error(occupy): could not find Fermi energy")')
  write(*,*)
  stop
  10 continue
! find the density of states at the Fermi surface in units of
! states/Hartree/spin/unit cell
  fermidos=0.d0
  do ik=1,nkpt
    do ist=1,nstsv
      if (evalsv(ist,ik).gt.evalmin) then
        x=(efermi-evalsv(ist,ik))*t1
        fermidos=fermidos+wkpt(ik)*sdelta(stype,x)*t1
      end if
    end do
  end do
  fermidos=fermidos*occ/2.d0
#ifdef TETRA
!<rga>
else
  !calculate the fermi energy and the density of states at the fermi energy
  call fermitet(nkpt,nstsv,evalsv,ntet,tnodes,wtet,tvol,chgval,spinpol, &
  efermi,fermidos,.false.)
  call tetiw(nkpt,ntet,nstsv,evalsv,tnodes,wtet,tvol,efermi,occsv)
  do ik=1,nkpt
    ! the "occsv" variable returned from "tetiw" already contains the
    ! weight "wkpt" and does not account for spin degeneracy - rescaling is
    ! necessary (Stephan Sagmeister)
    fact=occ/wkpt(ik)
    do ist=1,nstsv
      occsv(ist,ik)=fact*occsv(ist,ik)
    enddo  
  enddo
  
endif
!</rga>
#endif
return
end subroutine
!EOC
