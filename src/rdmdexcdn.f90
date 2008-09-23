
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rdmdexcdn(dedn)
! calculates derivative of exchange-correlation energy w.r.t. occupation numbers
use modmain
implicit none
! arguments
real(8), intent(inout) :: dedn(nstsv,nkpt)
! local variables
integer ik1,ik2,ik3
integer ist1,ist2,iv(3)
! parameter for calculating the functional derivatives
real(8), parameter :: eps=1.d-12
real(8) t1,t2,t3,t4
! external functions
real(8) r3taxi
external r3taxi
if (rdmxctype.eq.0) return
! calculate the pre-factor
if (rdmxctype.eq.1) then
  t1=1.d0/occmax
else if (rdmxctype.eq.2) then
  if (spinpol) then
    t1=rdmalpha
  else
    t1=2.d0*rdmalpha*(0.25d0)**rdmalpha
  end if
else
  write(*,*)
  write(*,'("Error(rdmdexcdn): rdmxctype not defined : ",I8)') rdmxctype
  write(*,*)
end if
do ik1=1,nkpt
  do ist1=1,nstsv
    do ik2=1,nkptnr
! find the equivalent reduced k-point
      iv(:)=ivknr(:,ik2)
      ik3=ikmap(iv(1),iv(2),iv(3))
      do ist2=1,nstsv
! Hartree-Fock functional
        if (rdmxctype.eq.1) then
          t2=t1*occsv(ist2,ik3)
! SDLG functional
        else if (rdmxctype.eq.2) then
          if ((ist1.eq.ist2).and. &
           (r3taxi(vkl(1,ik1),vklnr(1,ik2)).lt.epslat)) then
            t2=(1.d0/occmax)*occsv(ist2,ik3)
          else
            t3=max(occsv(ist1,ik1),eps)
            t4=max(occsv(ist2,ik3),eps)
            t2=t1*(t4**rdmalpha)/(t3**(1.d0-rdmalpha))
          end if
        end if
        dedn(ist1,ik1)=dedn(ist1,ik1)+t2*vnlrdm(ist1,ik1,ist2,ik2)
      end do
    end do
  end do
end do
return
end subroutine
