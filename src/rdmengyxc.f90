
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine rdmengyxc
! calculate the RDMFT exchange-correlation energy
use modmain
implicit none
! local variables
integer ik1,ik2,ik3
integer ist1,ist2,iv(3)
real(8) t1,t2,t3
! external functions
real(8) r3taxi,rfmtinp
external r3taxi,rfmtinp
! calculate the prefactor
if (rdmxctype.eq.0) then
  engyx=0.d0
  return
else if (rdmxctype.eq.1) then
  t1=0.5d0/occmax
else if (rdmxctype.eq.2) then
  if (spinpol) then
    t1=0.5d0
  else
    t1=(0.25d0)**rdmalpha
  end if
else
  write(*,*)
  write(*,'("Error(rdmengyxc): rdmxctype not defined : ",I8)') rdmxctype
  write(*,*)
  stop
end if
! exchange-correlation energy
engyx=0.d0
do ik1=1,nkpt
  do ist1=1,nstsv
    do ik2=1,nkptnr
! find the equivalent reduced k-point
      iv(:)=ivknr(:,ik2)
      ik3=ikmap(iv(1),iv(2),iv(3))
      do ist2=1,nstsv
! Hartree-Fock functional
        if (rdmxctype.eq.1) then
          t2=t1*wkpt(ik1)*occsv(ist1,ik1)*occsv(ist2,ik3)
! SDLG functional
        else if (rdmxctype.eq.2) then
          t3=occsv(ist1,ik1)*occsv(ist2,ik3)
          if ((ist1.eq.ist2).and. &
           (r3taxi(vkl(1,ik1),vklnr(1,ik2)).lt.epslat)) then
            t2=(0.5d0/occmax)*wkpt(ik1)*t3
          else
            if (t3.gt.0.d0) then
              t2=t1*wkpt(ik1)*t3**rdmalpha
            else
              t2=0.d0
            end if
          end if
        end if
        engyx=engyx-t2*vnlrdm(ist1,ik1,ist2,ik2)
      end do
    end do
  end do
end do
return
end subroutine

