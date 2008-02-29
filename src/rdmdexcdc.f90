
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rdmdexcdc(ikp,evecsv,dedc)
! calculate the derivative of exchange-correlation energy w.r.t. evecsv
use modmain
implicit none
! arguments
integer, intent(in) :: ikp
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(inout) :: dedc(nstsv,nstsv)
! local variables
integer ik,jk,iv(3)
integer ist1,ist2,ist3,ist4
real(8) t1,t2
! allocatable arrays
complex(8), allocatable :: vnl(:,:,:,:)
! external functions
real(8) r3taxi
external r3taxi
if (rdmxctype.eq.0) return
! calculate the prefactor
if (rdmxctype.eq.1) then
  t1=1.d0/occmax
else if (rdmxctype.eq.2) then
  if (spinpol) then
    t1=1.d0
  else
    t1=2.d0*(0.25d0)**rdmalpha
  end if
else
  write(*,*)
  write(*,'("Error(rdmdexcdc): rdmxctype not defined : ",I8)') rdmxctype
  write(*,*)
  stop
end if
allocate(vnl(nstsv,nstsv,nstsv,nkptnr))
! calculate non-local matrix elements of the type (l-jj-k)
call rdmvnlc(ikp,vnl)
! start loop over non-reduced k-points
do ik=1,nkptnr
! copy the matrix elements of the type i-jj-i to vnlrdm
  do ist1=1,nstsv
    do ist2=1,nstsv
      vnlrdm(ist1,ikp,ist2,ik)=dble(vnl(ist1,ist1,ist2,ik))
    end do
  end do
! find the equivalent reduced k-point
  iv(:)=ivknr(:,ik)
  jk=ikmap(iv(1),iv(2),iv(3))
  do ist1=1,nstsv
    do ist2=1,nstsv
      do ist3=1,nstsv
        do ist4=1,nstsv
          if (rdmxctype.eq.1) then
! Hartree-Fock functional
            t2=t1*occsv(ist3,ikp)*occsv(ist4,jk)
          else if (rdmxctype.eq.2) then
! SDLG functional
            t2=t1*(occsv(ist3,ikp)*occsv(ist4,jk))**rdmalpha
          end if
          dedc(ist2,ist3)=dedc(ist2,ist3)-t2*evecsv(ist2,ist1)* &
           vnl(ist1,ist3,ist4,ik)
        end do
      end do
    end do
  end do
! end loop over non-reduced k-points
end do
deallocate(vnl)
return
end subroutine

