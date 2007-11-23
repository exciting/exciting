
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfpack(tpack,n,lrstp,rfmt,rfir,nu)
use modmain
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(inout) :: n
integer, intent(in) :: lrstp
real(8), intent(inout) :: rfmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(inout) :: rfir(ngrtot)
real(8), intent(out) :: nu(*)
! local variables
integer is,ia,ias,ir,lm
if (tpack) then
! pack the function
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ir=1,nrmt(is),lrstp
        do lm=1,lmmaxvr
          n=n+1
          nu(n)=rfmt(lm,ir,ias)
        end do
      end do
    end do
  end do
  do ir=1,ngrtot
    n=n+1
    nu(n)=rfir(ir)
  end do
else
! unpack the function
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ir=1,nrmt(is),lrstp
        do lm=1,lmmaxvr
          n=n+1
          rfmt(lm,ir,ias)=nu(n)
        end do
      end do
    end do
  end do
  do ir=1,ngrtot
    n=n+1
    rfir(ir)=nu(n)
  end do
end if
return
end subroutine
