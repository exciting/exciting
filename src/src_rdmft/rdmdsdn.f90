
! Copyright (C) 2008 T. Baldsiefen, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine rdmdsdn(dedn)
use modmain
implicit none
! arguments
real(8), intent(inout) :: dedn(nstsv,nkpt)
! local variables
integer  ik,ist
real(8)  t1
do ik=1,nkpt
  do ist=1,nstsv
    t1=max(occsv(ist,ik),epsocc)
    t1=min(t1,occmax-epsocc)
    dedn(ist,ik)=dedn(ist,ik)-rdmtemp*kboltz*log(t1/(occmax-t1))
  end do
end do
return
end subroutine

