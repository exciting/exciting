


! Copyright (C) 2008 T. Baldsiefen, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.


subroutine rdmentropy
use modmain
use modinput
implicit none
! local variables
integer::ik, ist
real(8)::t1
rdmentrpy=0.d0
do ik=1, nkpt
  do ist=1, nstsv
    t1=max(occsv(ist, ik), input%groundstate%epsocc)
    t1=min(t1, occmax-input%groundstate%epsocc)
    rdmentrpy = rdmentrpy - wkpt(ik) * (t1 * log(t1/occmax) &
     +(occmax - t1) * log(1.d0 - t1/occmax))
  end do
end do	
rdmentrpy=kboltz*rdmentrpy
return
end subroutine
