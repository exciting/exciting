
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine writechg(fnum)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
integer is,ia,ias
! output charges
write(fnum,*)
write(fnum,'("Charges :")')
write(fnum,'(" core                    : ",G18.10)') chgcr
write(fnum,'(" core leakage            : ",G18.10)') chgcrlk
write(fnum,'(" valence                 : ",G18.10)') chgval
write(fnum,'(" interstitial            : ",G18.10)') chgir
write(fnum,'(" muffin-tins")')
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(fnum,'("  species ",I4," atom ",I4," : ",G18.10)') is,ia,chgmt(ias)
  end do
end do
write(fnum,'(" total in muffin-tins    : ",G18.10)') chgmttot
write(fnum,'(" total charge            : ",G18.10)') chgcalc
! output moments
if (spinpol) then
  write(fnum,*)
  write(fnum,'("Moments :")')
  write(fnum,'(" interstitial            : ",3G18.10)') momir(1:ndmag)
  write(fnum,'(" muffin-tins")')
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(fnum,'("  species ",I4," atom ",I4," : ",3G18.10)') is,ia, &
       mommt(1:ndmag,ias)
    end do
  end do
  write(fnum,'(" total in muffin-tins    : ",3G18.10)') mommttot(1:ndmag)
  write(fnum,'(" total moment            : ",3G18.10)') momtot(1:ndmag)
end if
call flushifc(fnum)
return
end subroutine

