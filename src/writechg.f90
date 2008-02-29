
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
write(fnum,'(" core",T30,": ",G18.10)') chgcr
write(fnum,'(" core leakage",T30,": ",G18.10)') chgcrlk
write(fnum,'(" valence",T30,": ",G18.10)') chgval
write(fnum,'(" interstitial",T30,": ",G18.10)') chgir
write(fnum,'(" muffin-tins")')
do is=1,nspecies
  write(fnum,'("  species : ",I4," (",A,")")') is,trim(spsymb(is))
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(fnum,'("   atom ",I4,T30,": ",G18.10)') ia,chgmt(ias)
  end do
end do
write(fnum,'(" total in muffin-tins",T30,": ",G18.10)') chgmttot
if (chgexs.ne.0.d0) then
  write(fnum,'(" excess",T30,": ",G18.10)') chgexs
end if
write(fnum,'(" total charge",T30,": ",G18.10)') chgcalc
! output moments
if (spinpol) then
  write(fnum,*)
  write(fnum,'("Moments :")')
  write(fnum,'(" interstitial",T30,": ",3G18.10)') momir(1:ndmag)
  write(fnum,'(" muffin-tins")')
  do is=1,nspecies
    write(fnum,'("  species : ",I4," (",A,")")') is,trim(spsymb(is))
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(fnum,'("   atom ",I4,T30,": ",3G18.10)') ia,mommt(1:ndmag,ias)
    end do
  end do
  write(fnum,'(" total in muffin-tins",T30,": ",3G18.10)') mommttot(1:ndmag)
  write(fnum,'(" total moment",T30,": ",3G18.10)') momtot(1:ndmag)
end if
call flushifc(fnum)
return
end subroutine

