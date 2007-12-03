
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine writeengy(fnum)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
write(fnum,*)
write(fnum,'("Energies :")')
write(fnum,'(" Fermi                   : ",G18.10)') efermi
write(fnum,'(" sum of eigenvalues      : ",G18.10)') evalsum
write(fnum,'(" electronic kinetic      : ",G18.10)') engykn
write(fnum,'(" core electron kinetic   : ",G18.10)') engykncr
write(fnum,'(" Coulomb                 : ",G18.10)') engycl
write(fnum,'(" Coulomb potential       : ",G18.10)') engyvcl
write(fnum,'(" nuclear-nuclear         : ",G18.10)') engynn
write(fnum,'(" electron-nuclear        : ",G18.10)') engyen
write(fnum,'(" Hartree                 : ",G18.10)') engyhar
write(fnum,'(" Madelung                : ",G18.10)') engymad
if (chgexs.ne.0.d0) then
  write(fnum,'(" comp. background charge : ",G18.10)') engycbc
end if
write(fnum,'(" xc potential            : ",G18.10)') engyvxc
if (spinpol) then
  write(fnum,'(" xc effective B-field    : ",G18.10)') engybxc
  write(fnum,'(" external B-field        : ",G18.10)') engybext
end if
write(fnum,'(" exchange                : ",G18.10)') engyx
write(fnum,'(" correlation             : ",G18.10)') engyc
if (ldapu.ne.0) then
  write(fnum,'(" LDA+U                   : ",G18.10)') engylu
end if
write(fnum,'(" total energy            : ",G18.10)') engytot
if (spinpol) then
  write(fnum,'(" (external B-field energy excluded from total)")')
end if
!call flushifc(fnum)
return
end subroutine

