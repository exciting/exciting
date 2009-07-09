


! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine putevecfv(ik, evecfv)
  use modmain
  use modmpi

  implicit none
  ! arguments
  integer, intent(in) :: ik
  complex(8), intent(in) :: evecfv(nmatmax, nstfv, nspnfv)
  ! local variables

  character(256) ::filetag
  character(256), external:: outfilenamestring
  integer::recl, koffset


  ! find the record length
  inquire(iolength=recl) vkl(:, ik), nmatmax, nstfv, nspnfv, evecfv
  !$OMP CRITICAL
  filetag='EVECFV'
  if (splittfile.or.(rank.eq.0)) then
  open(70, file = outfilenamestring(filetag, ik), action = 'WRITE', &
       form = 'UNFORMATTED', access = 'DIRECT', recl = recl)
if (splittfile) then
 koffset=ik-firstk(procofk(ik))+1
 else
 koffset =ik
 endif
write(70, rec=koffset)vkl(:, ik), nmatmax, nstfv, nspnfv, evecfv
  close(70)

  endif
   !$OMP END CRITICAL
  return
end subroutine putevecfv
