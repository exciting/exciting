
! Copyright (C) 2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ematqdealloc
  use modmain
  use modmpi
  use modxs
  implicit none
  deallocate(evecfv,evecfv0,evalsv0,xih,apwdlm,apwdlm0,lodlm,lodlm0)
end subroutine ematqdealloc
