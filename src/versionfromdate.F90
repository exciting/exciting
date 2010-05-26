
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public
! License. See the file COPYING for license details.

subroutine versionfromdate
#include "version.inc"
  use mod_misc
  implicit none
  version=(VERSIONFROMDATE)
  githash=GITHASH//GITHASH2
end subroutine
