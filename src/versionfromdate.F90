! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public
! License. See the file COPYING for license details.

!> Set code version, githash and compiler.
!> 
!> Note, this routine must be retained. If version.inc is directly included
!> in mod_misc, the whole code will rebuild each time `make` is issued.
subroutine versionfromdate
#include "version.inc"
    use mod_misc
    implicit none
    version=(VERSIONFROMDATE)
    githash=GITHASH//GITHASH2
    compiler_version=COMPILERVERSION
end subroutine