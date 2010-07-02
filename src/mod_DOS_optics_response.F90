
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_DOS_optics_response
!----------------------------------------------------------!
!     density of states, optics and response variables     !
!----------------------------------------------------------!
! number of energy intervals in the DOS/optics function
!replaced by inputstructureinteger::nwdos
! effective size of k/q-point grid for integrating the Brillouin zone
!replaced by inputstructureinteger::ngrdos
! smoothing level for DOS/optics function
!replaced by inputstructureinteger::nsmdos
! energy interval for DOS/optics function
      Real (8) :: wdos (2)
! scissors correction
!replaced by inputstructure!replaced by inputstructure!replaced by inputstructurereal(8)::scissor
! number of optical matrix components required
      Integer :: noptcomp
! required optical matrix components
!replaced by inputstructureinteger::optcomp(3, 27)
! usegdft is .true. if the generalised DFT correction is to be used
!replaced by inputstructurelogical::usegdft
! intraband is .true. if the intraband term is to be added to the optical matrix
!replaced by inputstructurelogical::intraband
! lmirep is .true. if the (l,m) band characters should correspond to the
! irreducible representations of the site symmetries
! Lorentzian lineshape in optics
      Logical :: optltz
! broadening for Lorentzian lineshape
      Real (8) :: optswidth
!replaced by inputstructurelogical::lmirep
! spin-quantisation axis in Cartesian coordinates used when plotting the
! spin-resolved DOS (z-axis by default)
!replaced by inputstructurereal(8)::sqados(3)
! q-vector in lattice coordinates for calculating the matrix elements
! < i,k+q | exp(iq.r) | j,k >
!replaced by inputstructurereal(8)::vecql(3)
End Module
!
