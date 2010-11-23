
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_eigenvalue_occupancy
      Use modinput
!--------------------------------------------!
!     eigenvalue and occupancy variables     !
!--------------------------------------------!
! number of empty states
!replaced by inputstructureinteger::nempty
! number of first-variational states
      Integer :: nstfv
! number of second-variational states
      Integer :: nstsv
! smearing type
!replaced by inputstructureinteger::stype
! smearing function description
      Character (256) :: sdescr
! smearing width
!replaced by inputstructurereal(8)::swidth
! maximum allowed occupancy (1 or 2)
      Real (8) :: occmax
! convergence tolerance for occupancies
!replaced by inputstructurereal(8)::epsocc
! second-variational occupation number array
      Real (8), Allocatable :: occsv (:, :)
! Fermi energy for second-variational states
      Real (8) :: efermi
! density of states at the Fermi energy
      Real (8) :: fermidos
! error tolerance for the first-variational eigenvalues
!replaced by inputstructurereal(8)::evaltol
! minimum allowed eigenvalue
!replaced by inputstructurereal(8)::evalmin
! second-variational eigenvalues
      Real (8), Allocatable :: evalsv (:, :)
! tevecsv is .true. if second-variational eigenvectors are calculated
!replaced by inputstructurelogical::tevecsv
End Module
!
