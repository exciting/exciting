
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_symmetry
!----------------------------!
!     symmetry variables     !
!----------------------------!
! nosym is .true. if no symmetry information should be used
!replaced by inputstructurelogical::nosym
! number of Bravais lattice point group symmetries
      Integer :: nsymlat
! Bravais lattice point group symmetries
      Integer :: symlat (3, 3, 48)
! determinants of lattice symmetry matrices (1 or -1)
      Integer :: symlatd (48)
! index to inverses of the lattice symmetries
      Integer :: isymlat (48)
! lattice point group symmetries in Cartesian coordinates
      Real (8) :: symlatc (3, 3, 48)
! tshift is .true. if atomic basis is allowed to be shifted
!replaced by inputstructurelogical::tshift
! maximum of symmetries allowed
      Integer, Parameter :: maxsymcrys = 192
! number of crystal symmetries
      Integer :: nsymcrys
! crystal symmetry translation vector in lattice coordinates
      Real (8) :: vtlsymc (3, maxsymcrys)
! spatial rotation element in lattice point group for each crystal symmetry
      Integer :: lsplsymc (maxsymcrys)
! global spin rotation element in lattice point group for each crystal symmetry
      Integer :: lspnsymc (maxsymcrys)
! equivalent atom index for each crystal symmetry
      Integer, Allocatable :: ieqatom (:, :, :)
! eqatoms(ia,ja,is) is .true. if atoms ia and ja are equivalent
      Logical, Allocatable :: eqatoms (:, :, :)
! number of site symmetries
      Integer, Allocatable :: nsymsite (:)
! site symmetry spatial rotation element in lattice point group
      Integer, Allocatable :: lsplsyms (:, :)
! site symmetry global spin rotation element in lattice point group
      Integer, Allocatable :: lspnsyms (:, :)
End Module
!
