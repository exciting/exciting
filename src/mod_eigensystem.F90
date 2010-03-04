
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_eigensystem
!-------------------------------------------!
!     overlap and Hamiltonian variables     !
!-------------------------------------------!
! order of overlap and Hamiltonian matrices for each k-point
      Integer, Allocatable :: nmat (:, :)
! maximum nmat over all k-points
      Integer :: nmatmax
! size of packed matrices
      Integer, Allocatable :: npmat (:, :)
! index to the position of the local-orbitals in the H and O matrices
      Integer, Allocatable :: idxlo (:, :, :)
! APW-local-orbital overlap integrals
      Real (8), Allocatable :: oalo (:, :, :)
! local-orbital-local-orbital overlap integrals
      Real (8), Allocatable :: ololo (:, :, :)
! APW-APW Hamiltonian integrals
      Real (8), Allocatable :: haa (:, :, :, :, :, :)
! local-orbital-APW Hamiltonian integrals
      Real (8), Allocatable :: hloa (:, :, :, :, :)
! local-orbital-local-orbital Hamiltonian integrals
      Real (8), Allocatable :: hlolo (:, :, :, :)
! complex Gaunt coefficient array
      Complex (8), Allocatable :: gntyry (:, :, :)
! tseqit is .true. if the first-variational secular equation is to be solved
! iteratively
      Logical :: tseqit
! number of secular equation iterations per self-consistent loop
      Integer :: nseqit
! iterative solver step length
      Real (8) :: tauseq
End Module
!
