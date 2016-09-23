
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
      Complex (8), Allocatable :: haaij(:,:,:)
      integer haaijSize
! local-orbital-APW Hamiltonian integrals
      Real (8), Allocatable :: hloa (:, :, :, :, :)
      Complex (8), Allocatable :: haloij(:,:,:)
      integer, allocatable :: haloijSize(:)
! local-orbital-local-orbital Hamiltonian integrals
      Real (8), Allocatable :: hlolo (:, :, :, :)
      Complex (8), Allocatable :: hloloij(:,:,:)
! The stuff for the linearized Koelling-Harmon
! APW-APW Hamiltonian integrals
      Real (8), Allocatable :: h1aa (:, :, :, :)
! local-orbital-APW Hamiltonian integrals
      Real (8), Allocatable :: h1loa (:, :, :)
! local-orbital-local-orbital Hamiltonian integrals
      Real (8), Allocatable :: h1lolo (:, :, :)
      logical h1on
! complex Gaunt coefficient array
      Complex (8), Allocatable :: gntyry (:, :, :),gntryy (:, :, :),gntnonz(:)
! list of non-zero Gaunt coefficients
      Integer, Allocatable :: gntnonzlm1(:),gntnonzlm2(:),gntnonzlm3(:),gntnonzlindex(:),gntnonzl2index(:,:)
!gntnonzlist (:, :, :)
! tseqit is .true. if the first-variational secular equation is to be solved
! iteratively
      Logical :: tseqit
! number of secular equation iterations per self-consistent loop
      Integer :: nseqit
! iterative solver step length
      Real (8) :: tauseq
! ARPACK seed vector
      complex(8),allocatable :: arpackseed(:,:)
End Module
!
