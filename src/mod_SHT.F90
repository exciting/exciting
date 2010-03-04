
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
#include "maxdefinitions.inc"
Module mod_SHT
!-----------------------------------------------------!
!     spherical harmonic transform (SHT) matrices     !
!-----------------------------------------------------!
! real backward SHT matrix for lmaxapw
      Real (8), Allocatable :: rbshtapw (:, :)
! real forward SHT matrix for lmmaxapw
      Real (8), Allocatable :: rfshtapw (:, :)
! real backward SHT matrix for lmaxvr
      Real (8), Allocatable :: rbshtvr (:, :)
! real forward SHT matrix for lmaxvr
      Real (8), Allocatable :: rfshtvr (:, :)
! complex backward SHT matrix for lmaxapw
      Complex (8), Allocatable :: zbshtapw (:, :)
! complex forward SHT matrix for lmaxapw
      Complex (8), Allocatable :: zfshtapw (:, :)
! complex backward SHT matrix for lmaxvr
      Complex (8), Allocatable :: zbshtvr (:, :)
! complex forward SHT matrix for lmaxvr
      Complex (8), Allocatable :: zfshtvr (:, :)
End Module
!
