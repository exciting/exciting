! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !MODULE: modxas
! !DESCRIPTION:
!   Contains additional global variables required for XAS calculations within the BSE EXCITING code.
!
! !REVISION HISTORY:
!
!   Created JUNE 2015 by Christian Vorwerk
!EOP   
!BOC
#include "maxdefinitions.inc"
module modxas
!----------------------------!
!     Core states            !
!----------------------------!
! WARNING: the definitions here are not identical to the one in modgw!
! See xasinit.f90 for how core treatment is done in XAS-calculations
! maximum number of core states per atom
      integer(4) :: ncmax
! Max. num of core states including lm
      integer(4) :: nclm 
! Max. L of core states
      integer(4) :: lcoremax      
! Total num of core states over all atoms including lm
      integer(4) :: ncg  
! number of core states of given species and atom
      integer(4) :: ncore  
! number of initial states for XAS calculation
       integer(4) :: nxas
! radial wavefunctions for given kappa and m
      real(8), allocatable :: ucore(:,:)
!  core energies for given kappa and m
	  real(8), allocatable :: ecore(:)
! unique lowest and highest core state for XAS calculation
    integer(4) :: xasstart, xasstop
! l-value of the considered XAS edge
    integer(4) :: lxas
! j-value of given core state
    real(8), allocatable :: spj(:)
! mj-value of a given core state
    real(8), allocatable :: mj(:)

!----------------------------!
!     general                !
!----------------------------!
! shortcut for atomic position array
      real(8) :: atposl(3,_MAXATOMS_,_MAXSPECIES_)
! muffin-tin volume, relative to cell volume
      real(8) :: vmt(_MAXSPECIES_)
!! plane wave matrix elements array (o-o part)
!      Complex (8), Allocatable :: xioo (:, :, :)
!! plane wave matrix elements array (u-u part)
!      Complex (8), Allocatable :: xiuu (:, :, :)
contains
  integer function mj2ml(m,i)
    implicit none
    real(8), intent(in) :: m
    integer, intent(in) :: i
    if (i==1) then
      mj2ml=Int(m-0.50d0)
    else if (i==2) then
      mj2ml=Int(m+0.50d0)
    end if
    if (mj2ml< -lxas .or. mj2ml> lxas) then
      mj2ml=0.0d0
      !write(*,*) 'MJ= ', m, 'inconsistent with lxas= ', lxas
    end if
    return
  end function mj2ml

  real function preml(l,j,m,i)
    implicit none
    integer, intent(in) :: l, i
    real(8), intent(in) :: j, m

    if (j==dble(l)+0.50d0) then ! j=l+1/2
      if (i==1) then
        preml=sqrt((dble(l)+m+0.50d0)/(2*dble(l)+1.0d0))
      else if (i==2) then
        preml=sqrt((dble(l)-m+0.50d0)/(2*dble(l)+1.0d0))
      end if
    else if (j==dble(l)-0.50d0) then ! j=l-1/2
      if (i==1) then
        preml=-sqrt((dble(l)-m+0.50d0)/(2*dble(l)+1.0d0))
      else if (i==2) then
        preml=sqrt((dble(l)+m+0.50d0)/(2*dble(l)+1.0d0))
      end if
    end if
    return
  end function preml
end module modxas
!EOC
