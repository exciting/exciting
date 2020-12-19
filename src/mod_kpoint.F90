
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!>  k-point set variables  
Module mod_kpoint
      implicit none 
! autokpt is .true. if the k-point set is determined automatically
!replaced by inputstructurelogical::autokpt
! radius of sphere used to determine k-point density when autokpt is .true.
!replaced by inputstructurereal(8)::radkpt
! k-point grid sizes
!replaced by inputstructureinteger::ngridk(3)
! total number of k-points
      integer, pointer :: nkpt_ptr
      Integer, target :: nkpt
! k-point offset
!replaced by inputstructurereal(8)::vkloff(3)
! reducek is .true. if k-points are to be reduced (with crystal symmetries)
!replaced by inputstructurelogical::reducek
! locations of k-points on integer grid
      Integer, Allocatable :: ivk (:, :)
! k-points in lattice coordinates
      real(8), pointer :: vkl_ptr(:,:)
      Real (8), Allocatable, target :: vkl(:,:)
! k-points in Cartesian coordinates
      Real (8), Allocatable :: vkc (:, :)
! k-point weights
      Real (8), Allocatable :: wkpt (:)
! map from non-reduced grid to reduced set
      Integer, Allocatable :: ikmap (:, :, :)
! total number of non-reduced k-points
      Integer :: nkptnr
! locations of non-reduced k-points on integer grid
      Integer, Allocatable :: ivknr (:, :)
! non-reduced k-points in lattice coordinates
      Real (8), Allocatable :: vklnr (:, :)
! non-reduced k-points in Cartesian coordinates
      Real (8), Allocatable :: vkcnr (:, :)
! non-reduced k-point weights
      Real (8), Allocatable :: wkptnr (:)
! map from non-reduced grid to non-reduced set
      Integer, Allocatable :: ikmapnr (:, :, :)
! k-point at which to determine effective mass tensor
!replaced by inputstructurereal(8)::vklem(3)
! displacement size for computing the effective mass tensor
!replaced by inputstructurereal(8)::deltaem
! number of displacements in each direction
!replaced by inputstructureinteger::ndspem
  !--------------------------------------!
  !     tetrahedron method variables     !
  !--------------------------------------!
  ! integer k-point offset
      Integer (4) :: ikloff (3)
  ! k-point offset divisor
      Integer (4) :: dkloff
  ! k-points common divisor
      Integer (4) :: dvk
  ! Number of tetrahedra
      Integer (4) :: ntet
  ! index of the k-points corresponding to the nodes of each tetrahedron
      Integer (4), Allocatable :: tnodes (:, :)
  ! weight of each tetrahedron.
      Integer (4), Allocatable :: wtet (:)
  ! volume of the tetrahedra relative to the BZ volume
      Real (8) :: tvol
  ! parameter specifying smalles diagonal in generic tetrahedron
      Integer (4) :: mnd
  ! mapping reducible -> irreducible k-point index
      integer(4), allocatable :: ik2ikp(:)
  ! mapping irreducible -> reducible k-point index
      integer(4), allocatable :: ikp2ik(:)
  ! integer weight of the irreducible k-point            
      integer(4), allocatable :: iwkp(:)
!
Contains
!
!BOP
! !ROUTINE: rtorat
! !INTERFACE:
      Subroutine rtorat (eps, n, x, k, div)
! !DESCRIPTION:
!   This subroutine factorizes the real coordinates of a vector {\bf x}.
!   The output is an integer vector {\bf k}, such that
!   $$ |x(i)-k(i)/{\rm div}| < {\rm eps} $$
!   for all $i=1,\ldots,n$.
!
! !REVISION HISTORY:
!   Created July 2008 by Sagmeister
!EOP
!BOC
         Implicit None
    ! arguments
         Real (8), Intent (In) :: eps
         Integer (4), Intent (In) :: n
         Real (8), Intent (In) :: x (n)
         Integer (4), Intent (Out) :: div
         Integer (4), Intent (Out) :: k (n)
    ! local variables
         Integer :: maxint
         Real (8) :: dx
         maxint = Nint (1.d0/eps)
         Do div = 1, maxint
            k (:) = Nint (dble(div)*x(:))
            dx = maxval (Abs(dble(k)/dble(div)-x))
            If (dx .Lt. eps) Exit
         End Do
         If (dx .Ge. eps) Then
            Write (*,*)
            Write (*, '("Error(modtetra:rtorat): factorization failed")&
           &')
            Write (*, '(" maximum integer :",i12)') maxint
            Write (*, '(" tolerance       :",g18.10)') eps
            Write (*, '(" deviation       :",g18.10)') dx
            Write (*,*)
            Stop
         End If
      End Subroutine rtorat
!EOC
End Module
