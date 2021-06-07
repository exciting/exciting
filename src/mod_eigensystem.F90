
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!> Overlap and Hamiltonian variables
Module mod_eigensystem
   implicit none
   !> order of overlap and Hamiltonian matrices for each k-point
   Integer, pointer :: nmat_ptr(:, :)
   Integer, Allocatable, target :: nmat(:, :)
   !> maximum nmat over all k-points
   integer, pointer :: nmatmax_ptr
   Integer, target :: nmatmax
   !> size of packed matrices
   Integer, Allocatable :: npmat(:, :)
   !> index to the position of the local-orbitals in the H and O matrices
   Integer, Allocatable :: idxlo(:, :, :)
   !> APW-local-orbital overlap integrals
   Real(8), Allocatable :: oalo(:, :, :)
   !> local-orbital-local-orbital overlap integrals
   Real(8), Allocatable :: ololo(:, :, :)
   !> APW-APW Hamiltonian integrals
   Real(8), Allocatable :: haa(:, :, :, :, :, :)
   Complex(8), Allocatable :: haaij(:, :, :)
   integer :: haaijSize
   !> local-orbital-APW Hamiltonian integrals
   Real(8), Allocatable :: hloa(:, :, :, :, :)
   Complex(8), Allocatable :: haloij(:, :, :)
   integer, allocatable :: haloijSize(:)
   !> local-orbital-local-orbital Hamiltonian integrals
   Real(8), Allocatable :: hlolo(:, :, :, :)
   Complex(8), Allocatable :: hloloij(:, :, :)
   ! The stuff for the linearized Koelling-Harmon
   !> APW-APW Hamiltonian integrals
   Real(8), Allocatable :: h1aa(:, :, :, :)
   !> local-orbital-APW Hamiltonian integrals
   Real(8), Allocatable :: h1loa(:, :, :)
   !> local-orbital-local-orbital Hamiltonian integrals
   Real(8), Allocatable :: h1lolo(:, :, :)
   logical :: h1on
   !> complex Gaunt coefficient array
   Complex(8), Allocatable :: gntyry(:, :, :), gntryy(:, :, :), gntnonz(:)
   !> list of non-zero Gaunt coefficients
   Integer, Allocatable :: gntnonzlm1(:), gntnonzlm2(:), gntnonzlm3(:), gntnonzlindex(:), gntnonzl2index(:, :)
   !> another compact list of Gaunt coefficients
   complex(8), allocatable :: listgnt(:, :, :)
   integer, allocatable :: indgnt(:, :, :)
   !gntnonzlist (:, :, :)
   !> tseqit is .true. if the first-variational secular equation is to be solved iteratively
   Logical :: tseqit
   !> number of secular equation iterations per self-consistent loop
   Integer :: nseqit
   !> iterative solver step length
   Real(8) :: tauseq
   !> ARPACK seed vector
   complex(8), allocatable :: arpackseed(:, :)

   !> Matrix-elements for muffin-tin functions
   Type MTHamiltonianType
      complex(8), pointer :: aa(:, :, :), alo(:, :, :), loa(:, :, :), lolo(:, :, :)
   End Type MTHamiltonianType

   !>  Muffin-tin Hamiltonian components 
   Type MTHamiltonianList
      Type(MTHamiltonianType) :: main
      Type(MTHamiltonianType) :: spinless
      Type(MTHamiltonianType) :: alpha
      Type(MTHamiltonianType) :: beta
      Type(MTHamiltonianType) :: ab
      Type(MTHamiltonianType) :: ba
      integer :: maxnlo
      integer :: maxaa
      integer, allocatable :: losize(:)
      contains
            procedure :: release => MTRelease    !! Release all memory 
   End Type MTHamiltonianList

   !> Muffin-tin Hamlitonian for SCF
   Type(MTHamiltonianList) :: mt_hscf

   !> Relativity settings
   integer :: level_nr, level_zora, level_iora
   parameter(level_nr=0, level_zora=1, level_iora=2)

Contains

!
!
!
!BOP
! !ROUTINE: MTNullify
! !INTERFACE:
!
!
   subroutine MTNullify(mt_h)
! !USES:
!      Use modinput
!      Use modmain
! !DESCRIPTION:
! Initialises all parts of the muffin-tin Hamiltonian or overlap.
!
! !REVISION HISTORY:
!   Created June 2019 (Andris)
!   Modified Dezember 2020 (Ronaldo)
!EOP
!BOC
      Implicit None
      Type(MTHamiltonianList), Intent(Inout) :: mt_h

      nullify (mt_h%main%aa)
      nullify (mt_h%main%alo)
      nullify (mt_h%main%loa)
      nullify (mt_h%main%lolo)

      nullify (mt_h%spinless%aa)
      nullify (mt_h%spinless%alo)
      nullify (mt_h%spinless%loa)
      nullify (mt_h%spinless%lolo)

      nullify (mt_h%alpha%aa)
      nullify (mt_h%alpha%alo)
      nullify (mt_h%alpha%loa)
      nullify (mt_h%alpha%lolo)

      nullify (mt_h%beta%aa)
      nullify (mt_h%beta%alo)
      nullify (mt_h%beta%loa)
      nullify (mt_h%beta%lolo)

      nullify (mt_h%ab%aa)
      nullify (mt_h%ab%alo)
      nullify (mt_h%ab%loa)
      nullify (mt_h%ab%lolo)

      nullify (mt_h%ba%aa)
      nullify (mt_h%ba%alo)
      nullify (mt_h%ba%loa)
      nullify (mt_h%ba%lolo)

   end subroutine MTNullify

!
!
!
!BOP
! !FUNCTION: MaxAPWs
! !INTERFACE:
!
!
   integer Function MaxAPWs()
! !USES:
      Use modinput
      Use mod_APW_LO
      Use mod_atoms
      Use mod_muffin_tin
! !DESCRIPTION:
! Finds the maximum number of APWs across all muffin-tins.
!
! !REVISION HISTORY:
!   Created January 2021 (Andris)
!EOP
!BOC
      Implicit None
      Type(MTHamiltonianList) :: mt_h
      integer ::  io, if1, l, m, lm, maxsize, is, ias

      maxsize = 0
      Do is = 1, nspecies
         if1 = 0
         Do l = 0, input%groundstate%lmaxmat
            Do m = -l, l
               lm = idxlm(l, m)
               Do io = 1, apword(l, is)
                  if1 = if1 + 1
               End Do
            End Do
         End Do
         maxsize = max(if1, maxsize)
      End do

      MaxAPWs = maxsize

   end function MaxAPWs

!
!
!
!BOP
! !ROUTINE: MTInitAll
! !INTERFACE:
!
!
   subroutine MTInitAll(mt_h)
! !USES:
      Use modinput
      Use mod_APW_LO
      Use mod_atoms
      Use mod_muffin_tin
! !DESCRIPTION:
! Initialises all parts of the muffin-tin Hamiltonian or overlap.
!
! !REVISION HISTORY:
!   Created October 2015 (Andris)
!EOP
!BOC
      Implicit None
      Type(MTHamiltonianList) :: mt_h
      integer ::  io, ilo, if1, l, m, lm, l1, m1, lm1, l3, m3, lm3, is, ias
!      integer, external :: MaxAPWs

      mt_h%maxaa = MaxAPWs()

      if (allocated(mt_h%losize)) deallocate (mt_h%losize)
      allocate (mt_h%losize(nspecies))
      mt_h%maxnlo = 0
      Do is = 1, nspecies
         ias = idxas(1, is)
         ilo = nlorb(is)
         if (ilo .gt. 0) then
            l1 = lorbl(ilo, is)
            lm1 = idxlm(l1, l1)
            l3 = lorbl(1, is)
            lm3 = idxlm(l3, -l3)
            mt_h%losize(is) = idxlo(lm1, ilo, ias) - idxlo(lm3, 1, ias) + 1
            if (mt_h%maxnlo .lt. mt_h%losize(is)) mt_h%maxnlo = mt_h%losize(is)
         else
            mt_h%losize(is) = 0
         end if
      End do

      call MTInit(mt_h%spinless, mt_h%maxaa, mt_h%maxnlo)
!      call MTRedirect(mt_h%main,mt_h%alpha)

   end subroutine MTInitAll

!
!
!
!BOP
! !ROUTINE: MTInit
! !INTERFACE:
!
!
   subroutine MTInit(mt_block, maxaa, maxnlo)
! !USES:
      Use mod_atoms
!      Use modmain
! !DESCRIPTION:
! Initialises a part of the muffin-tin Hamiltonian or overlap.
!
! !REVISION HISTORY:
!   Created October 2015 (Andris)
!EOP
!BOC
      Implicit None
      Type(MTHamiltonianType) :: mt_block
      integer, intent(in) :: maxaa, maxnlo

      if (.not. associated(mt_block%aa)) then
         allocate (mt_block%aa(maxaa, maxaa, natmtot))
         mt_block%aa = 0d0
         if (maxnlo .gt. 0) then
            allocate (mt_block%loa(maxnlo, maxaa, natmtot))
            allocate (mt_block%alo(maxaa, maxnlo, natmtot))
            allocate (mt_block%lolo(maxnlo, maxnlo, natmtot))
            mt_block%loa = 0d0
            mt_block%alo = 0d0
            mt_block%lolo = 0d0
         end if
      else
         mt_block%aa = 0d0
         if (maxnlo .gt. 0) then
            mt_block%loa = 0d0
            mt_block%alo = 0d0
            mt_block%lolo = 0d0
         end if

      end if

   end subroutine MTinit

!
!
!
!BOP
! !ROUTINE: MTRedirect
! !INTERFACE:
!
!
   subroutine MTRedirect(mt_blockA, mt_blockB)
! !USES:
      Use mod_atoms
!      Use modmain
! !DESCRIPTION:
! Redirects pointers of mt\_blockA to mt\_blockB of muffin-tin Hamiltonians or overlaps.
!
! !REVISION HISTORY:
!   Created October 2015 (Andris)
!EOP
!BOC
      Implicit None
      Type(MTHamiltonianType), intent(in)  :: mt_blockB
      Type(MTHamiltonianType), intent(out) :: mt_blockA

      if (associated(mt_blockB%aa)) mt_blockA%aa => mt_blockB%aa
      if (associated(mt_blockB%alo)) mt_blockA%alo => mt_blockB%alo
      if (associated(mt_blockB%loa)) mt_blockA%loa => mt_blockB%loa
      if (associated(mt_blockB%lolo)) mt_blockA%lolo => mt_blockB%lolo

   end subroutine MTRedirect
!
!
!
!BOP
! !ROUTINE: MTCopy
! !INTERFACE:
!
!
   subroutine MTCopy(mt_blockA, mt_blockB)
! !USES:
      use mod_atoms
      use modmpi, only: terminate_mpi_env, mpiglobal

!      Use modmain
! !DESCRIPTION:
! Copies the contents of mt\_blockA to mt\_blockB of muffin-tin Hamiltonians or overlaps.
!
! !REVISION HISTORY:
!   Created October 2015 (Andris)
!EOP
!BOC
      Implicit None
      Type(MTHamiltonianType) :: mt_blockA, mt_blockB
      integer :: maxaa, maxnlo

      ! Check if mt_blockA%aa is associated to some target
      if (.not. associated(mt_blockA%aa)) then
         call terminate_mpi_env(mpiglobal, &
           & 'Error (MTCopy): mt_blockA%aa is not allocated')
      end if
      ! Check if mt_blockB%aa is associated to some target
      if (.not. associated(mt_blockB%aa)) then
         call terminate_mpi_env(mpiglobal, &
           & 'Error (MTCopy): mt_blockB%aa is not allocated')
      end if

      ! Now we can safely make mt_blockB%aa points to the same target as mt_blockA%aa
      mt_blockB%aa = mt_blockA%aa

      ! We will repeat this procedure below

      ! This first check is to verify if there are any local orbitals
      if (associated(mt_blockA%lolo)) then
         if (.not. associated(mt_blockA%alo)) then
            call terminate_mpi_env(mpiglobal, &
              & 'Error (MTCopy): mt_blockA%alo is not allocated')
         end if
         if (.not. associated(mt_blockB%alo)) then
            call terminate_mpi_env(mpiglobal, &
              & 'Error (MTCopy): mt_blockB%alo is not allocated')
         end if
         mt_blockB%alo = mt_blockA%alo

         if (.not. associated(mt_blockA%loa)) then
            call terminate_mpi_env(mpiglobal, &
              & 'Error (MTCopy): mt_blockA%loa is not allocated')
         end if
         if (.not. associated(mt_blockB%loa)) then
            call terminate_mpi_env(mpiglobal, &
              & 'Error (MTCopy): mt_blockB%loa is not allocated')
         end if
         mt_blockB%loa = mt_blockA%loa
         if (.not. associated(mt_blockB%lolo)) then
            call terminate_mpi_env(mpiglobal, &
              & 'Error (MTCopy): mt_blockB%lolo is not allocated')
         end if
         mt_blockB%lolo = mt_blockA%lolo
      end if
   end subroutine MTCopy


   !> Release all memory used by the muffin-tin Hamiltonian or overlap
   subroutine MTRelease(this)
      implicit none
      class(MTHamiltonianList), Intent(Inout) :: this

      if (associated(this%spinless%aa)) deallocate (this%spinless%aa)
      if (associated(this%spinless%alo)) deallocate (this%spinless%alo)
      if (associated(this%spinless%loa)) deallocate (this%spinless%loa)
      if (associated(this%spinless%lolo)) deallocate (this%spinless%lolo)

      if (associated(this%alpha%aa)) deallocate (this%alpha%aa)
      if (associated(this%alpha%alo)) deallocate (this%alpha%alo)
      if (associated(this%alpha%loa)) deallocate (this%alpha%loa)
      if (associated(this%alpha%lolo)) deallocate (this%alpha%lolo)

      if (associated(this%beta%aa)) deallocate (this%beta%aa)
      if (associated(this%beta%alo)) deallocate (this%beta%alo)
      if (associated(this%beta%loa)) deallocate (this%beta%loa)
      if (associated(this%beta%lolo)) deallocate (this%beta%lolo)

      if (associated(this%ab%aa)) deallocate (this%ab%aa)
      if (associated(this%ab%alo)) deallocate (this%ab%alo)
      if (associated(this%ab%loa)) deallocate (this%ab%loa)
      if (associated(this%ab%lolo)) deallocate (this%ab%lolo)

      if (associated(this%ba%aa)) deallocate (this%ba%aa)
      if (associated(this%ba%alo)) deallocate (this%ba%alo)
      if (associated(this%ba%loa)) deallocate (this%ba%loa)
      if (associated(this%ba%lolo)) deallocate (this%ba%lolo)

      if (allocated(this%losize)) deallocate (this%losize)

   end subroutine MTRelease

End Module
