
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
      Integer, pointer :: nmat_ptr(:, :)
      Integer, Allocatable, target :: nmat (:, :)
! maximum nmat over all k-points
      integer, pointer :: nmatmax_ptr
      Integer, target :: nmatmax
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
      integer :: haaijSize
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
      logical :: h1on
! complex Gaunt coefficient array
      Complex (8), Allocatable :: gntyry (:, :, :),gntryy (:, :, :),gntnonz(:)
! list of non-zero Gaunt coefficients
      Integer, Allocatable :: gntnonzlm1(:),gntnonzlm2(:),gntnonzlm3(:),gntnonzlindex(:),gntnonzl2index(:,:)
! another compact list of Gaunt coefficients
      complex(8), allocatable :: listgnt(:,:,:)
      integer, allocatable :: indgnt(:,:,:)
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
! Matrix-elements for muffin-tin functions
      Type MTHamiltonianType
        complex(8), pointer :: aa(:,:,:),alo(:,:,:),loa(:,:,:),lolo(:,:,:)
      End Type MTHamiltonianType
      Type MTHamiltonianList
        Type (MTHamiltonianType) :: main
        Type (MTHamiltonianType) :: spinless
        Type (MTHamiltonianType) :: alpha
        Type (MTHamiltonianType) :: beta
        Type (MTHamiltonianType) :: ab
        Type (MTHamiltonianType) :: ba
        integer :: maxnlo
        integer :: maxaa
        integer,allocatable :: losize(:)
      End Type MTHamiltonianList
! Muffin-tin Hamlitonian for SCF
      Type (MTHamiltonianList) :: mt_hscf
! Relativity settings
      integer :: level_nr, level_zora, level_iora
      parameter (level_nr=0, level_zora=1, level_iora=2)
      


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
!EOP
!BOC
      Implicit None
      Type (MTHamiltonianList), Intent (Inout) :: mt_h


      nullify(mt_h%main%aa)
      nullify(mt_h%main%alo)
      nullify(mt_h%main%loa)
      nullify(mt_h%main%lolo)

      nullify(mt_h%spinless%aa)
      nullify(mt_h%spinless%alo)
      nullify(mt_h%spinless%loa)
      nullify(mt_h%spinless%lolo)

      nullify(mt_h%alpha%aa)
      nullify(mt_h%alpha%alo)
      nullify(mt_h%alpha%loa)
      nullify(mt_h%alpha%lolo)

      nullify(mt_h%beta%aa)
      nullify(mt_h%beta%alo)
      nullify(mt_h%beta%loa)
      nullify(mt_h%beta%lolo)

      nullify(mt_h%ab%aa)
      nullify(mt_h%ab%alo)
      nullify(mt_h%ab%loa)
      nullify(mt_h%ab%lolo)

      nullify(mt_h%ba%aa)
      nullify(mt_h%ba%alo)
      nullify(mt_h%ba%loa)
      nullify(mt_h%ba%lolo)

      end subroutine MTNullify

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
      Type (MTHamiltonianList) :: mt_h
      integer ::  io, ilo, if1, l, m, lm, l1, m1, lm1, l3, m3, lm3, haaijsize, is, ias

      haaijSize=0
      Do is = 1, nspecies
        if1=0
        Do l = 0, input%groundstate%lmaxmat
          Do m = - l, l
            lm = idxlm (l, m)
            Do io = 1, apword (l, is)
              if1=if1+1
            End Do
          End Do
        End Do
        if (if1.gt.haaijSize) haaijSize=if1
      Enddo
      mt_h%maxaa=haaijSize


      if (allocated(mt_h%losize)) deallocate(mt_h%losize)
      allocate(mt_h%losize(nspecies))
      mt_h%maxnlo=0
      Do is = 1, nspecies
        ias=idxas (1, is)
        ilo=nlorb (is)
        if (ilo.gt.0) then
          l1 = lorbl (ilo, is)
          lm1 = idxlm (l1, l1)
          l3 = lorbl (1, is)
          lm3 = idxlm (l3, -l3)
          mt_h%losize(is)=idxlo (lm1, ilo, ias)- idxlo (lm3, 1, ias)+1
          if (mt_h%maxnlo.lt.mt_h%losize(is)) mt_h%maxnlo=mt_h%losize(is)
        else
          mt_h%losize(is)=0
        endif
      Enddo

      call MTInit(mt_h%spinless,mt_h%maxaa,mt_h%maxnlo)
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
      subroutine MTInit(mt_block,maxaa,maxnlo)
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
     Type (MTHamiltonianType) :: mt_block
     integer :: maxaa,maxnlo

     if (.not.associated(mt_block%aa)) then
       allocate(mt_block%aa(maxaa,maxaa,natmtot))
       mt_block%aa=0d0
       if (maxnlo.gt.0) then
         allocate(mt_block%loa(maxnlo,maxaa,natmtot))
         allocate(mt_block%alo(maxaa,maxnlo,natmtot))
         allocate(mt_block%lolo(maxnlo,maxnlo,natmtot))
         mt_block%loa=0d0
         mt_block%alo=0d0
         mt_block%lolo=0d0
       endif
     else
       mt_block%aa=0d0
       if (maxnlo.gt.0) then
         mt_block%loa=0d0
         mt_block%alo=0d0
         mt_block%lolo=0d0
       endif
    
     endif

     end subroutine MTinit

!
!
!
!BOP
! !ROUTINE: MTRedirect
! !INTERFACE:
!
!
      subroutine MTRedirect(mt_blockA,mt_blockB)
! !USES:
      Use mod_atoms
!      Use modmain
! !DESCRIPTION:
! Redirects pointers of mt_blockA to mt_blockB of muffin-tin Hamiltonians or overlaps. 
!
! !REVISION HISTORY:
!   Created October 2015 (Andris)
!EOP
!BOC
     Implicit None
     Type (MTHamiltonianType) :: mt_blockA,mt_blockB
     integer :: maxaa,maxnlo

     if (associated(mt_blockB%aa)) mt_blockA%aa=>mt_blockB%aa
     if (associated(mt_blockB%alo)) mt_blockA%alo=>mt_blockB%alo
     if (associated(mt_blockB%loa)) mt_blockA%loa=>mt_blockB%loa
     if (associated(mt_blockB%lolo)) mt_blockA%lolo=>mt_blockB%lolo

     end subroutine MTRedirect
!
!
!
!BOP
! !ROUTINE: MTCopy
! !INTERFACE:
!
!
      subroutine MTCopy(mt_blockA,mt_blockB)
! !USES:
      Use mod_atoms
!      Use modmain
! !DESCRIPTION:
! Copies the contents of mt_blockA to mt_blockB of muffin-tin Hamiltonians or overlaps. 
!
! !REVISION HISTORY:
!   Created October 2015 (Andris)
!EOP
!BOC
     Implicit None
     Type (MTHamiltonianType) :: mt_blockA,mt_blockB
     integer :: maxaa,maxnlo

     if (associated(mt_blockA%aa).and.associated(mt_blockB%aa)) then
       mt_blockB%aa=mt_blockA%aa
     else
       write(*,*) 'Error (MTCopy): mt_blockA%aa or mt_blockB%aa is not allocated'
       stop
     endif
     if (associated(mt_blockA%lolo)) then
       if (associated(mt_blockA%alo).and.associated(mt_blockB%alo)) then 
         mt_blockB%alo=mt_blockA%alo
       else
         write(*,*) 'Error (MTCopy): mt_blockA%alo or mt_blockB%alo is not allocated'
         stop
       endif
       if (associated(mt_blockB%loa).and.associated(mt_blockB%loa)) then 
         mt_blockB%loa=mt_blockA%loa
       else
         write(*,*) 'Error (MTCopy): mt_blockA%loa or mt_blockB%loa is not allocated'
         stop
       endif
       if (associated(mt_blockB%lolo).and.associated(mt_blockA%lolo)) then 
         mt_blockB%lolo=mt_blockA%lolo
       else
         write(*,*) 'Error (MTCopy): mt_blockA%lolo or mt_blockB%lolo is not allocated'
         stop
       endif
     endif
     end subroutine MTCopy


!
!
!
!BOP
! !ROUTINE: MTRelease
! !INTERFACE:
!
!
     subroutine MTRelease(mt_h)
! !USES:
! !DESCRIPTION:
! Release all memory used by the muffin-tin Hamiltonian or overlap. 
!
! !REVISION HISTORY:
!   Created June 2019 (Andris)
!EOP
!BOC
     implicit none
     type (MTHamiltonianList), Intent (Inout) :: mt_h

     if (associated(mt_h%spinless%aa)) deallocate(mt_h%spinless%aa)
     if (associated(mt_h%spinless%alo)) deallocate(mt_h%spinless%alo)
     if (associated(mt_h%spinless%loa)) deallocate(mt_h%spinless%loa)
     if (associated(mt_h%spinless%lolo)) deallocate(mt_h%spinless%lolo)

     if (associated(mt_h%alpha%aa)) deallocate(mt_h%alpha%aa)
     if (associated(mt_h%alpha%alo)) deallocate(mt_h%alpha%alo)
     if (associated(mt_h%alpha%loa)) deallocate(mt_h%alpha%loa)
     if (associated(mt_h%alpha%lolo)) deallocate(mt_h%alpha%lolo)

     if (associated(mt_h%beta%aa)) deallocate(mt_h%beta%aa)
     if (associated(mt_h%beta%alo)) deallocate(mt_h%beta%alo)
     if (associated(mt_h%beta%loa)) deallocate(mt_h%beta%loa)
     if (associated(mt_h%beta%lolo)) deallocate(mt_h%beta%lolo)

     if (associated(mt_h%ab%aa)) deallocate(mt_h%ab%aa)
     if (associated(mt_h%ab%alo)) deallocate(mt_h%ab%alo)
     if (associated(mt_h%ab%loa)) deallocate(mt_h%ab%loa)
     if (associated(mt_h%ab%lolo)) deallocate(mt_h%ab%lolo)

     if (associated(mt_h%ba%aa)) deallocate(mt_h%ba%aa)
     if (associated(mt_h%ba%alo)) deallocate(mt_h%ba%alo)
     if (associated(mt_h%ba%loa)) deallocate(mt_h%ba%loa)
     if (associated(mt_h%ba%lolo)) deallocate(mt_h%ba%lolo)

!     if (associated(mt_h%main%aa)) deallocate(mt_h%main%aa)
!     if (associated(mt_h%main%alo)) deallocate(mt_h%main%alo)
!     if (associated(mt_h%main%loa)) deallocate(mt_h%main%loa)
!     if (associated(mt_h%main%lolo)) deallocate(mt_h%main%lolo)

     if (allocated(mt_h%losize)) deallocate(mt_h%losize)

     end subroutine MTRelease

End Module
!
