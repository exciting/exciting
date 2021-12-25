
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!> Potential and density variables 
Module mod_potential_and_density
      Use modinput
      implicit none 
! exchange-correlation functional type
 	integer::xctype(3)
! exchange-correlation functional description
      Character (512) :: xcdescr
! exchange-correlation functional spin treatment
      Integer :: xcspin
! exchange-correlation functional density gradient treatment
      Integer :: xcgrad
! exchange mixing parameter for hybrid functionals
      Real (8) :: ex_coef
      Real (8) :: ec_coef
! ?
      Integer :: maxncv
! muffin-tin charge density
      Real (8), Allocatable :: rhomt (:, :, :)
! interstitial real-space charge density
      Real (8), Allocatable :: rhoir (:)
! muffin-tin magnetisation vector field
      Real (8), Allocatable :: magmt (:, :, :, :)
! interstitial magnetisation vector field
      Real (8), Allocatable :: magir (:, :)
! muffin-tin Coulomb potential
      Real (8), Allocatable :: vclmt (:, :, :)
! Madulung potential (excluding the on-ste nuclear contribution) in the innermost radial point
      Real (8), Allocatable :: vmad (:)
! interstitial real-space Coulomb potential
      Real (8), Allocatable :: vclir (:)
! order of polynomial for pseudocharge density
!replaced by inputstructureinteger::npsden
! muffin-tin exchange-correlation potential
      Real (8), Allocatable :: vxcmt (:, :, :)
      !real (8), allocatable :: vxmt(:,:,:)
      !real (8), allocatable :: vcmt(:,:,:)
! interstitial real-space exchange-correlation potential
      Real (8), Allocatable :: vxcir (:)
! muffin-tin exchange-correlation magnetic field
      Real (8), Allocatable :: bxcmt (:, :, :, :)
! interstitial exchange-correlation magnetic field
      Real (8), Allocatable :: bxcir (:, :)
! nosource is .true. if the field is to be made source-free
!replaced by inputstructurelogical::nosource
! muffin-tin effective potential
      Real (8), Allocatable :: veffmt (:, :, :)
! interstitial effective potential
      Real (8), Allocatable :: veffir (:)

! muffin-tin 'vhalf' potential (needed for DFT-1/2 calculations)
      Real (8), Allocatable :: vhalfmt (:, :, :)
! interstitial 'vhalf' potential (needed for DFT-1/2 calculations)
      Real (8), Allocatable :: vhalfir (:)
! intermediate 'vhalf' potential (needed for DFT-1/2 calculations)
      Real (8), Allocatable :: vhalfsph (:,:)
! partial densities (used in DFT-1/2 calculations)
	  Real (8), Allocatable :: nalphamt(:,:,:,:,:), nalphair(:,:,:)

! G-space interstitial effective potential
      Complex (8), Allocatable :: veffig (:),meffig (:), m2effig (:)
! interstitial effective potential at G=0
      real(8) :: veffig0
      real(8) :: meffig0
      real(8) :: m2effig0

! muffin-tin exchange energy density
      Real (8), Allocatable :: exmt (:, :, :)
! interstitial real-space exchange energy density
      Real (8), Allocatable :: exir (:)
! muffin-tin correlation energy density
      Real (8), Allocatable :: ecmt (:, :, :)
! interstitial real-space correlation energy density
      Real (8), Allocatable :: ecir (:)
! type of mixing to use for the potential
!replaced by inputstructureinteger::mixtype
! adaptive mixing parameters
!replaced by inputstructurereal(8)::beta0
!replaced by inputstructurereal(8)::betainc
!replaced by inputstructurereal(8)::betadec

! type for the density matrix in the muffin-tin region
! It is similar to the Hamiltonian type. Do we need a separate type here?
      Type MTDensityMatrixType
        complex(8), pointer :: ff(:,:,:)
      End Type MTDensityMatrixType
      Type MTDensityMatrixList
        Type (MTDensityMatrixType) :: main
        Type (MTDensityMatrixType) :: alpha
        Type (MTDensityMatrixType) :: beta
        Type (MTDensityMatrixType) :: ab
        Type (MTDensityMatrixType) :: ba
        integer :: maxnlo
        integer :: maxaa
        integer, pointer :: losize(:)
      End Type MTDensityMatrixList

      Type(MTDensityMatrixList) :: mt_dm

Contains
!
!
!
!BOP
! !ROUTINE: DMNullify
! !INTERFACE:
!
!
      subroutine DMNullify(mt_dm)
! !USES:
! !DESCRIPTION:
! Nullifies all pointers of the muffin-tin density matrix. 
!
! !REVISION HISTORY:
!   Created June 2019 (Andris)
!EOP
!BOC
      Implicit None
      Type (MTDensityMatrixList), Intent (Inout) :: mt_dm


      nullify(mt_dm%main%ff)
      nullify(mt_dm%alpha%ff)
      nullify(mt_dm%beta%ff)
      nullify(mt_dm%ab%ff)
      nullify(mt_dm%ba%ff)

      end subroutine DMNullify
      

!
!
!
!BOP
! !ROUTINE: DMInit
! !INTERFACE:
!
!
      subroutine DMInit(mt_block,maxaa,maxnlo)
! !USES:
      Use mod_atoms
! !DESCRIPTION:
! Initialises a part of the muffin-tin density matrix. 
!
! !REVISION HISTORY:
!   Created July 2019 (Andris)
!EOP
!BOC
     Implicit None
     Type (MTDensityMatrixType) :: mt_block
     integer :: maxaa,maxnlo

     if (.not.associated(mt_block%ff)) allocate(mt_block%ff(maxaa+maxnlo,maxaa+maxnlo,natmtot))
     mt_block%ff=0d0

     end subroutine DMinit


!
!
!
!BOP
! !ROUTINE: DMInitAll
! !INTERFACE:
!
!
      subroutine DMInitAll(mt_dm)
! !USES:
      Use modinput
      Use mod_APW_LO
      Use mod_atoms
      Use mod_muffin_tin
      Use mod_eigensystem
      Use mod_spin, only: ncmag
! !DESCRIPTION:
! Initialises all parts of the muffin-tin density matrix. 
!
! !REVISION HISTORY:
!   Created July 2019 (Andris)
!EOP
!BOC
      Implicit None
      Type (MTDensityMatrixList), Intent (Inout) :: mt_dm
      integer ::  io, ilo, if1, l, m, lm, l1, m1, lm1, l3, m3, lm3, maxaa, is, ias

      maxaa=0
      Do is = 1, nspecies
        if1=0
        Do l = 0, input%groundstate%lmaxapw
          Do m = - l, l
            lm = idxlm (l, m)
            Do io = 1, apword (l, is)
              if1=if1+1
            End Do
          End Do
        End Do
        if (if1.gt.maxaa) maxaa=if1
      Enddo
      mt_dm%maxaa=maxaa


      if (associated(mt_dm%losize)) deallocate(mt_dm%losize)
      allocate(mt_dm%losize(nspecies))
      mt_dm%maxnlo=0
      Do is = 1, nspecies
        ias=idxas (1, is)
        ilo=nlorb (is)
        if (ilo.gt.0) then
          l1 = lorbl (ilo, is)
          lm1 = idxlm (l1, l1)
          l3 = lorbl (1, is)
          lm3 = idxlm (l3, -l3)
          mt_dm%losize(is)=idxlo (lm1, ilo, ias)- idxlo (lm3, 1, ias)+1
          if (mt_dm%maxnlo.lt.mt_dm%losize(is)) mt_dm%maxnlo=mt_dm%losize(is)
        else
          mt_dm%losize(is)=0
        endif
      Enddo

      call DMInit(mt_dm%alpha,mt_dm%maxaa,mt_dm%maxnlo)
      if (associated(input%groundstate%spin)) then
        call DMInit(mt_dm%beta,mt_dm%maxaa,mt_dm%maxnlo)
        call DMInit(mt_dm%ba,mt_dm%maxaa,mt_dm%maxnlo)
        if (ncmag) then
          call DMInit(mt_dm%ab,mt_dm%maxaa,mt_dm%maxnlo)
        endif
      endif
      
      end subroutine DMInitAll


!
!
!
!BOP
! !ROUTINE: DMRelease
! !INTERFACE:
!
!
     subroutine DMRelease(mt_dm)
! !USES:
! !DESCRIPTION:
! Release all memory used by the density matrix. 
!
! !REVISION HISTORY:
!   Created July 2019 (Andris)
!EOP
!BOC
     implicit none
     Type (MTDensityMatrixList), Intent (Inout) :: mt_dm

!     write(*,*) "DM_Release " ,associated(mt_dm%alpha%ff)
!     read(*,*)
     if (associated(mt_dm%alpha%ff)) deallocate(mt_dm%alpha%ff)
     if (associated(mt_dm%beta%ff)) deallocate(mt_dm%beta%ff)
     if (associated(mt_dm%ab%ff)) deallocate(mt_dm%ab%ff)
     if (associated(mt_dm%ba%ff)) deallocate(mt_dm%ba%ff)
     if (associated(mt_dm%losize)) deallocate(mt_dm%losize)

     end subroutine DMRelease
      
      
End Module
!
