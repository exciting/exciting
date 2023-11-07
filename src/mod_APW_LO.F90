
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
Module mod_APW_LO
      use constants, only: maxspecies, maxlapw
      implicit none

!-----------------------------------------!
!     APW and local-orbital variables     !
!-----------------------------------------!
! maximum allowable APW order
      Integer, Parameter :: maxapword = 4
! APW order
      Integer :: apword (0:maxlapw, maxspecies)
! maximum of apword over all angular momenta and species
      Integer :: apwordmax
! APW initial linearisation energies
      Real (8) :: apwe0 (maxapword, 0:maxlapw, maxspecies)
! APW linearisation energies
      Real (8), Allocatable :: apwe (:, :, :)
! APW derivative order
      Integer :: apwdm (maxapword, 0:maxlapw, maxspecies)
! APW principal quantum number
! If a principal quantum number is speciefied for a custom (L)APW in the species file,
! it is used to calculate the according trial energy automatically.
! If it is not specified, the principal quantum number is set to -1 by default.
      Integer :: apwn (maxapword, 0:maxlapw, maxspecies)
! Default (L)APW principal quantum number
      Integer, parameter :: default_apwn = -1
! apwve is .true. if the linearisation energies are allowed to vary
      Logical :: apwve (maxapword, 0:maxlapw, maxspecies)
! APW radial functions
      Real (8), target, Allocatable :: apwfr (:, :, :, :, :)
! derivate of radial functions at the muffin-tin surface
      Real (8), Allocatable :: apwdfr (:, :, :)
! maximum number of local-orbitals
      Integer, Parameter :: maxlorb = 100
! maximum allowable local-orbital order
      Integer, Parameter :: maxlorbord = 4
! number of local-orbitals
      Integer :: nlorb (maxspecies)
! maximum nlorb over all species
      Integer :: nlomax
! total number of local-orbitals
      Integer :: nlotot
! local-orbital order
      Integer :: lorbord (maxlorb, maxspecies)
! local-orbital angular momentum
      Integer :: lorbl (maxlorb, maxspecies)
! local-orbital principle quantum number                                                               
      Integer :: lorbn (maxlorbord, maxlorb, maxspecies)
! Default local orbital principal quantum number
      Integer, parameter :: default_lorbn = -1
! local-orbital relativistic quantum number kappa = (l-j)(2j+1)                                                                                         
      Integer :: lorbk (maxlorb,maxspecies)
! wave-function relativistic quantum number kappa = (l-j)(2j+1)     
      Integer :: wfkappa (maxlorbord, maxlorb, maxspecies)
! maximum lorbl over all species
      Integer :: lolmax
! (lolmax+1)^2
      Integer :: lolmmax
! local-orbital initial energies
      Real (8) :: lorbe0 (maxlorbord, maxlorb, maxspecies)
! local-orbital energies
      Real (8), Allocatable :: lorbe (:, :, :)
! local-orbital derivative order
      Integer :: lorbdm (maxlorbord, maxlorb, maxspecies)
! lorbve is .true. if the linearisation energies are allowed to vary
      Logical :: lorbve (maxlorbord, maxlorb, maxspecies)
! lorbwfproj is .true. if the local-orbital is used as a Wannier-projector for bandstructure interpolation
      Logical :: lorbwfproj (maxlorb, maxspecies)
! local-orbital radial functions
      Real (8), target, Allocatable :: lofr (:, :, :, :)
! energy step size for locating the band energy
!replaced by inputstructurereal(8)::deband
! minimum of the default linearisation energy over all APW and local-orbitals
! functions
      real(8) :: mine0
! Data structure for storing muffin-tin basis functions
      Type apw_lo_basis_type
        Real (8), pointer :: apwfr (:, :, :, :, :),lofr (:, :, :, :)
      end type
!      type (apw_lo_basis_type) :: mt_basis !,mt_basis_alpha,mt_basis_beta

Contains

!
!
!
!BOP
! !ROUTINE: MTBasisInit
! !INTERFACE:
!
!
      subroutine MTBasisInit(mt_basis)
! !USES:
      Use modinput
      Use mod_muffin_tin
      Use mod_atoms
! !DESCRIPTION:
! Initialises storage for basis functions in the muffin-tin region. 
!
! !REVISION HISTORY:
!   Created June 2019 (Andris)
!EOP
!BOC
      implicit none 
      type (apw_lo_basis_type) :: mt_basis

      nullify(mt_basis%apwfr)
      allocate (mt_basis%apwfr(nrmtmax, 2, apwordmax, 0:input%groundstate%lmaxapw, natmtot))
      nullify(mt_basis%lofr)
      allocate (mt_basis%lofr(nrmtmax, 2, nlomax, natmtot))
      
      end subroutine MTBasisInit

!
!
!
!BOP
! !ROUTINE: MTBasisRelease
! !INTERFACE:
!
!
      subroutine MTBasisRelease(mt_basis)
! !USES:
! !DESCRIPTION:
! Releases storage for basis functions in the muffin-tin region. 
!
! !REVISION HISTORY:
!   Created June 2019 (Andris)
!EOP
!BOC
      implicit none
      type (apw_lo_basis_type) :: mt_basis

      deallocate (mt_basis%apwfr)
      deallocate (mt_basis%lofr)
      nullify(mt_basis%apwfr)
      nullify(mt_basis%lofr)

      end subroutine MTBasisRelease


End Module
!
