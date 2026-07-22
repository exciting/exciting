!> Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma,
!> C. Meisenbichler and C. Ambrosch-Draxl.
!> This file is distributed under the terms of the
!> GNU General Public License.
!> See the file COPYING for license details.

!>  APW and local-orbital variables 
module mod_apw_lo

      use constants, only : maxspecies, maxlapw
      use precision, only : i32, dp 
   
      implicit none

      private

      public :: mtbasisinit, mtbasisrelease, apw_lo_basis_type, &
                load_apwlo, save_apwlo, maxlapw
   
      !> Maximum allowable APW order
      integer(i32), public, parameter :: maxapword = 4
   
      !> APW order
      integer(i32), public :: apword(0:maxlapw, maxspecies)
   
      !> Maximum APW order over all angular momenta and species
      integer(i32), public :: apwordmax
   
      !> APW initial linearisation energies
      real(dp), public :: apwe0(maxapword, 0:maxlapw, maxspecies)
   
      !> APW linearisation energies
      real(dp), allocatable, public :: apwe(:, :, :)
   
      !> APW derivative order
      integer(i32), public :: apwdm(maxapword, 0:maxlapw, maxspecies)
   
      !> APW principal quantum number
      !
      !> If a principal quantum number is specified for a custom (L)APW
      !> in the species file, it is used to calculate the corresponding
      !> trial energy automatically.
      !
      !> If it is not specified, the principal quantum number is set to -1.
      integer(i32), public :: apwn(maxapword, 0:maxlapw, maxspecies)
   
      !> Default (L)APW principal quantum number
      integer(i32), parameter, public :: default_apwn = -1
   
      !> True if the linearisation energies are allowed to vary
      logical(i32), public :: apwve(maxapword, 0:maxlapw, maxspecies)
   
      !> APW radial functions
      real(dp), target, allocatable, public :: apwfr(:, :, :, :, :)
   
      !> Derivative of radial functions at the muffin-tin surface
      !> Note (mrm): This is currently unused
      real(dp), allocatable, public :: apwdfr(:, :, :)
   
      !> Maximum number of local orbitals
      integer(i32), parameter, public :: maxlorb = 100
   
      !> Maximum allowable local-orbital order
      integer(i32), parameter, public :: maxlorbord = 4
   
      !> Number of local orbitals
      integer(i32), public :: nlorb(maxspecies)
   
      !> Maximum nlorb over all species
      integer(i32), public :: nlomax
   
      !> Total number of local orbitals
      integer(i32), public :: nlotot
   
      !> Local-orbital order
      integer(i32), public :: lorbord(maxlorb, maxspecies)
   
      !> Local-orbital angular momentum
      integer(i32), public :: lorbl(maxlorb, maxspecies)
   
      !> Local-orbital principal quantum number
      integer(i32), public :: lorbn(maxlorbord, maxlorb, maxspecies)
   
      !> Default local-orbital principal quantum number
      integer(i32), parameter, public :: default_lorbn = -1
   
      !> Local-orbital relativistic quantum number:
      !> kappa = (l - j)(2j + 1)
      integer(i32), public :: lorbk(maxlorb, maxspecies)
   
      !> Wave-function relativistic quantum number:
      !> kappa = (l - j)(2j + 1)
      integer(i32), public :: wfkappa(maxlorbord, maxlorb, maxspecies)
   
      !> Maximum lorbl over all species
      integer(i32), public :: lolmax
   
      !> (lolmax + 1)^2
      integer(i32), public :: lolmmax
   
      !> Local-orbital initial energies
      real(dp), public :: lorbe0(maxlorbord, maxlorb, maxspecies)
   
      !> Local-orbital energies
      real(dp), allocatable, public :: lorbe(:, :, :)
   
      !> Local-orbital derivative order
      integer(i32), public :: lorbdm(maxlorbord, maxlorb, maxspecies)
   
      !> True if the linearisation energies are allowed to vary
      logical(i32), public :: lorbve(maxlorbord, maxlorb, maxspecies)
   
      !> True if the local orbital is used as a Wannier projector
      !> for band-structure interpolation
      logical(i32), public :: lorbwfproj(maxlorb, maxspecies)
   
      !> Local-orbital radial functions
      real(dp), target, allocatable, public :: lofr(:, :, :, :)
   
      !> Minimum default linearisation energy over all APW and
      !> local-orbital functions
      real(dp), public :: mine0
   

      !> Muffin-tin basis container   
      type :: apw_lo_basis_type
         real(dp), pointer :: apwfr(:, :, :, :, :)
         real(dp), pointer :: lofr(:, :, :, :)
      end type apw_lo_basis_type


      interface
            !> Save APW and local-orbital data 
            module subroutine save_apwlo()
            end subroutine save_apwlo
            !> Load APW and local-orbital data
            module subroutine load_apwlo()
            end subroutine load_apwlo
      end interface
   
   contains
   
      !--------------------------------------------------------------------!
      !> Initialise muffin-tin basis storage                                !
      !--------------------------------------------------------------------!
   
      subroutine mtbasisinit(mt_basis)
   
         use modinput
         use mod_muffin_tin
         use mod_atoms
   
         implicit none
   
         type(apw_lo_basis_type) :: mt_basis
   
         nullify(mt_basis%apwfr)
   
         allocate(mt_basis%apwfr( &
            nrmtmax,                         &
            2,                               &
            apwordmax,                       &
            0:input%groundstate%lmaxapw,     &
            natmtot))
   
         nullify(mt_basis%lofr)
   
         allocate(mt_basis%lofr( &
            nrmtmax, &
            2,       &
            nlomax,  &
            natmtot))
   
      end subroutine mtbasisinit
   
      !--------------------------------------------------------------------!
      !> Release muffin-tin basis storage                                   !
      !--------------------------------------------------------------------!
   
      subroutine mtbasisrelease(mt_basis)
   
         implicit none
   
         type(apw_lo_basis_type) :: mt_basis
   
         deallocate(mt_basis%apwfr)
         deallocate(mt_basis%lofr)
   
         nullify(mt_basis%apwfr)
         nullify(mt_basis%lofr)
   
      end subroutine mtbasisrelease
   
end module mod_apw_lo
