!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine rdmminc
! minimises the total energy w.r.t. evecsv using steepest descent
      Use modinput
      Use modmain
      Implicit None
      Integer :: it, ik, idm
      Real (8) :: sum, sp, ds
! parameter to check energy convergence
      Real (8), Parameter :: eps = 1.d-10
! allocatable arrays
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
! allocate arrays
      Allocate (evecfv(nmatmax, nstfv))
      Allocate (evecsv(nstsv, nstsv))
      sp = 0.d0
      ds = 0.d0
      Open (61, File='RDMC_ENERGY.OUT', Action='WRITE', Form='FORMATTED&
     &')
      If (associated(input%groundstate%spin)) Then
         Open (62, File='RDMC_MOMENT.OUT', Action='WRITE', Form='FORMAT&
        &TED')
      End If
      Write (*,*)
! begin iteration loop
      Do it = 1, input%groundstate%RDMFT%maxitc
         Write (*, '("Info(rdmminc): iteration ", I4, " of ", I4)') it, &
        & input%groundstate%RDMFT%maxitc
! vary evecsv and orthogonalise it
         Call rdmvaryc (sum)
! zero the density
         rhomt (:, :, :) = 0.d0
         rhoir (:) = 0.d0
! zero the magnetisation
         If (associated(input%groundstate%spin)) Then
            magmt (:, :, :, :) = 0.d0
            magir (:, :) = 0.d0
         End If
         Do ik = 1, nkpt
! get the eigenvectors and values from file
            Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
            Call getevecsv (vkl(:, ik), evecsv)
! calculate the density
            Call rhovalk (ik, evecfv, evecsv)
         End Do
! symmetrise the density
         Call symrf (input%groundstate%lradstep, rhomt, rhoir)
! convert the muffin-tin density from coarse to a fine grid
         Call rfmtctof (rhomt)
         If (associated(input%groundstate%spin)) Then
! symmetrise the magnetisation
            Call symrvf (input%groundstate%lradstep, magmt, magir)
! convert the magnetisation from a coarse to a fine radial mesh
            Do idm = 1, ndmag
               Call rfmtctof (magmt(:, :, :, idm))
            End Do
         End If
! add core density to the valence density
         Call addrhocr
! calculate the charges
         Call charge
! calculate the magnetic moment
         If (associated(input%groundstate%spin)) Then
            Call moment
            Write (62, '(I6, 3G18.10)') it, momtot (1:ndmag)
            Call flushifc (62)
         End If
! normalise the density
         Call rhonorm
! calculate the Coulomb potential
         Call potcoul
! calculate Coulomb matrix elements
         Call genvmat (vclmt, vclir, vclmat)
! calculate derivative of kinetic energy w.r.t. evecsv
         Call rdmdkdc
! calculate the energy
         Call rdmenergy
! check for convergence of derivative of energy w.r.t. evecsv
         If (it .Gt. 1) Then
            ds = sp - sum
            sp = sum
         End If
! write energy and convergence factor to a file
         Write (61, '(I6, 2G18.10)') it, engytot, ds
         Call flushifc (61)
! end iteration loop
      End Do
      Close (61)
      Close (62)
      Deallocate (evecfv, evecsv)
      Return
End Subroutine
