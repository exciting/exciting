!
!
!
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine rdmminn
! minimise the total energy w.r.t. occupation numbers
      Use modinput
      Use modmain
      Implicit None
! allocatable arrays
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Integer :: ik, it, idm
      Real (8) :: ep, de
! parameter to check energy convergence
      Real (8), Parameter :: eps = 1.d-8
! allocate arrays
      Allocate (evecfv(nmatmax, nstfv))
      Allocate (evecsv(nstsv, nstsv))
      Open (61, File='RDMN_ENERGY.OUT', Action='WRITE', Form='FORMATTED&
     &')
      If (associated(input%groundstate%spin)) Then
         Open (62, File='RDMN_MOMENT.OUT', Action='WRITE', Form='FORMAT&
        &TED')
      End If
! calculate the non-local matrix elements (i-jj-i)
      If ((input%groundstate%RDMFT%rdmxctype .Ne. 0) .And. &
     & (input%groundstate%RDMFT%maxitc .Lt. 1)) Then
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
         Do ik = 1, nkpt
!$OMP CRITICAL
            Write (*, '("Info(rdmminn): ", I6, " of ", I6, " k-points")&
           &') ik, nkpt
!$OMP END CRITICAL
            Call rdmvnln (ik)
         End Do
!$OMP END DO
!$OMP END PARALLEL
      End If
      ep = 0.d0
! begin iteration loop
      Do it = 1, input%groundstate%RDMFT%maxitn
         Write (*, '("Info(rdmminn): iteration ", I4, " of ", I4)') it, &
        & input%groundstate%RDMFT%maxitn
! vary the occupation numbers
         Call rdmvaryn
! zero the density
         rhomt (:, :, :) = 0.d0
         rhoir (:) = 0.d0
! zero the magnetisation
         If (associated(input%groundstate%spin)) Then
            magmt (:, :, :, :) = 0.d0
            magir (:, :) = 0.d0
         End If
! compute the charge density and magnetisation with the new occupancies
         Do ik = 1, nkpt
! get the eigenvectors from file
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
! calculate Coulomb potential matrix elements (RDM states)
         Call genvmat (vclmt, vclir, vclmat)
! calculate the energy
         Call rdmenergy
! check for convergence
         de = ep - engytot
         If (it .Gt. 1) Then
            If (Abs(de) .Lt. eps) Go To 10
         End If
         ep = engytot
! write energy and convergence factor to a file
         Write (61, '(I6, 2G18.10)') it, engytot, de
         Call flushifc (61)
! end iteration loop
      End Do
10    Continue
      Close (61)
      If (associated(input%groundstate%spin)) close (62)
      Deallocate (evecfv, evecsv)
      Return
End Subroutine
