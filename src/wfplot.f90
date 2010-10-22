!
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine wfplot (dostm)
      Use modmain
      Use modinput
      Implicit None
      Logical, Intent (In) :: dostm
! local variables
      Integer :: ik, ist
      Real (8) :: x, t1
! allocatable arrays
      Complex (8), Allocatable :: evecfv (:, :)
      Complex (8), Allocatable :: evecsv (:, :)
! external functions
      Real (8) :: sdelta
      External sdelta
! initialise universal variables
      Call init0
      Call init1
      Allocate (evecfv(nmatmax, nstfv))
      Allocate (evecsv(nstsv, nstsv))
! read the density and potentials from file
      Call readstate
! read Fermi energy from file
      Call readfermi
! find the new linearisation energies
      Call linengy
! generate the APW radial functions
      Call genapwfr
! generate the local-orbital radial functions
      Call genlofr
! set the occupancies
      If ( .Not. dostm) Then
         ik = kstlist (1, 1)
         ist = kstlist (2, 1)
         If ((ik .Lt. 1) .Or. (ik .Gt. nkpt)) Then
            Write (*,*)
            Write (*, '("Error(wfplot): k-point out of range : ", I8)') &
           & ik
            Write (*,*)
            Stop
         End If
         If ((ist .Lt. 1) .Or. (ist .Gt. nstsv)) Then
            Write (*,*)
            Write (*, '("Error(wfplot): state out of range : ", I8)') &
           & ist
            Write (*,*)
            Stop
         End If
! plotting a single wavefunction
         occsv (:, :) = 0.d0
         occsv (ist, ik) = 1.d0
      Else
! plotting an STM image by setting occupancies to be a delta function at the
! Fermi energy
         t1 = 1.d0 / input%groundstate%swidth
         Do ik = 1, nkpt
! get the eigenvalues from file
            Call getevalsv (vkl(:, ik), evalsv(:, ik))
            Do ist = 1, nstsv
               x = (efermi-evalsv(ist, ik)) * t1
               occsv (ist, ik) = occmax * wkpt (ik) * sdelta &
              & (input%groundstate%stypenumber, x) * t1
            End Do
         End Do
      End If
! set the charge density to zero
      rhomt (:, :, :) = 0.d0
      rhoir (:) = 0.d0
! compute the charge density with the new occupancies
      Do ik = 1, nkpt
! get the eigenvectors from file
         Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
         Call getevecsv (vkl(:, ik), evecsv)
         Call rhovalk (ik, evecfv, evecsv)
      End Do
! symmetrise the density for the STM plot
      If (dostm) Then
         Call symrf (input%groundstate%lradstep, rhomt, rhoir)
      End If
! convert the density from a coarse to a fine radial mesh
      Call rfmtctof (rhomt)
! write the wavefunction modulus squared plot to file
      If (associated(input%properties%wfplot%plot1d)) Then
         Call plot1d ("WF", 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rhomt, rhoir, input%properties%wfplot%plot1d)
         Write (*,*)
         Write (*, '("Info(wfplot):")')
         Write (*, '(" 1D wavefunction modulus squared written to WF1D.&
        &OUT")')
         Write (*, '(" vertex location lines written to WFLINES.OUT")')
      End If
      If (associated(input%properties%wfplot%plot2d)) Then
         Call plot2d ("WF", 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rhomt, rhoir, input%properties%wfplot%plot2d)
         Write (*,*)
         Write (*, '("Info(wfplot):")')
         Write (*, '(" 2D wavefunction modulus squared written to WF2D.&
        &OUT")')
      End If
      If (dostm) Then
         Call plot2d ("STM", 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rhomt, rhoir, input%properties%STM%plot2d)
         Write (*,*)
         Write (*, '("Info(wfplot):")')
         Write (*, '(" 2D STM image written to STM2d.xml")')
      End If
      If (associated(input%properties%wfplot%plot3d)) Then
         Call plot3d ("WF", 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rhomt, rhoir, input%properties%wfplot%plot3d)
         Write (*,*)
         Write (*, '("Info(wfplot):")')
         Write (*, '(" 3D wavefunction modulus squared written to WF3D.&
        &OUT")')
      End If
      If ( .Not. dostm) Then
         Write (*, '(" for k-point ", I6, " and state ", I6)') kstlist &
        & (1, 1), kstlist (2, 1)
      End If
      Write (*,*)
      Deallocate (evecfv, evecsv)
      Return
End Subroutine
