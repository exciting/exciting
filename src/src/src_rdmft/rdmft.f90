!
!
!
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine rdmft
! 1-reduced density matrix functional theory
      Use modinput
      Use modmain
      Implicit None
! local variables
      Integer :: ik
      Call init0
      Call init1
! generate q-point set and wiq2 array
      Call init2
! read density and potentials from file
      Call readstate
! generate the core wavefunctions and densities
      Call gencore
! find the new linearisation energies
      Call linengy
! generate the APW radial functions
      Call genapwfr
! generate the local-orbital radial functions
      Call genlofr
! compute the overlap radial integrals
      Call olprad
! compute the Hamiltonian radial integrals
      Call hmlrad
! compute the kinetic energy of the core
      Call energykncr
! generate the kinetic matrix elements
      Call genkinmat
! read in the occupancies
      Do ik = 1, nkpt
         Call getoccsv (vkl(:, ik), occsv(:, ik))
      End Do
! calculate Coulomb potential matrix elements
      Call genvmat (vclmt, vclir, vclmat)
! derivative of kinetic energy w.r.t. evecsv
      Call rdmdkdc
! open information files
      Open (60, File='RDM_INFO.OUT', Action='WRITE', Form='FORMATTED')
! write out general information to RDM_INFO.OUT
      Call writeinfo (60)
! begin main self-consistent loop
      Do iscl = 1, input%groundstate%rdmft%rdmmaxscl
         Write (60,*)
         Write (60, '("+-------------------------+")')
         Write (60, '("| Iteration number : ", I4, " |")') iscl
         Write (60, '("+-------------------------+")')
         Call flushifc (60)
! minimisation over natural orbitals
         If (input%groundstate%rdmft%maxitc .Ge. 1) Then
            Call rdmminc
            Write (60,*)
            Write (60, '("Natural orbital minimisation done")')
            Call rdmwriteengy (60)
         End If
! minimisation over occupation number
         If (input%groundstate%rdmft%maxitn .Ge. 1) Then
            Call rdmminn
            Write (60,*)
            Write (60, '("Occupation number minimisation done")')
            Call rdmwriteengy (60)
         End If
! end loop over iscl
      End Do
! write density to STATE.OUT
      Call writestate
! write occupation numbers for restart
      Do ik = 1, nkpt
         Call putoccsv (ik, occsv(:, ik))
      End Do
! close RDM_INFO.OUT file
      Close (60)
      Return
End Subroutine
