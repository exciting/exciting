!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
!BOP
! !ROUTINE: hartfock
! !INTERFACE:
Subroutine hartfock
! !USES:
      Use modmain
      Use modinput
! !DESCRIPTION:
!  Computes the self-consistent Hartree Fock ground state.
!
! !REVISION HISTORY:
!  ... 
!EOP
!BOC
      Implicit None
! local variables
      Logical :: exist
      Integer :: ik, is, ia, idm, Recl
      Real(8) :: etp, de
      character*(77) :: string
! allocatable arrays
      Complex(8), Allocatable :: evecfv (:, :, :)
      Complex(8), Allocatable :: evecsv (:, :)
! initialise universal variables
      Call init0
      Call init1
      Call init2
! open INFO.OUT file
      Open (60, File='INFO'//trim(filext), Action='WRITE', Form='FORMAT&
     &TED')
! open TOTENERGY.OUT
      Open (61, File='TOTENERGY'//trim(filext), Action='WRITE', Form='F&
     &ORMATTED')
! open FERMIDOS.OUT
      Open (62, File='FERMIDOS'//trim(filext), Action='WRITE', Form='FO&
     &RMATTED')
! open MOMENT.OUT if required
      If (associated(input%groundstate%spin)) open (63, file='MOMENT'//&
     & trim(filext), action='WRITE', form='FORMATTED')
! open FORCEMAX.OUT if required
      If (input%groundstate%tforce) open (64, file='FORCEMAX'//&
     & trim(filext), action='WRITE', form='FORMATTED')
! open DENERGY.OUT
      Open (65, File='DENERGY'//trim(filext), Action='WRITE', Form='FOR&
     &MATTED')
! write out general information to INFO.OUT
      Call writeinfo (60)
! read the charge density and potentials from file
      Call readstate
! compute the effective potential
      Call poteff
! Fourier transform effective potential to G-space
      Call genveffig
      Call genmeffig
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
      Call hmlint
! compute the Hamiltonian radial integrals
      Call hmlrad
! generate the kinetic matrix elements
      Call genkinmat
! find the occupation numbers and Fermi energy
      Call occupy
10    Continue
! set last iteration flag
      tlast = .False.
      etp = 0.d0
! begin the self-consistent loop
      Write (60,*)
      Write (60, '("+------------------------------+")')
      Write (60, '("| Self-consistent loop started |")')
      Write (60, '("+------------------------------+")')
      Do iscl = 1, input%groundstate%maxscl
         Write (60,*)
         Write (60, '("+-------------------------+")')
         Write (60, '("| Iteration number : ", I4, " |")') iscl
         Write (60, '("+-------------------------+")')
         Call flushifc (60)
         If (iscl .Ge. input%groundstate%maxscl) Then
            Write (60,*)
            Write (60, '("Reached self-consistent loops maximum")')
            tlast = .True.
         End If
!x$OMP PARALLEL DEFAULT(SHARED) &
!x$OMP PRIVATE(evecsv)
!x$OMP DO
         Do ik = 1, nkpt
            Allocate (evecsv(nstsv, nstsv))
            Call getevecsv (vkl(:, ik), evecsv)
! solve the Hartree-Fock secular equation
            Call seceqnhf (ik, evecsv)
! write the eigenvalues/vectors to file
            Call putevalsv (ik, evalsv(:, ik))
            Call putevecsv (ik, evecsv)
            Deallocate (evecsv)
         End Do
!x$OMP END DO
!x$OMP END PARALLEL
! find the occupation numbers and Fermi energy
         Call occupy
! write out the eigenvalues and occupation numbers
         Call writeeval
! write the Fermi energy to file
         Call writefermi
! set the charge density and magnetisation to zero
         rhomt (:, :, :) = 0.d0
         rhoir (:) = 0.d0
         If (associated(input%groundstate%spin)) Then
            magmt (:, :, :, :) = 0.d0
            magir (:, :) = 0.d0
         End If
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecfv,evecsv)
!$OMP DO
         Do ik = 1, nkpt
            Allocate (evecfv(nmatmax, nstfv, nspnfv))
            Allocate (evecsv(nstsv, nstsv))
! write the occupancies to file
            Call putoccsv (ik, occsv(:, ik))
! get the eigenvectors from file
            Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
            Call getevecsv (vkl(:, ik), evecsv)
! add to the density and magnetisation
            Call rhovalk (ik, evecfv, evecsv)
            Deallocate (evecfv, evecsv)
         End Do
!$OMP END DO
!$OMP END PARALLEL
! symmetrise the density
         Call symrf (input%groundstate%lradstep, rhomt, rhoir)
! symmetrise the magnetisation
         If (associated(input%groundstate%spin)) Call symrvf &
        & (input%groundstate%lradstep, magmt, magir)
! convert the density from a coarse to a fine radial mesh
         Call rfmtctof (rhomt)
! convert the magnetisation from a coarse to a fine radial mesh
         Do idm = 1, ndmag
            Call rfmtctof (magmt(:, :, :, idm))
         End Do
! add the core density to the total density
         Call addrhocr
! calculate the charges
         Call charge
! calculate the moments
         If (associated(input%groundstate%spin)) Call moment
! normalise the density
         Call rhonorm
! compute the Coulomb potential
         Call potcoul
! compute the energy components
         Call energy
! output energy components
         Call writeengy (60)
         Write (60,*)
         Write (60, '("Density of states at Fermi energy : ", G18.10)') &
        & fermidos
         Write (60, '(" (states/Hartree/unit cell)")')
! write total energy to TOTENERGY.OUT and flush
         Write (61, '(G22.12)') engytot
         Call flushifc (61)
! write DOS at Fermi energy to FERMIDOS.OUT and flush
         Write (62, '(G18.10)') fermidos
         Call flushifc (62)
! output charges and moments
         Call writechg (60,input%groundstate%outputlevelnumber)
! write total moment to MOMENT.OUT and flush
         If (associated(input%groundstate%spin)) Then
            Write (63, '(3G18.10)') momtot (1:ndmag)
            Call flushifc (63)
         End If
         If (tlast) Go To 20
! compute the change in total energy and check for convergence
         If (iscl .Ge. 2) Then
            de = Abs (engytot-etp) / (Abs(engytot)+1.d0)
            Write (60,*)
            Write (60, '("Relative change in total energy (target) : ",&
           & G18.10, " (", G18.10, ")")') de, &
           & input%groundstate%epsengy
            If (de .Lt. input%groundstate%epsengy) Then
               Write (60,*)
               Write (60, '("Energy convergence target achieved")')
               tlast = .True.
            End If
!            Write (65, '(G18.10)') de
!            Call flushifc (65)
         End If
         etp = engytot
! check for STOP file
         Inquire (File='STOP', Exist=Exist)
         If (exist) Then
            Write (60,*)
            Write (60, '("STOP file exists - stopping self-consistent l&
           &oop")')
            tlast = .True.
            Open (50, File='STOP')
            Close (50, Status='DELETE')
         End If
      End Do
20    Continue
      Write (60,*)
      Write (60, '("+------------------------------+")')
      Write (60, '("| Self-consistent loop stopped |")')
      Write (60, '("+------------------------------+")')
      If (input%groundstate%maxscl .Gt. 1) Then
         Call writestate
         Write (60,*)
         Write (60, '("Wrote STATE.OUT")')
      End If
!-----------------------!
!     compute forces    !
!-----------------------!
      If (( .Not. tstop) .And. (input%groundstate%tforce)) Then
         Call force
! output forces to INFO.OUT
         Call writeforce(60,input%groundstate%outputlevelnumber)
! write maximum force magnitude to FORCEMAX.OUT
!         Write (64, '(G18.10)') forcemax
!         Call flushifc (64)
      End If
      write(string,'("EXCITING ", a, " stopped")') trim(versionname)
      call printbox(60,"=",string)
! close the INFO.OUT file
      Close (60)
! close the TOTENERGY.OUT file
      Close (61)
! close the FERMIDOS.OUT file
      Close (62)
! close the MOMENT.OUT file
      If (associated(input%groundstate%spin)) close (63)
! close the FORCEMAX.OUT file
      If (input%groundstate%tforce) close (64)
! close the DENERGY.OUT file
      Close (65)

!----------------------------------------
! Save HF energies into binary file
!----------------------------------------

      Inquire (IoLength=Recl) nkpt, nstsv, vkl(:,1), evalsv(:,1)
      Open (70, File='EVALHF.OUT', Action='WRITE', Form='UNFORMATTED', &
     &   Access='DIRECT', status='REPLACE', Recl=Recl)
      do ik = 1, nkpt
        write(70, Rec=ik) nkpt, nstsv, vkl(:,ik), evalsv(:,ik)-efermi
      end do ! ik
      Close(70)
      
      Return
End Subroutine
