!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_bandgap
      Implicit None
Contains
!
!BOP
! !ROUTINE: bandgap
! !INTERFACE:
!
!
      Subroutine bandgap (n, e, ef, egf, ego, ikgf, ikgo, istho)
! !USES:
         Use modmain
! !DESCRIPTION:
!   Determines the fundamental and optical band gap if present.
!
! !REVISION HISTORY:
!   Created July 2007 (Sagmeister)
!EOP
!BOC
         Implicit None
    ! arguments
         Integer, Intent (In) :: n
         Real (8), Intent (In) :: e (:, :), ef
         Real (8), Intent (Out) :: egf, ego
         Integer, Intent (Out) :: ikgf (2), ikgo, istho
    ! local variables
         Integer :: ik
         Integer :: klu1 (1), kho1 (1), kluho1 (1)
    ! allocatable arrays
         Real (8), Allocatable :: de (:), eho (:), elu (:)
!
         Allocate (de(nkpt), eho(nkpt), elu(nkpt))
         Do ik = 1, n
            istho = count (e(:, ik) <= ef)
            If (istho == nstfv) Go To 10
            eho (ik) = e (istho, ik)
            elu (ik) = e (istho+1, ik)
            de (ik) = elu (ik) - eho (ik)
         End Do
         kho1 = maxloc (eho)
         klu1 = minloc (elu)
         kluho1 = minloc (elu-eho)
         ikgf (1) = kho1 (1)
         ikgf (2) = klu1 (1)
         ikgo = kluho1 (1)
         egf = elu (ikgf(2)) - eho (ikgf(1))
         ego = elu (ikgo) - eho (ikgo)
         Return
10       Continue
    ! all states occupied
         egf = 0.0
         ego = 0.0
         ikgf = 1
         ikgo = 1
         istho = nstfv
         Deallocate (de, eho, elu)
!
      End Subroutine bandgap
!EOC
End Module m_bandgap
!
!
Subroutine writebandgap
      Use modmain
      Use modinput
      Use modxs
      Use m_bandgap
      Implicit None
  ! local variables
      Real (8) :: egf, ego, eho, elu, de, v1 (3), v2 (3)
      Integer :: ikgf (2), ikgo, istho, iv, ik
      Logical, Allocatable :: done (:)
      Real (8), External :: r3dist
!
  ! calculate bandgap
      Call bandgap (nkpt, evalsv, efermi, egf, ego, ikgf, ikgo, istho)
!
  ! write band gap to file
      If (egf == 0.d0) Go To 10
      If (task == 23) Then
         Open (50, File='BANDGAP_GRID.OUT', Form='formatted', Action='w&
        &rite', Status='replace')
         Write (50, '(a)') 'Band gaps determined from energies on grid'
         Write (50, '(a, 3i6)') ' k-point grid:', &
        & input%groundstate%ngridk
         Write (50, '(a, 3f12.6)') ' k-point offset:', &
        & input%groundstate%vkloff
         Write (50, '(a, 3f12.6)') ' k - point shift :', &
        & input%groundstate%vkloff / dble (input%groundstate%ngridk)
      Else If (task == 20) Then
         Open (50, File='BANDGAP.OUT', Form='formatted', Action='write',&
        &  Status='replace')
         Write (50, '(a)') 'Band gaps determined from energies on k-poi&
        &nt path'
         Write (50, '(a)') 'energies and differences on vertex location&
        &s'
         Write (50, '(a)') 'iv, ik, vkl, homo, lumo, de, de[eV] below'
         Allocate (done(nvp1d))
         done (:) = .False.
         Do iv = 1, nvp1d
            v1 (:) = vvlp1d (:, iv)
            Do ik = 1, nkpt
               v2 (:) = vkl (:, ik)
               If ( .Not. done(iv) .And. (r3dist(v1, v2) < &
              & input%structure%epslat)) Then
                  eho = evalsv (istho, ik)
                  elu = evalsv (istho+1, ik)
                  de = elu - eho
                  Write (50, '(2i6, 7f12.3)') iv, ik, vkl (:, ik), eho, &
                 & elu, de, h2ev * de
                  done (iv) = .True.
               End If
            End Do
         End Do
      End If
      Write (50, '(a, g16.8, a, g16.8, a)') 'fundamental gap: ', egf, '&
     & (', h2ev * egf, ' eV )'
      Write (50, '(a, i9, 4f12.6)') ' k - point (homo), energy: ', ikgf &
     & (1), vkl (:, ikgf(1)), evalsv (istho, ikgf(1))
      Write (50, '(a, i9, 4f12.6)') ' k - point (lumo), energy: ', ikgf &
     & (2), vkl (:, ikgf(2)+1), evalsv (istho+1, ikgf(2))
      Write (50, '(a, g16.8, a, g16.8, a)') 'optical gap    : ', ego, '&
     & (', h2ev * ego, ' eV )'
      Write (50, '(a, i9, 4f12.6)') ' k - point (homo), energy: ', &
     & ikgo, vkl (:, ikgo), evalsv (istho, ikgo)
      Write (50, '(a, i9, 4f12.6)') ' k - point (lumo), energy: ', &
     & ikgo, vkl (:, ikgo), evalsv (istho+1, ikgo)
      Close (50)
      Return
!
10    Continue
      Write (50, '(a)') 'No bandgaps found.'
!
End Subroutine writebandgap
!
!
Subroutine writebandgapgrid
      Use modmain
      Use modxs
      Use m_genfilname
      Implicit None
  ! local variables
      Integer :: ik
  ! initialise universal variables
      Call init0
  ! file extension for q-point
      Call genfilname (iqmt=0, setfilext=.True.)
      Call init1
  ! read Fermi energy from file
      Call readfermi
      Do ik = 1, nkpt
     ! get the eigenvectors and values from EIGVEC.OUT
         Call getevalsv (vkl(1, ik), evalsv(1, ik))
         Call getoccsv (vkl(1, ik), occsv(1, ik))
      End Do
      Call writebandgap
      Call genfilname (revertfilext=.True.)
End Subroutine writebandgapgrid
