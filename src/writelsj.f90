!
!
!
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma, C. Ambrosch-Draxl and
! F. Cricchio. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.
!
!
Subroutine writelsj
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: kst, ik, ist, lm
      Integer :: ispn, is, ia, ias
      Real (8) :: xl (3), xs (3), t1
! allocatable arrays
      Complex (8), Allocatable :: apwalm (:, :, :, :, :)
      Complex (8), Allocatable :: evecfv (:, :, :)
      Complex (8), Allocatable :: evecsv (:, :)
      Complex (8), Allocatable :: dmat1 (:, :, :, :, :)
      Complex (8), Allocatable :: dmat2 (:, :, :, :, :)
      Complex (8), Allocatable :: zlflm (:, :)
      Character(256) :: string
! initialise universal variables
      Call init0
      Call init1
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
      Allocate (evecfv(nmatmax, nstfv, nspnfv))
      Allocate (evecsv(nstsv, nstsv))
      Allocate (dmat1(lmmaxapw, lmmaxapw, nspinor, nspinor, natmtot))
      Allocate (dmat2(lmmaxapw, lmmaxapw, nspinor, nspinor, nstsv))
      Allocate (zlflm(lmmaxapw, 3))
! read density and potentials from file
        If (associated(input%groundstate%Hybrid)) Then
           If (input%groundstate%Hybrid%exchangetypenumber == 1) Then
! in case of HF hybrids use PBE potential
            string=filext
            filext='_PBE.OUT'
            Call readstate
            filext=string
           Else
               Call readstate
           End If
        Else         
           Call readstate
        End If 
! find the new linearisation energies
      Call linengy
! generate the APW radial functions
      Call genapwfr
! generate the local-orbital radial functions
      Call genlofr
! update potential in case if HF Hybrids
        If (associated(input%groundstate%Hybrid)) Then
           If (input%groundstate%Hybrid%exchangetypenumber == 1) Then
               Call readstate
           End If
        End If 
      If (task .Eq. 15) Then
! compute total L, S and J
         dmat1 (:, :, :, :, :) = 0.d0
         Do ik = 1, nkpt
! get the eigenvectors and occupancies from file
            Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
            Call getevecsv (vkl(:, ik), evecsv)
            Call getoccsv (vkl(:, ik), occsv(:, ik))
! find the matching coefficients
            Do ispn = 1, nspnfv
               Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, &
              & ispn, ik), sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, &
              & ispn))
            End Do
! loop over species and atoms
            Do is = 1, nspecies
               Do ia = 1, natoms (is)
                  ias = idxas (ia, is)
! generate the density matrix
                  Call gendmat (.False., .False., 0, &
                 & input%groundstate%lmaxapw, is, ia, ngk(:, ik), &
                 & apwalm, evecfv, evecsv, lmmaxapw, dmat2)
                  Do ist = 1, nstsv
                     t1 = wkpt (ik) * occsv (ist, ik)
                     dmat1 (:, :, :, :, ias) = dmat1 (:, :, :, :, ias) &
                    & + t1 * dmat2 (:, :, :, :, ist)
                  End Do
               End Do
            End Do
! end loop over k-points
         End Do
! symmetrise the density matrix
         Call symdmat (input%groundstate%lmaxapw, lmmaxapw, dmat1)
         Open (50, File='LSJ.OUT', Action='WRITE', Form='FORMATTED')
         Write (50,*)
         Write (50, '("Expectation values are computed only over the mu&
        &ffin-tin")')
! loop over species and atoms
         Do is = 1, nspecies
            Write (50,*)
            Write (50, '("Species : ", I4, " (", A, ")")') is, trim &
           & (input%structure%speciesarray(is)%species%chemicalSymbol)
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
! compute tr(LD)
               xl (:) = 0.d0
               Do ispn = 1, nspinor
                  Do lm = 1, lmmaxapw
                     Call lopzflm (input%groundstate%lmaxapw, dmat1(:, &
                    & lm, ispn, ispn, ias), lmmaxapw, zlflm)
                     xl (:) = xl (:) + dble (zlflm(lm, :))
                  End Do
               End Do
! compute tr(sigma D)
               xs (:) = 0.d0
               If (associated(input%groundstate%spin)) Then
                  Do lm = 1, lmmaxapw
                     xs (1) = xs (1) + dble (dmat1(lm, lm, 2, 1, &
                    & ias)+dmat1(lm, lm, 1, 2, ias))
                     xs (2) = xs (2) + dble (-zi*dmat1(lm, lm, 2, 1, &
                    & ias)+zi*dmat1(lm, lm, 1, 2, ias))
                     xs (3) = xs (3) + dble (dmat1(lm, lm, 1, 1, &
                    & ias)-dmat1(lm, lm, 2, 2, ias))
                  End Do
               End If
! S = 1/2 sigma
               xs (:) = 0.5d0 * xs (:)
               Write (50, '(" atom : ", I4)') ia
               Write (50, '("  L : ", 3G18.10)') xl (:)
               Write (50, '("  S : ", 3G18.10)') xs (:)
               Write (50, '("  J : ", 3G18.10)') xl (:) + xs (:)
! end loop over atoms and species
            End Do
         End Do
         Close (50)
         Write (*,*)
         Write (*, '("Info(writelsj):")')
         Write (*, '(" total L, S and J expectation values written to LSJ.OUT")')
         Write (*,*)
      Else
! compute L, S and J for all states in kstlist
         Open (50, File='LSJ_KST.OUT', Action='WRITE', Form='FORMATTED')
         Write (50,*)
         Write (50, '("Expectation values are computed only over the mu&
        &ffin-tin")')
         Do kst = 1, size(input%properties%LSJ%kstlist%pointstatepair,2)
            ik = input%properties%LSJ%kstlist%pointstatepair(1,kst)
            ist = input%properties%LSJ%kstlist%pointstatepair(2,kst)
            If ((ik .Le. 0) .Or. (ik .Gt. nkpt)) Then
               Write (*,*)
               Write (*, '("Error(writelsj): k-point out of range : ", &
              &I8)') ik
               Write (*,*)
               Stop
            End If
            If ((ist .Le. 0) .Or. (ist .Gt. nstsv)) Then
               Write (*,*)
               Write (*, '("Error(writelsj): state out of range : ", I8&
              &)') ist
               Write (*,*)
               Stop
            End If
! get the eigenvectors and occupancies from file
            Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
            Call getevecsv (vkl(:, ik), evecsv)
            Call getoccsv (vkl(:, ik), occsv(:, ik))
! find the matching coefficients
            Do ispn = 1, nspnfv
               Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, &
              & ispn, ik), sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, &
              & ispn))
            End Do
! loop over species and atoms
            Do is = 1, nspecies
               Do ia = 1, natoms (is)
                  ias = idxas (ia, is)
! generate the density matrix
                  Call gendmat (.False., .False., 0, &
                 & input%groundstate%lmaxapw, is, ia, ngk(:, ik), &
                 & apwalm, evecfv, evecsv, lmmaxapw, dmat2)
! compute tr(LD)
                  xl (:) = 0.d0
                  Do ispn = 1, nspinor
                     Do lm = 1, lmmaxapw
                        Call lopzflm (input%groundstate%lmaxapw, &
                       & dmat2(:, lm, ispn, ispn, ist), lmmaxapw, &
                       & zlflm)
                        xl (:) = xl (:) + dble (zlflm(lm, :))
                     End Do
                  End Do
! compute tr(sigma D)
                  xs (:) = 0.d0
                  If (associated(input%groundstate%spin)) Then
                     Do lm = 1, lmmaxapw
                        xs (1) = xs (1) + dble (dmat2(lm, lm, 2, 1, &
                       & ist)+dmat2(lm, lm, 1, 2, ist))
                        xs (2) = xs (2) + dble (-zi*dmat2(lm, lm, 2, 1, &
                       & ist)+zi*dmat2(lm, lm, 1, 2, ist))
                        xs (3) = xs (3) + dble (dmat2(lm, lm, 1, 1, &
                       & ist)-dmat2(lm, lm, 2, 2, ist))
                     End Do
                  Else
                     xs (3) = 1.d0
                  End If
! S = 1/2 sigma
                  xs (:) = 0.5d0 * xs (:)
                  Write (50,*)
                  Write (50, '("k-point : ", I6, 3G18.10)') ik, vkl (:, &
                 & ik)
                  Write (50, '("state : ", I6)') ist
                  Write (50, '("species : ", I4, " (", A, "), atom : ", I4)') is, trim &
                 & (input%structure%speciesarray(is)%species%chemicalSymbol), ia
                  Write (50, '(" L : ", 3G18.10)') xl (:)
                  Write (50, '(" S : ", 3G18.10)') xs (:)
                  Write (50, '(" J : ", 3G18.10)') xl (:) + xs (:)
               End Do
            End Do
         End Do
         Close (50)
         Write (*,*)
         Write (*, '("Info(writelsj):")')
         Write (*, '(" L, S and J expectation values for each k-point a&
        &nd state")')
         Write (*, '("	in kstlist written to LSJ_KST.OUT")')
         Write (*,*)
      End If
      Deallocate (apwalm, evecfv, evecsv, dmat1, dmat2, zlflm)
      Return
End Subroutine
