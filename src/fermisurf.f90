!
!
!
!
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and
!C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine fermisurf
      Use modmain
      Use modinput
      Use FoX_wxml
      Implicit None
  ! local variables
      Integer :: ik, jk, ist, i
      Integer :: ist0, ist1, nst
      Real (8) :: prd1, prd2
      Character (128) :: buffer
      Type (xmlf_t), Save :: xf
!
  ! allocatable arrays
      Real (8), Allocatable :: evalfv (:, :)
      Complex (8), Allocatable :: evecfv (:, :, :)
      Complex (8), Allocatable :: evecsv (:, :)
  ! initialise universal variables
!
      Call init0
      Call init1
  ! read density and potentials from file
      Call readstate
  ! read Fermi energy from file
      Call readfermi
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
  ! begin parallel loop over reduced k-points set
  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(evalfv,evecfv,evecsv)
  !$OMP DO
      Do ik = 1, nkpt
         Allocate (evalfv(nstfv, nspnfv))
         Allocate (evecfv(nmatmax, nstfv, nspnfv))
         Allocate (evecsv(nstsv, nstsv))
     !$OMP CRITICAL
         Write (*, '("Info(fermisurf): ", I6, " of ", I6, " k-points")') ik, nkpt
     !$OMP END CRITICAL
     ! solve the first- and second-variational secular equations
         Call seceqn (ik, evalfv, evecfv, evecsv)
         Deallocate (evalfv, evecfv, evecsv)
     ! end loop over reduced k-points set
      End Do
  !$OMP END DO
  !$OMP END PARALLEL
      Call xml_OpenFile ("fermisurface.xml", xf, replace=.True., &
     & pretty_print=.True.)
      Call xml_NewElement (xf, "fermisurface")
      Call xml_NewElement (xf, "runitcell")
      Write (buffer, '(3I6)') np3d (:)
      Call xml_addAttribute (xf, "grid", trim(adjustl(buffer)))
      Do i = 1, 3
         Call xml_NewElement (xf, "bvec")
         Write (buffer, '(4G18.10)') bvec (:, i)
         Call xml_AddCharacters (xf, trim(adjustl(buffer)))
         Call xml_endElement (xf, "bvec")
      End Do
      Call xml_endElement (xf, "runitcell")
      If (ndmag .Eq. 1) Then

     ! special case of collinear magnetism
         Open (50, File='FERMISURF_UP.OUT', Action='WRITE', Form='FORMA&
        &TTED')
         Open (51, File='FERMISURF_DN.OUT', Action='WRITE', Form='FORMA&
        &TTED')
         If (task .Eq. 100) Then

        ! write product of eigenstates minus the Fermi energy
            Write (50, '(3I6, " : grid size")') np3d (:)
            Write (51, '(3I6, " : grid size")') np3d (:)
            Do ik = 1, nkptnr
!
               jk = ikmap (ivknr(1, ik), ivknr(2, ik), ivknr(3, ik))
               prd1 = 1.d0
               prd2 = 1.d0
               Do ist = 1, nstfv
                  prd1 = prd1 * (evalsv(ist, jk)-efermi)
                  prd2 = prd2 * (evalsv(nstfv+ist, jk)-efermi)
               End Do
               Call xml_NewElement (xf, "point")
               Write (buffer, '(4G18.10)') vkcnr (1, ik)
               Call xml_addAttribute (xf, "x", trim(adjustl(buffer)))
               Write (buffer, '(4G18.10)') vkcnr (2, ik)
               Call xml_addAttribute (xf, "y", trim(adjustl(buffer)))
               Write (buffer, '(4G18.10)') vkcnr (3, ik)
               Call xml_addAttribute (xf, "z", trim(adjustl(buffer)))
               Write (buffer, '(4G18.10)') prd1
               Call xml_addAttribute (xf, "up", trim(adjustl(buffer)))
               Write (buffer, '(4G18.10)') prd2
               Call xml_addAttribute (xf, "down", &
              & trim(adjustl(buffer)))
               Call xml_endElement (xf, "point")
               Write (50, '(4G18.10)') vkcnr (:, ik), prd1
               Write (51, '(4G18.10)') vkcnr (:, ik), prd2
            End Do
         Else
        ! write the eigenvalues minus the Fermi energy separately

!
            ist = nstfv - input%groundstate%nempty
            ist0 = Max (ist-input%properties%fermisurfaceplot%nstfsp/2, &
           & 1)
            ist1 = Min (ist+input%properties%fermisurfaceplot%nstfsp/2, &
           & nstfv)
            nst = ist1 - ist0 + 1
            Write (50, '(4I6, " : grid size, number of states")') np3d &
           & (:), nst
            Write (51, '(4I6, " : grid size, number of states")') np3d &
           & (:), nst
!

            Do ik = 1, nkptnr
               jk = ikmap (ivknr(1, ik), ivknr(2, ik), ivknr(3, ik))
               Write (50, '(G18.10)', Advance='NO') vkcnr (:, ik)
               Write (51, '(G18.10)', Advance='NO') vkcnr (:, ik)
               Call xml_NewElement (xf, "point")
               Write (buffer, '(4G18.10)') vkcnr (1, ik)
               Call xml_addAttribute (xf, "x", trim(adjustl(buffer)))
               Write (buffer, '(4G18.10)') vkcnr (2, ik)
               Call xml_addAttribute (xf, "y", trim(adjustl(buffer)))
               Write (buffer, '(4G18.10)') vkcnr (3, ik)
               Call xml_addAttribute (xf, "z", trim(adjustl(buffer)))
               Do ist = ist0, ist1
                  Call xml_NewElement (xf, "band")
                  Write (buffer, '(4G18.10)') evalsv (ist, jk) - efermi
                  Call xml_addAttribute (xf, "evalup", &
                 & trim(adjustl(buffer)))
                  Write (buffer, '(4G18.10)') evalsv (nstfv+ist, jk) - &
                 & efermi
                  Call xml_addAttribute (xf, "evaldown", &
                 & trim(adjustl(buffer)))
                  Call xml_endElement (xf, "band")
                  Write (50, '(F14.8)', Advance='NO') evalsv (ist, jk) &
                 & - efermi
                  Write (51, '(F14.8)', Advance='NO') evalsv &
                 & (nstfv+ist, jk) - efermi
               End Do
               Call xml_endElement (xf, "point")
               Write (50,*)
               Write (51,*)
            End Do
         End If
         Close (50)
         Close (51)
      Else
     ! spin-unpolarised and non-collinear cases
         Open (50, File='FERMISURF.OUT', Action='WRITE', Form='FORMATTE&
        &D')
         If (task .Eq. 100) Then
        ! write product of eigenstates minus the Fermi energy
            Write (50, '(3I6, " : grid size")') np3d (:)
!
            Do ik = 1, nkptnr
               jk = ikmap (ivknr(1, ik), ivknr(2, ik), ivknr(3, ik))
               prd1 = 1.d0
               Do ist = 1, nstsv
                  prd1 = prd1 * (evalsv(ist, jk)-efermi)
               End Do
               Write (50, '(4G18.10)') vkcnr (:, ik), prd1
               Call xml_NewElement (xf, "point")
               Write (buffer, '(4G18.10)') vkcnr (1, ik)
               Call xml_addAttribute (xf, "x", trim(adjustl(buffer)))
               Write (buffer, '(4G18.10)') vkcnr (2, ik)
               Call xml_addAttribute (xf, "y", trim(adjustl(buffer)))
               Write (buffer, '(4G18.10)') vkcnr (3, ik)
               Call xml_addAttribute (xf, "z", trim(adjustl(buffer)))
               Write (buffer, '(4G18.10)') prd1
               Call xml_addAttribute (xf, "product", &
              & trim(adjustl(buffer)))
               Call xml_endElement (xf, "point")
            End Do
         Else
        ! write the eigenvalues minus the Fermi energy separately
            ist = (nstfv-input%groundstate%nempty) * nspinor
            ist0 = Max (ist-input%properties%fermisurfaceplot%nstfsp/2, &
           & 1)
            ist1 = Min (ist+input%properties%fermisurfaceplot%nstfsp/2, &
           & nstsv)
            nst = ist1 - ist0 + 1
            Write (50, '(4I6, " : grid size, number of states")') np3d &
           & (:), nst
!
            Do ik = 1, nkptnr
               Call xml_NewElement (xf, "point")
               Write (buffer, '(4G18.10)') vkcnr (1, ik)
               Call xml_addAttribute (xf, "x", trim(adjustl(buffer)))
               Write (buffer, '(4G18.10)') vkcnr (2, ik)
               Call xml_addAttribute (xf, "y", trim(adjustl(buffer)))
               Write (buffer, '(4G18.10)') vkcnr (3, ik)
               Call xml_addAttribute (xf, "z", trim(adjustl(buffer)))
!
               jk = ikmap (ivknr(1, ik), ivknr(2, ik), ivknr(3, ik))
               Write (50, '(3G18.10)', Advance='NO') vkcnr (:, ik)
               Do ist = ist0, ist1
                  Call xml_NewElement (xf, "band")
                  Write (buffer, '(4G18.10)') evalsv (ist, jk) - efermi
                  Call xml_addAttribute (xf, "eval", &
                 & trim(adjustl(buffer)))
                  Write (buffer,*) ist
                  Call xml_addAttribute (xf, "nr", &
                 & trim(adjustl(buffer)))
                  Call xml_endElement (xf, "band")
                  Write (50, '(F14.8)', Advance='NO') evalsv (ist, jk) &
                 & - efermi
               End Do
               Call xml_endElement (xf, "point")
               Write (50,*)
            End Do
         End If
         Close (50)
      End If
      Call xml_close (xf)
      Write (*,*)
      Write (*, '("Info(fermisurf):")')
      If (ndmag .Eq. 1) Then
         Write (*, '(" 3D Fermi surface data written to FERMISURF_UP.OU&
        &T and FERMISURF_DN.OUT")')
      Else
         Write (*, '(" 3D Fermi surface data written to FERMISURF.OUT")&
        &')
      End If
      If (task .Eq. 100) Then
         Write (*, '(" in terms of the product of eigenvalues minus the&
        & Fermi energy")')
      Else
         Write (*, '(" in terms of separate eigenvalues minus the Fermi&
        & energy")')
      End If
      Write (*,*)
      Return
End Subroutine fermisurf
