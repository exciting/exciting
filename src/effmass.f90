!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine effmass
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer, Parameter :: lwork = 10
      Integer :: ik, ik0, ist, info
      Integer :: i, j, k, l, m, n
      Integer :: i1, i2, i3, j1, j2, j3
      Real (8) :: d (3, 3), em (3, 3)
      Real (8) :: v1 (3), v2 (3)
      Real (8) :: w (3), work (lwork)
! allocatable arrays
      Integer, Allocatable :: ipiv (:)
      Real (8), Allocatable :: a (:, :)
      Real (8), Allocatable :: b (:, :, :, :)
      Real (8), Allocatable :: c (:, :, :)
      Real (8), Allocatable :: evalfv (:, :)
      Complex (8), Allocatable :: evecfv (:, :, :)
      Complex (8), Allocatable :: evecsv (:, :)
! initialise universal variables
      Call init0
      Call init1
      Allocate (ipiv(nkpt))
      Allocate (a(nkpt, nkpt))
      n = 2 * input%properties%masstensor%ndspem + 1
      Allocate (b(0:n-1, 0:n-1, 0:n-1, nstsv))
      Allocate (c(0:n-1, 0:n-1, 0:n-1))
! read density and potentials from file
      Call readstate
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
      ik0 = 0
! begin parallel loop over k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv) &
!$OMP PRIVATE(i1,i2,i3,j1,j2,j3,ist)
!$OMP DO
      Do ik = 1, nkpt
         Allocate (evalfv(nstfv, nspnfv))
         Allocate (evecfv(nmatmax, nstfv, nspnfv))
         Allocate (evecsv(nstsv, nstsv))
         i1 = ivk (1, ik)
         i2 = ivk (2, ik)
         i3 = ivk (3, ik)
         If ((i1 .Eq. 0) .And. (i2 .Eq. 0) .And. (i3 .Eq. 0)) ik0 = ik
! solve the first- and second-variational secular equations
         Call seceqn (ik, evalfv, evecfv, evecsv)
! copy eigenvalues to new array
         j1 = i1 + input%properties%masstensor%ndspem
         j2 = i2 + input%properties%masstensor%ndspem
         j3 = i3 + input%properties%masstensor%ndspem
         Do ist = 1, nstsv
            b (j1, j2, j3, ist) = evalsv (ist, ik)
         End Do
         Deallocate (evalfv, evecfv, evecsv)
      End Do
!$OMP END DO
!$OMP END PARALLEL
! set up polynomial matrix
      i = 0
      Do i3 = - input%properties%masstensor%ndspem, &
     & input%properties%masstensor%ndspem
         Do i2 = - input%properties%masstensor%ndspem, &
        & input%properties%masstensor%ndspem
            Do i1 = - input%properties%masstensor%ndspem, &
           & input%properties%masstensor%ndspem
               i = i + 1
               v1 (1) = dble (i1)
               v1 (2) = dble (i2)
               v1 (3) = dble (i3)
               v1 (:) = v1 (:) * input%properties%masstensor%deltaem
               j = 0
               v2 (3) = 1.d0
               Do j3 = 0, n - 1
                  v2 (2) = 1.d0
                  Do j2 = 0, n - 1
                     v2 (1) = 1.d0
                     Do j1 = 0, n - 1
                        j = j + 1
                        a (i, j) = v2 (1) * v2 (2) * v2 (3)
                        v2 (1) = v2 (1) * v1 (1)
                     End Do
                     v2 (2) = v2 (2) * v1 (2)
                  End Do
                  v2 (3) = v2 (3) * v1 (3)
               End Do
            End Do
         End Do
      End Do
! solve for the polynomial coefficients
      Call dgesv (nkpt, nstsv, a, nkpt, ipiv, b, nkpt, info)
      If (info .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(effmass): could not determine polynomial co&
        &efficients")')
         Write (*, '(" DGESV returned INFO = ", I8)') info
         Write (*,*)
         Stop
      End If
      Open (50, File='EFFMASS.OUT', Action='WRITE', Form='FORMATTED')
      Write (50,*)
      Write (50, '("(effective mass matrices are in Cartesian coordinat&
     &es)")')
      Write (50,*)
      Write (50, '("k-point (lattice coordinates) :")')
      Write (50, '(3G18.10)') input%properties%masstensor%vklem
      Write (50,*)
      Write (50, '("k-point (Cartesian coordinates) :")')
      Call r3mv (bvec, input%properties%masstensor%vklem, v1)
      Write (50, '(3G18.10)') v1
! begin loop over states
      Do ist = 1, nstsv
! compute matrix of derivatives with respect to k-vector
         Do k = 1, 3
            Do l = 1, 3
               c (:, :, :) = b (:, :, :, ist)
               Do i = 1, 2
                  If (i .Eq. 1) Then
                     m = k
                  Else
                     m = l
                  End If
                  If (m .Eq. 1) Then
                     Do j = 0, n - 2
                        c (j, :, :) = dble (j+1) * c (j+1, :, :)
                     End Do
                     c (n-1, :, :) = 0.d0
                  Else If (m .Eq. 2) Then
                     Do j = 0, n - 2
                        c (:, j, :) = dble (j+1) * c (:, j+1, :)
                     End Do
                     c (:, n-1, :) = 0.d0
                  Else If (m .Eq. 3) Then
                     Do j = 0, n - 2
                        c (:, :, j) = dble (j+1) * c (:, :, j+1)
                     End Do
                     c (:, :, n-1) = 0.d0
                  End If
               End Do
! derivative evaluated at zero
               d (k, l) = c (0, 0, 0)
            End Do
         End Do
         Write (50,*)
         Write (50,*)
         Write (50, '("State, eigenvalue : ", I6, G18.10)') ist, evalsv &
        & (ist, ik0)
         Write (50,*)
         Write (50, '(" matrix of eigenvalue derivatives with respect t&
        &o k :")')
         Do i = 1, 3
            Write (50, '(3G18.10)') (d(i, j), j=1, 3)
         End Do
         Write (50, '(" trace : ", G18.10)') d (1, 1) + d (2, 2) + d &
        & (3, 3)
! invert derivative matrix
         Call r3minv (d, em)
         Write (50,*)
         Write (50, '(" effective mass tensor (inverse of derivative matrix) :")')
         Do i = 1, 3
            Write (50, '(3G18.10)') (em(i, j), j=1, 3)
         End Do
         Write (50, '(" trace : ", G18.10)') em (1, 1) + em (2, 2) + em &
        & (3, 3)
! find the eigenvalues
         Call dsyev ('N', 'U', 3, em, 3, w, work, lwork, info)
         Write (50, '(" eigenvalues :")')
         Write (50, '(3G18.10)') w
! end loop over states
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(effmass):")')
      Write (*, '(" effective mass tensor for each state written to EFF&
     &MASS.OUT")')
      Write (*, '(" for k-point (lattice) ", 3G18.10)') &
     & input%properties%masstensor%vklem
      Write (*,*)
      Deallocate (ipiv, a, b, c)
      Return
End Subroutine
