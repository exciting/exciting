!
!
!
!
Subroutine writelambda (wq, gq)
      Use modmain
      Implicit None
! arguments
      Real (8), Intent (In) :: wq (3*natmtot, nqpt)
      Real (8), Intent (In) :: gq (3*natmtot, nqpt)
! local variables
      Integer :: ik, iq, i
      Real (8) :: t1, t2
! get the eigenvalues and occupancies from file
      Do ik = 1, nkpt
         Call getevalsv (vkl(:, ik), evalsv(:, ik))
         Call getoccsv (vkl(:, ik), occsv(:, ik))
      End Do
! compute the density of states at the Fermi energy
      Call occupy
      Open (50, File='LAMBDAQ.OUT', Action='WRITE', Form='FORMATTED')
      Write (50,*)
      Write (50, '(I4, " : total number of atoms")') natmtot
      Write (50, '(I6, " : number of q-points")') nqpt
      Write (50,*)
      Do iq = 1, nqpt
         Write (50, '(I6, " : q-point")') iq
         Write (50, '(3G18.10, " : q-vector (lattice coordinates)")') &
        & vql (:, iq)
         Write (50, '(3G18.10, " : q-vector (Cartesian coordinates)")') &
        & vqc (:, iq)
         Do i = 1, 3 * natmtot
            t1 = pi * fermidos * wq (i, iq) ** 2
            If (t1 .Gt. 1.d-8) Then
               t2 = gq (i, iq) / t1
            Else
               t2 = 0.d0
            End If
            Write (50, '(I4, G18.10)') i, t2
         End Do
         Write (50,*)
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(writelambda):")')
      Write (*, '(" wrote electron-phonon coupling constants for all q-&
     &points to LAMBDAQ.OUT")')
      Write (*,*)
      Return
End Subroutine
