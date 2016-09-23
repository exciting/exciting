
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
! module to switch the different scl solver modes
!
Module sclcontroll
      Use modmain, Only: iscl, currentconvergence
      Use modinput
      Implicit None
 !scl index
      Integer :: diiscounter !! counter for DIIS iterations
      Logical :: packedmatrixstorage
      Logical :: tarpack, tlapack, tdiis, tjdqz
      Integer :: diisfirstscl
      Integer, Parameter :: diismax = 35, maxdiisspace = 15
      Integer :: iseed (4) = 1
      Real (8) :: lowesteval
      Real (8) :: epsarpack
      Real (8) :: epsresid
!
      Real (8), Parameter :: diisthreshould = 1
      Real (8) :: lastresnorm
      Integer, Parameter :: jacofidavidsonfirstscl = 1
      Integer :: idamax
      External idamax
      Logical :: recalculate_preconditioner
Contains
!
!
      Function calculate_preconditioner ()
         Logical :: calculate_preconditioner
         calculate_preconditioner = .False.
         If (diiscounter .Eq. 1) Then
            calculate_preconditioner = .True.
!
            Write (*,*) "precond"
         Else If (Mod(diiscounter, 5) .Eq. 0) Then
            If (currentconvergence .Gt. 5.0e-4) Then
               calculate_preconditioner = .True.
            End If
            Write (*,*) "precon"
         End If
      End Function calculate_preconditioner
!
!
      Function doDIIScycle ()
         Logical :: doDIIScycle
         doDIIScycle = .False.
         If (tdiis) Then
       !this may get more advanced:
            If (iscl .Ge. diisfirstscl) doDIIScycle = .True.
            If (currentconvergence .Gt. 1.0) doDIIScycle = .False.
            If (doDIIScycle) write (*,*) "DIIS"
         End If
         lastresnorm = 1.e10
      End Function doDIIScycle
!
!
      Function doprerotate_preconditioner ()
         Logical :: doprerotate_preconditioner
         doprerotate_preconditioner = .False.
         If (diiscounter .Ge. 3) doprerotate_preconditioner = .True.
      End Function doprerotate_preconditioner
!
!
      Function doARPACKiteration ()
         Logical :: doARPACKiteration
         doARPACKiteration = .False.
         If (associated(input%groundstate%solver)) Then
            If (input%groundstate%solver%typenumber .Eq. 2) Then
               doARPACKiteration = .True.
     ! write(*,*)"ARPACK"
               diiscounter = 1
            End If
         End If
!
      End Function doARPACKiteration
!
!
      Function doLAPACKsolver ()
         Logical :: doLAPACKsolver
         doLAPACKsolver = .False.
         If (associated(input%groundstate%solver)) Then
            If ((input%groundstate%solver%typenumber .Eq. 1)) Then
               doLAPACKsolver = .True.
!!!       write(*,*)"LAPACK hevx"
               diiscounter = 1
            End If
         Else
            doLAPACKsolver = .True.
         End If
!
      End Function doLAPACKsolver
!
!
      Function allconverged (n, rnorms)
         Logical :: allconverged
         Integer, Intent (In) :: n
         Real (8), Intent (In) :: rnorms (n)
         Real (8) :: rnormmax
         rnormmax = rnorms (idamax(n, rnorms, 1))
         If (rnormmax .Lt. epsresid) Then
            allconverged = .True.
            Write (*,*) " converged", rnorms (idamax(n, rnorms, 1))
         Else
            allconverged = .False.
            Write (*,*) "not converged", rnorms (idamax(n, rnorms, 1)), &
           & idamax (n, rnorms, 1)
         End If
!
         If (rnormmax/lastresnorm .Gt. 1.1) Then
 		!allconverged=.true.
            Write (*,*) "warning: error is gettig larger again", rnorms &
           & (idamax(n, rnorms, 1))
            If (rnormmax .Gt. .5e-6) Then
               allconverged = .False.
            Else
               Write (*,*) "error:error is gettig larger again"
      !  stop
            End If
!
         End If
         lastresnorm = rnormmax
      End Function allconverged
!
!
      Function dojacobdavidson ()
         Logical :: dojacobdavidson
         dojacobdavidson = .False.
         If (tjdqz) Then
            dojacobdavidson = .True.
            Write (*,*) "JDQZ"
         End If
      End Function
End Module sclcontroll
