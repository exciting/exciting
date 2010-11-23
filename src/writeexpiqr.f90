!
!
!
!
Subroutine writeexpiqr
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: nk, ik, jk, i, j
      Real (8) :: vecqc (3), a, b
! allocatable arrays
      Complex (8), Allocatable :: emat (:, :)
! initialise universal variables
      Call init0
      Call init1
! allocate the matrix elements array for < i,k+G+q | exp(iq.r) | j,k >
      Allocate (emat(nstsv, nstsv))
! read in the density and potentials from file
      Call readstate
! find the new linearisation energies
      Call linengy
! generate the APW radial functions
      Call genapwfr
! generate the local-orbital radial functions
      Call genlofr
! number of k-points to write out
      nk = nkpt
      if (associated(input%properties%expiqr%kstlist)) &
         nk = size(input%properties%expiqr%kstlist%pointstatepair,2)
      Open (50, File='EXPIQR.OUT', Action='WRITE', Form='FORMATTED')
      Write (50,*)
      Write (50, '("q-vector (lattice coordinates) :")')
      Write (50, '(3G18.10)') input%properties%elnes%vecql
      Write (50, '("q-vector (Cartesian coordinates) :")')
      Call r3mv (bvec, input%properties%elnes%vecql, vecqc)
      Write (50, '(3G18.10)') vecqc
      Write (50,*)
      Write (50, '(I8," : number of k-points")') nk
      Write (50, '(I6," : number of states per k-point")') nstsv
      Do jk = 1, nk
         If (associated(input%properties%expiqr%kstlist)) Then
            ik = input%properties%expiqr%kstlist%pointstatepair(1, jk)
         Else
            ik = jk
         End If
         If ((ik .Le. 0) .Or. (ik .Gt. nkpt)) Then
            Write (*,*)
            Write (*, '("Error(writeexpiqr): k-point out of range : ",I&
           &8)') ik
            Write (*,*)
            Stop
         End If
         Write (50,*)
         Write (50, '(" k-point (lattice coordinates) :")')
         Write (50, '(3G18.10)') vkl (:, ik)
         Write (50,*)
         Write (50, '(" k-point (Cartesian coordinates) :")')
         Write (50, '(3G18.10)') vkc (:, ik)
         Call genexpiqr (ik, emat)
         Do i = 1, nstsv
            Write (50,*)
            Write (50, '(I6," : state i; state j, <...>, |<...>|^2 belo&
           &w")') i
            Do j = 1, nstsv
               a = dble (emat(i, j))
               b = aimag (emat(i, j))
               Write (50, '(I6,3G18.10)') j, a, b, a ** 2 + b ** 2
            End Do
         End Do
! end loop over k-points
      End Do
      Close (50)
      Write (*,*)
      Write (*, '("Info(writeexpiqr)")')
      Write (*, '(" < i,k+q | exp(iq.r) | j,k > matrix elements written&
     & to EXPIQR.OUT")')
      Write (*, '(" for the q-vector in vecql and all k-points in kstli&
     &st")')
      Write (*,*)
      Deallocate (emat)
      Return
End Subroutine
