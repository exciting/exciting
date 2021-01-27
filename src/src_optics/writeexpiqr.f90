!
!
!
!
Subroutine writeexpiqr
      Use modmain
      Use modinput
      Use FoX_wxml
      Implicit None
! local variables
      Integer :: nk, ik, jk, i, j
      Real (8) :: vecqc (3), a, b
      Type (xmlf_t), Save :: xf
! allocatable arrays
      Complex (8), Allocatable :: emat (:, :)
      Character(256) :: string
! initialise universal variables

      Call init0
      Call init1
! allocate the matrix elements array for < i,k+G+q | exp(iq.r) | j,k >
      Allocate (emat(nstsv, nstsv))
! read in the density and potentials from file
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
      Call xml_OpenFile("expiqr.xml", xf, replace=.True., &
            &pretty_print=.True.)
      Call xml_NewElement(xf, "expiqr")
      Call xml_NewElement(xf, "q-vector")
      Call xml_AddAttribute(xf, "vecql", input%properties%elnes%vecql)
      Call xml_AddAttribute(xf, "vecqc", vecqc)
      Call xml_EndElement(xf, "q-vector")
      Call xml_NewElement(xf, "k-grid")
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
         Call xml_NewElement(xf, "k-vector")
         Call xml_AddAttribute(xf, "vkl", vkl (:, ik))
         Call xml_AddAttribute(xf, "vkc", vkc (:, ik))
         Call genexpiqr (ik, emat)
         Do i = 1, nstsv
            Write (50,*)
            Write (50, '(I6," : state i; state j, <...>, |<...>|^2 belo&
           &w")') i
           Call xml_NewElement(xf, "state")
           Call xml_AddAttribute(xf, "i", i)
            Do j = 1, nstsv
               a = dble (emat(i, j))
               b = aimag (emat(i, j))
               Write (50, '(I6,3G18.10)') j, a, b, a ** 2 + b ** 2
               Call xml_NewElement(xf, "state")
               Call xml_AddAttribute(xf, "j", j)
               Call xml_AddAttribute(xf, "Re", a)
               Call xml_AddAttribute(xf, "Im", b)
               Call xml_AddAttribute(xf, "norm", a ** 2 + b ** 2)
               Call xml_EndElement(xf, "state")
            End Do
            Call xml_EndElement(xf, "state")
         End Do 
         Call xml_EndElement(xf, "k-vector")
! end loop over k-points
      End Do
      Call xml_EndElement(xf, "k-grid")
      Call xml_EndElement(xf, "expiqr")
      Call xml_Close(xf)
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
