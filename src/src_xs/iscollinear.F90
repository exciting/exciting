
! Crated Jul 2013 SR

Logical Function iscollinear ()
      Use modmain
  ! local variables
      iscollinear = .False.
      If (ndmag .Eq. 1 .And. ( associated(input%groundstate%spin) .Or. (ldapu .Ne. 0)) ) Then
        iscollinear = .True.
      End If
End Function iscollinear
