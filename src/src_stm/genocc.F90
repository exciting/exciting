Subroutine genocc(emin,emax)
  Use modmain
  Use modinput
  ! arguments
  Real (8), Intent (In) :: emin, emax
  ! external functions
  Real (8) :: sdelta, stheta
  External sdelta, stheta
  ! local variables
  Real (8) :: x,y,t1
  Integer :: ik, ist

  If (emin.eq.emax) Then
     t1 = 1.d0 / input%groundstate%swidth
     Do ik = 1, nkpt
        ! get the eigenvalues from file
        Call getevalsv (vkl(:, ik), evalsv(:, ik))
        Do ist = 1, nstsv
           x = (emin-evalsv(ist, ik)) * t1
           occsv (ist, ik) = occmax * wkpt (ik) * sdelta &
                & (input%groundstate%stypenumber, x) * t1
        End Do
     End Do
  Else
     Write(*,*) "emin,emax", emin, emax
     t1 = 1.d0 / input%groundstate%swidth
     Do ik = 1, nkpt
        Call getevalsv (vkl(:, ik), evalsv(:, ik))
        Do ist = 1, nstsv
           x = (emax - evalsv(ist, ik)) * t1
           y = (evalsv(ist, ik) - emin) * t1
           occsv (ist, ik) = occmax * wkpt (ik) * stheta &
                & (input%groundstate%stypenumber, x)*&
                & stheta(input%groundstate%stypenumber, y)
           !write(*,*) x,y,ist,ik,occsv(ist,ik)
        End Do
     End Do

  End If
End Subroutine genocc
