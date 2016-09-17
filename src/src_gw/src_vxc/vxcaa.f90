
subroutine vxcaa(is, ia, ngp, apwalm, v, h)

    use modinput
    use modmain, only : idxas, idxlm, apword, gntyry, &
    &                   ngkmax, apwordmax, lmmaxapw, natmtot, zone
    use mod_vxc, only : vxcraa
    implicit none
    ! arguments
    Integer, Intent(In) :: is
    Integer, Intent(In) :: ia
    Integer, Intent(In) :: ngp
    complex(8), Intent(In) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
    complex(8), Intent(In) :: v(*)
    complex(8), Intent(Inout) :: h(*)
    ! local variables
    Integer :: ias, io1, io2
    Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
    Complex(8) :: zt1, zsum
    ! automatic arrays
    Complex(8) :: zv (ngp)

    ias = idxas (ia, is)
    Do l1 = 0, input%groundstate%lmaxmat
      Do m1 = - l1, l1
        lm1 = idxlm (l1, m1)
        Do io1 = 1, apword (l1, is)
          zv (:) = 0.d0
          Do l3 = 0, input%groundstate%lmaxmat
            Do m3 = - l3, l3
              lm3 = idxlm (l3, m3)
                If (lm1 .Ge. lm3) Then
                  Do io2 = 1, apword (l3, is)
                    zsum = 0.d0
                    Do l2 = 0, input%groundstate%lmaxvr
                      If (Mod(l1+l2+l3, 2) .Eq. 0) Then
                        Do m2 = - l2, l2
                          lm2 = idxlm (l2, m2)
                          If ((l2 .Eq. 0) .Or. (l1 .Ge. l3)) Then
                            zt1 = gntyry(lm1, lm2, lm3) * &
                            &     vxcraa(io1, l1, io2, l3, lm2, ias)
                          Else
                            zt1 = gntyry(lm1, lm2, lm3) * &
                            &     vxcraa(io1, l3, io2, l1, lm2, ias)
                          End If
                          zsum = zsum + zt1
                        End Do
                      End If
                    End Do
                    If (lm1 .Eq. lm3) zsum = zsum * 0.5d0
                    If (Abs(dble(zsum))+Abs(aimag(zsum)) > 1.d-14) Then
                      Call zaxpy (ngp, zsum, apwalm(:, io2, lm3, ias), 1, zv, 1)
                    End If
                  End Do
                End If
              End Do
            End Do
            Call zmatinp (.true., ngp, zone, apwalm(:, io1, lm1, ias), zv, v, h)
          End Do
        End Do
      End Do

      Return
End Subroutine
