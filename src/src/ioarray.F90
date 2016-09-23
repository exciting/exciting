
! Copyright (C) 2007-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

Module ioarray
      Implicit None
      Integer, Parameter :: kndi = 4
      Integer, Parameter :: kndr = 8
      Integer, Parameter :: kndc = 8
      Integer, Parameter :: fmtlen = 1024
      Integer, Parameter :: und = 1
      Logical, Parameter :: tparend = .True.
      Character (*), Parameter :: fmtdl = 'l2, 1x'
      Character (*), Parameter :: fmtdi = 'i12'
      Character (*), Parameter :: fmtdr = 'es23.15E3'
      Character (*), Parameter :: fmtdc = 'es23.15E3'
      Private :: kndi, kndr, kndc, fmtlen, und
      Private :: fmtdl, fmtdi, fmtdr, fmtdc
      Private :: ioaction

Contains

      Subroutine ioarr (arr1dl, arr1di, arr1dr, arr1dc, arr2dl, arr2di, &
     & arr2dr, arr2dc, arr3dl, arr3di, arr3dr, arr3dc, arr4dl, arr4di, &
     & arr4dr, arr4dc, arr5dl, arr5di, arr5dr, arr5dc, arr6dl, arr6di, &
     & arr6dr, arr6dc, arr7dl, arr7di, arr7dr, arr7dc, ioa, un, fmt, &
     & fmtidx, tparen, header)
         Implicit None
    ! arguments
         Character (*), Intent (In) :: ioa
    ! * 1D arrays *
         Logical, Optional :: arr1dl (:)
         Integer (kndi), Optional :: arr1di (:)
         Real (kndr), Optional :: arr1dr (:)
         Complex (kndc), Optional :: arr1dc (:)
    ! * 2D arrays *
         Logical, Optional :: arr2dl (:, :)
         Integer (kndi), Optional :: arr2di (:, :)
         Real (kndr), Optional :: arr2dr (:, :)
         Complex (kndc), Optional :: arr2dc (:, :)
    ! * 3D arrays *
         Logical, Optional :: arr3dl (:, :, :)
         Integer (kndi), Optional :: arr3di (:, :, :)
         Real (kndr), Optional :: arr3dr (:, :, :)
         Complex (kndc), Optional :: arr3dc (:, :, :)
    ! * 4D arrays *
         Logical, Optional :: arr4dl (:, :, :, :)
         Integer (kndi), Optional :: arr4di (:, :, :, :)
         Real (kndr), Optional :: arr4dr (:, :, :, :)
         Complex (kndc), Optional :: arr4dc (:, :, :, :)
    ! * 5D arrays *
         Logical, Optional :: arr5dl (:, :, :, :, :)
         Integer (kndi), Optional :: arr5di (:, :, :, :, :)
         Real (kndr), Optional :: arr5dr (:, :, :, :, :)
         Complex (kndc), Optional :: arr5dc (:, :, :, :, :)
    ! * 6D arrays *
         Logical, Optional :: arr6dl (:, :, :, :, :, :)
         Integer (kndi), Optional :: arr6di (:, :, :, :, :, :)
         Real (kndr), Optional :: arr6dr (:, :, :, :, :, :)
         Complex (kndc), Optional :: arr6dc (:, :, :, :, :, :)
    ! * 7D arrays *
         Logical, Optional :: arr7dl (:, :, :, :, :, :, :)
         Integer (kndi), Optional :: arr7di (:, :, :, :, :, :, :)
         Real (kndr), Optional :: arr7dr (:, :, :, :, :, :, :)
         Complex (kndc), Optional :: arr7dc (:, :, :, :, :, :, :)
    ! * other arguments *
         Integer, Optional, Intent (In) :: un
         Character (*), Optional, Intent (In) :: fmt
         Character (*), Optional, Intent (In) :: fmtidx
         Logical, Optional, Intent (In) :: tparen
         Character (*), Optional, Intent (In) :: header
    ! local variables
         Integer :: npr, unt, ndim
         Integer :: i1, i2, i3, i4, i5, i6, i7
         Integer :: j1, j2, j3, j4, j5, j6, j7
         Integer :: lb (7), ub (7)
         Character (1) :: iot
         Character (fmtlen) :: frmt, str
         Logical :: tpr (28), twrite, tparent, tfmtidx
         Logical :: lt
         Integer (kndi) :: it
         Real (kndr) :: rt
         Complex (kndc) :: zt
    ! check if only one array is present
         tpr = (/ present (arr1dl), present (arr1di), present (arr1dr), &
        & present (arr1dc), present (arr2dl), present (arr2di), present &
        & (arr2dr), present (arr2dc), present (arr3dl), present &
        & (arr3di), present (arr3dr), present (arr3dc), present &
        & (arr4dl), present (arr4di), present (arr4dr), present &
        & (arr4dc), present (arr5dl), present (arr5di), present &
        & (arr5dr), present (arr5dc), present (arr6dl), present &
        & (arr6di), present (arr6dr), present (arr6dc), present &
        & (arr7dl), present (arr7di), present (arr7dr), present &
        & (arr7dc) /)
         npr = count (tpr)
         If (npr .Ne. 1) Then
            Write (*,*)
            Write (*,*) 'Error(ioarray::ioarr): more than one array spe&
           &cified'
            Write (*,*)
            Stop
         End If
    ! determine dimension
         If (any(tpr(1:4))) ndim = 1
         If (any(tpr(5:8))) ndim = 2
         If (any(tpr(9:12))) ndim = 3
         If (any(tpr(13:16))) ndim = 4
         If (any(tpr(17:20))) ndim = 5
         If (any(tpr(21:24))) ndim = 6
         If (any(tpr(25:28))) ndim = 7
    ! determine type
         If (any(tpr(1::4))) iot = 'l'
         If (any(tpr(2::4))) iot = 'i'
         If (any(tpr(3::4))) iot = 'r'
         If (any(tpr(4::4))) iot = 'c'
    ! determine IO action (read or write)
         twrite = ioaction (ioa)
    ! optional parentheses for complex numbers
         tparent = tparend
         If (present(tparen)) tparent = tparen
         If (iot .Ne. 'c') tparent = .False.
    ! determine file unit
         unt = und
         If (present(un)) unt = un
    ! determine format
         Select Case (iot)
         Case ('l')
            frmt = trim (fmtdl)
         Case ('i')
            frmt = trim (fmtdi)
         Case ('r')
            frmt = trim (fmtdr)
         Case ('c')
            frmt = trim (fmtdc)
         End Select
         If (present(fmt)) frmt = trim (adjustl(fmt))
         tfmtidx = .False.
         If (present(fmtidx)) tfmtidx = .True.
    ! write header
         If (present(header)) write (unit=unt, fmt='(a)') trim (header)
    ! select by dimension
         Select Case (ndim)
         Case (1)
       !---------------------!
       !     1 dimension     !
       !---------------------!
            Select Case (iot)
            Case ('l')
               lb (1:1) = lbound (arr1dl)
               ub (1:1) = ubound (arr1dl)
               str = '(' // fmtdi // ', ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  If (twrite) Then
                     lt = arr1dl (i1)
                     Write (Unit=unt, Fmt=trim(str)) i1, lt
                  Else
                     Read (Unit=unt, Fmt=*) j1, lt
                     arr1dl (j1) = lt
                  End If
               End Do
            Case ('i')
               lb (1:1) = lbound (arr1di)
               ub (1:1) = ubound (arr1di)
               str = '(' // fmtdi // ', ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  If (twrite) Then
                     it = arr1di (i1)
                     Write (Unit=unt, Fmt=trim(str)) i1, it
                  Else
                     Read (Unit=unt, Fmt=*) j1, it
                     arr1di (j1) = it
                  End If
               End Do
            Case ('r')
               lb (1:1) = lbound (arr1dr)
               ub (1:1) = ubound (arr1dr)
               str = '(' // fmtdi // ', " ", ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', " ", ' // trim &
              & (frmt) // ')'
               Do i1 = lb (1), ub (1)
                  If (twrite) Then
                     rt = arr1dr (i1)
                     Write (Unit=unt, Fmt=trim(str)) i1, rt
                  Else
                     Read (Unit=unt, Fmt=*) j1, rt
                     arr1dr (j1) = rt
                  End If
               End Do
            Case ('c')
               lb (1:1) = lbound (arr1dc)
               ub (1:1) = ubound (arr1dc)
               str = '(' // fmtdi // ', " ", ' // trim (frmt) // ', " "&
              &, ' // trim (frmt) // ')'
               If (tparent) str = '(' // fmtdi // ', " (", ' // trim &
              & (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               If (tfmtidx) Then
                  str = '(' // fmtidx // ', " ", ' // trim (frmt) // ',&
                 & " ", ' // trim (frmt) // ')'
                  If (tparent) str = '(' // fmtidx // ', " (", ' // &
                 & trim (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               End If
               Do i1 = lb (1), ub (1)
                  If (twrite) Then
                     zt = arr1dc (i1)
                     Write (Unit=unt, Fmt=trim(str)) i1, zt
                  Else
                     Read (Unit=unt, Fmt=*) j1, zt
                     arr1dc (j1) = zt
                  End If
               End Do
            End Select
         Case (2)
       !----------------------!
       !     2 dimensions     !
       !----------------------!
            Select Case (iot)
            Case ('l')
               lb (1:2) = lbound (arr2dl)
               ub (1:2) = ubound (arr2dl)
               str = '(2' // fmtdi // ', ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     If (twrite) Then
                        lt = arr2dl (i1, i2)
                        Write (Unit=unt, Fmt=trim(str)) i1, i2, lt
                     Else
                        Read (Unit=unt, Fmt=*) j1, j2, lt
                        arr2dl (j1, j2) = lt
                     End If
                  End Do
               End Do
            Case ('i')
               lb (1:2) = lbound (arr2di)
               ub (1:2) = ubound (arr2di)
               str = '(2' // fmtdi // ', ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     If (twrite) Then
                        it = arr2di (i1, i2)
                        Write (Unit=unt, Fmt=trim(str)) i1, i2, it
                     Else
                        Read (Unit=unt, Fmt=*) j1, j2, it
                        arr2di (j1, j2) = it
                     End If
                  End Do
               End Do
            Case ('r')
               lb (1:2) = lbound (arr2dr)
               ub (1:2) = ubound (arr2dr)
               str = '(2' // fmtdi // ', " ", ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     If (twrite) Then
                        rt = arr2dr (i1, i2)
                        Write (Unit=unt, Fmt=trim(str)) i1, i2, rt
                     Else
                        Read (Unit=unt, Fmt=*) j1, j2, rt
                        arr2dr (j1, j2) = rt
                     End If
                  End Do
               End Do
            Case ('c')
               lb (1:2) = lbound (arr2dc)
               ub (1:2) = ubound (arr2dc)
               str = '(2' // fmtdi // ', " ", ' // trim (frmt) // ', " &
              &", ' // trim (frmt) // ')'
               If (tparent) str = '(2' // fmtdi // ', " (", ' // trim &
              & (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               If (tfmtidx) Then
                  str = '(' // fmtidx // ', " ", ' // trim (frmt) // ',&
                 & " ", ' // trim (frmt) // ')'
                  If (tparent) str = '(' // fmtidx // ', " (", ' // &
                 & trim (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               End If
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     If (twrite) Then
                        zt = arr2dc (i1, i2)
                        If (tparent) Then
                           Write (Unit=unt, Fmt=trim(str)) i1, i2, dble &
                          & (zt), aimag (zt)
                        Else
                           Write (Unit=unt, Fmt=trim(str)) i1, i2, zt
                        End If
                     Else
                        Read (Unit=unt, Fmt=*) j1, j2, zt
                        arr2dc (j1, j2) = zt
                     End If
                  End Do
               End Do
            End Select
         Case (3)
       !----------------------!
       !     3 dimensions     !
       !----------------------!
            Select Case (iot)
            Case ('l')
               lb (1:3) = lbound (arr3dl)
               ub (1:3) = ubound (arr3dl)
               str = '(3' // fmtdi // ', ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        If (twrite) Then
                           lt = arr3dl (i1, i2, i3)
                           Write (Unit=unt, Fmt=trim(str)) i1, i2, i3, &
                          & lt
                        Else
                           Read (Unit=unt, Fmt=*) j1, j2, j3, lt
                           arr3dl (j1, j2, j3) = lt
                        End If
                     End Do
                  End Do
               End Do
            Case ('i')
               lb (1:3) = lbound (arr3di)
               ub (1:3) = ubound (arr3di)
               str = '(3' // fmtdi // ', ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        If (twrite) Then
                           it = arr3di (i1, i2, i3)
                           Write (Unit=unt, Fmt=trim(str)) i1, i2, i3, &
                          & it
                        Else
                           Read (Unit=unt, Fmt=*) j1, j2, j3, it
                           arr3di (j1, j2, j3) = it
                        End If
                     End Do
                  End Do
               End Do
            Case ('r')
               lb (1:3) = lbound (arr3dr)
               ub (1:3) = ubound (arr3dr)
               str = '(3' // fmtdi // ', " ", ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        If (twrite) Then
                           rt = arr3dr (i1, i2, i3)
                           Write (Unit=unt, Fmt=trim(str)) i1, i2, i3, &
                          & rt
                        Else
                           Read (Unit=unt, Fmt=*) j1, j2, j3, rt
                           arr3dr (j1, j2, j3) = rt
                        End If
                     End Do
                  End Do
               End Do
            Case ('c')
               lb (1:3) = lbound (arr3dc)
               ub (1:3) = ubound (arr3dc)
               str = '(3' // fmtdi // ', " ", ' // trim (frmt) // ', " &
              &", ' // trim (frmt) // ')'
               If (tparent) str = '(3' // fmtdi // ', " (", ' // trim &
              & (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               If (tfmtidx) Then
                  str = '(' // fmtidx // ', " ", ' // trim (frmt) // ',&
                 & " ", ' // trim (frmt) // ')'
                  If (tparent) str = '(' // fmtidx // ', " (", ' // &
                 & trim (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               End If
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        If (twrite) Then
                           zt = arr3dc (i1, i2, i3)
                           Write (Unit=unt, Fmt=trim(str)) i1, i2, i3, &
                          & zt
                        Else
                           Read (Unit=unt, Fmt=*) j1, j2, j3, zt
                           arr3dc (j1, j2, j3) = zt
                        End If
                     End Do
                  End Do
               End Do
            End Select
         Case (4)
       !----------------------!
       !     4 dimensions     !
       !----------------------!
            Select Case (iot)
            Case ('l')
               lb (1:4) = lbound (arr4dl)
               ub (1:4) = ubound (arr4dl)
               str = '(4' // fmtdi // ', ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           If (twrite) Then
                              lt = arr4dl (i1, i2, i3, i4)
                              Write (Unit=unt, Fmt=trim(str)) i1, i2, &
                             & i3, i4, lt
                           Else
                              Read (Unit=unt, Fmt=*) j1, j2, j3, j4, lt
                              arr4dl (j1, j2, j3, j4) = lt
                           End If
                        End Do
                     End Do
                  End Do
               End Do
            Case ('i')
               lb (1:4) = lbound (arr4di)
               ub (1:4) = ubound (arr4di)
               str = '(4' // fmtdi // ', ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           If (twrite) Then
                              it = arr4di (i1, i2, i3, i4)
                              Write (Unit=unt, Fmt=trim(str)) i1, i2, &
                             & i3, i4, it
                           Else
                              Read (Unit=unt, Fmt=*) j1, j2, j3, j4, it
                              arr4di (j1, j2, j3, j4) = it
                           End If
                        End Do
                     End Do
                  End Do
               End Do
            Case ('r')
               lb (1:4) = lbound (arr4dr)
               ub (1:4) = ubound (arr4dr)
               str = '(4' // fmtdi // ', " ", ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           If (twrite) Then
                              rt = arr4dr (i1, i2, i3, i4)
                              Write (Unit=unt, Fmt=trim(str)) i1, i2, &
                             & i3, i4, rt
                           Else
                              Read (Unit=unt, Fmt=*) j1, j2, j3, j4, rt
                              arr4dr (j1, j2, j3, j4) = rt
                           End If
                        End Do
                     End Do
                  End Do
               End Do
            Case ('c')
               lb (1:4) = lbound (arr4dc)
               ub (1:4) = ubound (arr4dc)
               str = '(4' // fmtdi // ', " ", ' // trim (frmt) // ', " &
              &", ' // trim (frmt) // ')'
               If (tparent) str = '(4' // fmtdi // ', " (", ' // trim &
              & (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               If (tfmtidx) Then
                  str = '(' // fmtidx // ', " ", ' // trim (frmt) // ',&
                 & " ", ' // trim (frmt) // ')'
                  If (tparent) str = '(' // fmtidx // ', " (", ' // &
                 & trim (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               End If
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           If (twrite) Then
                              zt = arr4dc (i1, i2, i3, i4)
                              Write (Unit=unt, Fmt=trim(str)) i1, i2, &
                             & i3, i4, zt
                           Else
                              Read (Unit=unt, Fmt=*) j1, j2, j3, j4, zt
                              arr4dc (j1, j2, j3, j4) = zt
                           End If
                        End Do
                     End Do
                  End Do
               End Do
            End Select
         Case (5)
       !----------------------!
       !     5 dimensions     !
       !----------------------!
            Select Case (iot)
            Case ('l')
               lb (1:5) = lbound (arr5dl)
               ub (1:5) = ubound (arr5dl)
               str = '(5' // fmtdi // ', ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           Do i5 = lb (5), ub (5)
                              If (twrite) Then
                                 lt = arr5dl (i1, i2, i3, i4, i5)
                                 Write (Unit=unt, Fmt=trim(str)) i1, &
                                & i2, i3, i4, i5, lt
                              Else
                                 Read (Unit=unt, Fmt=*) j1, j2, j3, j4, &
                                & j5, lt
                                 arr5dl (j1, j2, j3, j4, j5) = lt
                              End If
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            Case ('i')
               lb (1:5) = lbound (arr5di)
               ub (1:5) = ubound (arr5di)
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               str = '(5' // fmtdi // ', ' // trim (frmt) // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           Do i5 = lb (5), ub (5)
                              If (twrite) Then
                                 it = arr5di (i1, i2, i3, i4, i5)
                                 Write (Unit=unt, Fmt=trim(str)) i1, &
                                & i2, i3, i4, i5, it
                              Else
                                 Read (Unit=unt, Fmt=*) j1, j2, j3, j4, &
                                & j5, it
                                 arr5di (j1, j2, j3, j4, j5) = it
                              End If
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            Case ('r')
               lb (1:5) = lbound (arr5dr)
               ub (1:5) = ubound (arr5dr)
               str = '(5' // fmtdi // ', " ", ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           Do i5 = lb (5), ub (5)
                              If (twrite) Then
                                 rt = arr5dr (i1, i2, i3, i4, i5)
                                 Write (Unit=unt, Fmt=trim(str)) i1, &
                                & i2, i3, i4, i5, rt
                              Else
                                 Read (Unit=unt, Fmt=*) j1, j2, j3, j4, &
                                & j5, rt
                                 arr5dr (j1, j2, j3, j4, j5) = rt
                              End If
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            Case ('c')
               lb (1:5) = lbound (arr5dc)
               ub (1:5) = ubound (arr5dc)
               str = '(5' // fmtdi // ', " ", ' // trim (frmt) // ', " &
              &", ' // trim (frmt) // ')'
               If (tparent) str = '(5' // fmtdi // ', " (", ' // trim &
              & (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               If (tfmtidx) Then
                  str = '(' // fmtidx // ', " ", ' // trim (frmt) // ',&
                 & " ", ' // trim (frmt) // ')'
                  If (tparent) str = '(' // fmtidx // ', " (", ' // &
                 & trim (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               End If
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           Do i5 = lb (5), ub (5)
                              If (twrite) Then
                                 zt = arr5dc (i1, i2, i3, i4, i5)
                                 Write (Unit=unt, Fmt=trim(str)) i1, &
                                & i2, i3, i4, i5, zt
                              Else
                                 Read (Unit=unt, Fmt=*) j1, j2, j3, j4, &
                                & j5, zt
                                 arr5dc (j1, j2, j3, j4, j5) = zt
                              End If
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Select
         Case (6)
       !----------------------!
       !     6 dimensions     !
       !----------------------!
            Select Case (iot)
            Case ('l')
               lb (1:6) = lbound (arr6dl)
               ub (1:6) = ubound (arr6dl)
               str = '(6' // fmtdi // ', ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           Do i5 = lb (5), ub (5)
                              Do i6 = lb (6), ub (6)
                                 If (twrite) Then
                                    lt = arr6dl (i1, i2, i3, i4, i5, &
                                   & i6)
                                    Write (Unit=unt, Fmt=trim(str)) i1, &
                                   & i2, i3, i4, i5, i6, lt
                                 Else
                                    Read (Unit=unt, Fmt=*) j1, j2, j3, &
                                   & j4, j5, j6, lt
                                    arr6dl (j1, j2, j3, j4, j5, j6) = &
                                   & lt
                                 End If
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            Case ('i')
               lb (1:6) = lbound (arr6di)
               ub (1:6) = ubound (arr6di)
               str = '(6' // fmtdi // ', ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           Do i5 = lb (5), ub (5)
                              Do i6 = lb (6), ub (6)
                                 If (twrite) Then
                                    it = arr6di (i1, i2, i3, i4, i5, &
                                   & i6)
                                    Write (Unit=unt, Fmt=trim(str)) i1, &
                                   & i2, i3, i4, i5, i6, it
                                 Else
                                    Read (Unit=unt, Fmt=*) j1, j2, j3, &
                                   & j4, j5, j6, it
                                    arr6di (j1, j2, j3, j4, j5, j6) = &
                                   & it
                                 End If
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            Case ('r')
               lb (1:6) = lbound (arr6dr)
               ub (1:6) = ubound (arr6dr)
               str = '(6' // fmtdi // ', " ", ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           Do i5 = lb (5), ub (5)
                              Do i6 = lb (6), ub (6)
                                 If (twrite) Then
                                    rt = arr6dr (i1, i2, i3, i4, i5, &
                                   & i6)
                                    Write (Unit=unt, Fmt=trim(str)) i1, &
                                   & i2, i3, i4, i5, i6, rt
                                 Else
                                    Read (Unit=unt, Fmt=*) j1, j2, j3, &
                                   & j4, j5, j6, rt
                                    arr6dr (j1, j2, j3, j4, j5, j6) = &
                                   & rt
                                 End If
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            Case ('c')
               lb (1:6) = lbound (arr6dc)
               ub (1:6) = ubound (arr6dc)
               str = '(6' // fmtdi // ', " ", ' // trim (frmt) // ', " &
              &", ' // trim (frmt) // ')'
               If (tparent) str = '(6' // fmtdi // ', " (", ' // trim &
              & (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               If (tfmtidx) Then
                  str = '(' // fmtidx // ', " ", ' // trim (frmt) // ',&
                 & " ", ' // trim (frmt) // ')'
                  If (tparent) str = '(' // fmtidx // ', " (", ' // &
                 & trim (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               End If
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           Do i5 = lb (5), ub (5)
                              Do i6 = lb (6), ub (6)
                                 If (twrite) Then
                                    zt = arr6dc (i1, i2, i3, i4, i5, &
                                   & i6)
                                    Write (Unit=unt, Fmt=trim(str)) i1, &
                                   & i2, i3, i4, i5, i6, zt
                                 Else
                                    Read (Unit=unt, Fmt=*) j1, j2, j3, &
                                   & j4, j5, j6, zt
                                    arr6dc (j1, j2, j3, j4, j5, j6) = &
                                   & zt
                                 End If
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Select
         Case (7)
       !----------------------!
       !     7 dimensions     !
       !----------------------!
            Select Case (iot)
            Case ('l')
               lb (1:7) = lbound (arr7dl)
               ub (1:7) = ubound (arr7dl)
               str = '(7' // fmtdi // ', ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           Do i5 = lb (5), ub (5)
                              Do i6 = lb (6), ub (6)
                                 Do i7 = lb (7), ub (7)
                                    If (twrite) Then
                                       lt = arr7dl (i1, i2, i3, i4, i5, &
                                      & i6, i7)
                                       Write (Unit=unt, Fmt=trim(str)) &
                                      & i1, i2, i3, i4, i5, i6, i7, lt
                                    Else
                                       Read (Unit=unt, Fmt=*) j1, j2, &
                                      & j3, j4, j5, j6, j7, lt
                                       arr7dl (j1, j2, j3, j4, j5, j6, &
                                      & j7) = lt
                                    End If
                                 End Do
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            Case ('i')
               lb (1:7) = lbound (arr7di)
               ub (1:7) = ubound (arr7di)
               str = '(7' // fmtdi // ', ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           Do i5 = lb (5), ub (5)
                              Do i6 = lb (6), ub (6)
                                 Do i7 = lb (7), ub (7)
                                    If (twrite) Then
                                       it = arr7di (i1, i2, i3, i4, i5, &
                                      & i6, i7)
                                       Write (Unit=unt, Fmt=trim(str)) &
                                      & i1, i2, i3, i4, i5, i6, i7, it
                                    Else
                                       Read (Unit=unt, Fmt=*) j1, j2, &
                                      & j3, j4, j5, j6, j7, it
                                       arr7di (j1, j2, j3, j4, j5, j6, &
                                      & j7) = it
                                    End If
                                 End Do
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            Case ('r')
               lb (1:7) = lbound (arr7dr)
               ub (1:7) = ubound (arr7dr)
               str = '(7' // fmtdi // ', " ", ' // trim (frmt) // ')'
               If (tfmtidx) str = '(' // fmtidx // ', ' // trim (frmt) &
              & // ')'
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           Do i5 = lb (5), ub (5)
                              Do i6 = lb (6), ub (6)
                                 Do i7 = lb (7), ub (7)
                                    If (twrite) Then
                                       rt = arr7dr (i1, i2, i3, i4, i5, &
                                      & i6, i7)
                                       Write (Unit=unt, Fmt=trim(str)) &
                                      & i1, i2, i3, i4, i5, i6, i7, rt
                                    Else
                                       Read (Unit=unt, Fmt=*) j1, j2, &
                                      & j3, j4, j5, j6, j7, rt
                                       arr7dr (j1, j2, j3, j4, j5, j6, &
                                      & j7) = rt
                                    End If
                                 End Do
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            Case ('c')
               lb (1:7) = lbound (arr7dc)
               ub (1:7) = ubound (arr7dc)
               str = '(7' // fmtdi // ', " ", ' // trim (frmt) // ', " &
              &", ' // trim (frmt) // ')'
               If (tparent) str = '(7' // fmtdi // ', " (", ' // trim &
              & (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               If (tfmtidx) Then
                  str = '(' // fmtidx // ', " ", ' // trim (frmt) // ',&
                 & " ", ' // trim (frmt) // ')'
                  If (tparent) str = '(' // fmtidx // ', " (", ' // &
                 & trim (frmt) // ', ", ", ' // trim (frmt) // ', ")")'
               End If
               Do i1 = lb (1), ub (1)
                  Do i2 = lb (2), ub (2)
                     Do i3 = lb (3), ub (3)
                        Do i4 = lb (4), ub (4)
                           Do i5 = lb (5), ub (5)
                              Do i6 = lb (6), ub (6)
                                 Do i7 = lb (7), ub (7)
                                    If (twrite) Then
                                       zt = arr7dc (i1, i2, i3, i4, i5, &
                                      & i6, i7)
                                       Write (Unit=unt, Fmt=trim(str)) &
                                      & i1, i2, i3, i4, i5, i6, i7, zt
                                    Else
                                       Read (Unit=unt, Fmt=*) j1, j2, &
                                      & j3, j4, j5, j6, j7, zt
                                       arr7dc (j1, j2, j3, j4, j5, j6, &
                                      & j7) = zt
                                    End If
                                 End Do
                              End Do
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Select
         End Select
      End Subroutine ioarr

  ! determine IO (read/write) action
      Logical Function ioaction (ioa)
         implicit none
         Character (*), Intent (In) :: ioa
         Select Case (trim(adjustl(ioa)))
         Case ('r', 'R', 'read', 'Read', 'READ')
            ioaction = .False.
         Case ('w', 'W', 'write', 'Write', 'WRITE')
            ioaction = .True.
         Case Default
            Write (*,*) 'Error(ioarray::ioaction): unknown IO action: ' &
           & // trim (adjustl(ioa))
            Stop
         End Select
      End Function ioaction
End Module ioarray
