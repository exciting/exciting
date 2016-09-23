!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_puttetcw
      Implicit None
Contains
!
!
      Subroutine puttetcw (iq, ik, i1, i2, n1, n2, filnam, cw, cwa, &
     & cwsurf)
         Use modmain
         Use modxs
         Use modmpi
         Use m_getunit
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq, ik, i1, i2, n1, n2
         Character (*), Intent (In) :: filnam
         Real (8), Intent (In) :: cw (:), cwa (:), cwsurf (:)
    ! local variables
         Integer :: un, recl, irec, iqt, err
    ! check input parameters
         err = 0
         If ((n1 .Lt. 1) .Or. (n1 .Gt. nstsv)) Then
            Write (unitout,*)
            Write (unitout, '("Error(puttetcw): n1 < 1 or n1 > nstsv")')
            Write (unitout, '(" n1	:", i6)') n1
            Write (unitout, '(" nstsv :", i6)') nstsv
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If ((n2 .Lt. 1) .Or. (n2 .Gt. nstsv)) Then
            Write (unitout,*)
            Write (unitout, '("Error(puttetcw): n2 < 1 or n2 > nstsv")')
            Write (unitout, '(" n1	:", i6)') n2
            Write (unitout, '(" nstsv :", i6)') nstsv
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If (err .Ne. 0) Call terminate
    ! q-point
         iqt = iq
    ! record position
         irec = (ik-1) * n1 * n2 + (i1-1) * n2 + i2
    ! I/O record length
         Inquire (IoLength=Recl) vql (:, iq), vkl (:, ik), nstsv, n1, &
        & n2, cw, cwa, cwsurf
         Call getunit (un)
         Open (Unit=un, File=trim(filnam), Form='unformatted', Action='&
        &write', Access='direct', Recl=Recl)
         Write (un, Rec=irec) vql (:, iq), vkl (:, ik), nstsv, n1, n2, &
        & cw, cwa, cwsurf
         Close (un)
      End Subroutine puttetcw
!
End Module m_puttetcw
