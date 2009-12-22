!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_putx0
      Implicit None
Contains
!
!
      Subroutine putx0 (tp0, iq, iw, filnam, filxt, ch0, ch0wg, ch0hd)
         Use modmain
         Use modxs
         Use m_getunit
         Implicit None
    ! arguments
         Logical, Intent (In) :: tp0
         Integer, Intent (In) :: iq, iw
         Character (*), Intent (In) :: filnam, filxt
         Complex (8), Intent (In) :: ch0 (:, :)
         Complex (8), Intent (In), Optional :: ch0wg (:, :, :), ch0hd &
        & (:, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'putx0'
         Integer :: un, recl
    ! q=0 but head or wings missing
         If (tp0 .And. (( .Not. present(ch0wg)) .Or. ( .Not. &
        & present(ch0wg)))) Then
            Write (*,*) 'Error(' // trim (thisnam) // '): q=0 but head &
           &or wings missing'
            Call terminate
         End If
         Call getunit (un)
         If (tp0) Then
       ! I/O record length
            Inquire (IoLength=Recl) ngq (iq), vql (:, iq), ch0, ch0wg, &
           & ch0hd
            Open (Unit=un, File=trim(filnam)//trim(filxt), Form='unform&
           &atted', Action='write', Access='direct', Recl=Recl)
            Write (un, Rec=iw) ngq (iq), vql (:, iq), ch0, ch0wg, ch0hd
         Else
       ! I/O record length
            Inquire (IoLength=Recl) ngq (iq), vql (:, iq), ch0
            Open (Unit=un, File=trim(filnam)//trim(filxt), Form='unform&
           &atted', Action='write', Access='direct', Recl=Recl)
            Write (un, Rec=iw) ngq (iq), vql (:, iq), ch0
         End If
         Close (un)
      End Subroutine putx0
!
End Module m_putx0
