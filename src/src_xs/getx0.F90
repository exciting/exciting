!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_getx0
      Implicit None
Contains
!
!
      Subroutine getx0 (tp0, iq, iw, filnam, filxt, ch0, ch0wg, ch0hd)
         Use modmain
         use modmpi
         Use modxs
         use mod_large_io, only: inquire_large, open_direct_unformatted_large
         use precision, only: dp, i32, long_int
         Implicit None
    ! arguments
         Logical, Intent (In) :: tp0
         integer(i32), Intent (In) :: iq, iw
         Character (*), Intent (In) :: filnam, filxt
         complex(dp), Intent (Out) :: ch0 (:, :)
         complex(dp), Intent (Out), Optional :: ch0wg (:, :, :), ch0hd &
        & (:, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'getx0'
         integer(i32) :: un, ngq_
         integer(long_int) :: recl
         real(dp) :: vql_ (3)
         Logical :: existent
    ! check if file exists
         Inquire (File=trim(filnam)//trim(filxt), Exist=existent)
         If ( .Not. existent) Then
            Write (unitout, '(a)') 'Error(' // trim (thisnam) // '): fi&
           &le does not exist:' // trim (filnam) // trim (filxt)
            Call terminate
         End If
    ! q=0 but head or wings missing
         If (tp0 .And. (( .Not. present(ch0wg)) .Or. ( .Not. &
        & present(ch0wg)))) Then
            Write (*,*) 'Error(' // trim (thisnam) // '): q=0 but head &
           &or wings missing'
            Call terminate
         End If
    ! I/O record length
         If (tp0) Then
            call inquire_large( recl, [ngq(iq)], vql(:, iq), ch0, ch0wg, ch0hd )
            call open_direct_unformatted_large( un, trim(filnam)//trim(filxt), "read", recl, "old" )
            Read (un, Rec=iw) ngq_, vql_, ch0, ch0wg, ch0hd
         Else
            call inquire_large( recl, [ngq(iq)], vql(:, iq), ch0 )
            call open_direct_unformatted_large( un, trim(filnam)//trim(filxt), "read", recl, "old" )
            Read (un, Rec=iw) ngq_, vql_, ch0
         End If
         Close (un)
         If ((ngq_ .Ne. ngq(iq)) .Or. (any(vql_ .Ne. vql(:, iq)))) Then
            Write (unitout, '(a)') 'Error(' // trim (thisnam) // '): di&
           &fferring parameters for matrix elements (current/file): '
            Write (unitout, '(a, 2i6)') 'ngq', ngq (iq), ngq_
            Write (unitout, '(a, 3f12.6, a, 3f12.6)') 'vql', vql (:, &
           & iq), ', ', vql_
            Write (unitout, '(a, i6)') 'for q-point :', iq
            Write (unitout, '(a, i6)') 'for w-point :', iw
            Write (unitout, '(a)') ' file: ', trim (filnam) // trim &
           & (filxt)
            Call terminate
         End If
      End Subroutine getx0
!
End Module m_getx0
