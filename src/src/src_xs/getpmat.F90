!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_getpmat
      Implicit None
Contains
!
!
      Subroutine getpmat (ik, vklt, i1, f1, i2, f2, tarec, filnam, pm)
         Use modmain
         Use modinput
         Use modxs
         Use modmpi
         Use m_getunit
         Implicit None
    ! arguments
         Integer, Intent (In) :: ik, i1, f1, i2, f2
         Real (8), Intent (In) :: vklt (:, :)
         Logical, intent(in) :: tarec
         Character (*), intent(in) :: filnam
         Complex (8), Intent (Out) :: pm (:, :, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'getpmat'
         Integer :: recl, un, ikr, nstsv_, err
         Real (8) :: vkl_ (3)
         Logical :: existent
         Complex (8), Allocatable :: pmt (:, :, :)
    ! functions
         Real (8), External :: r3dist
    ! check if file exists
         Inquire (File=trim(filnam), Exist=existent)
         If ( .Not. existent) Then
            Write (unitout, '(a)') 'Error(' // thisnam // '): file does&
           & not exist: ' // trim (filnam)
            Call terminate
         End If
    ! record position for k-point
         ikr = ik
         If ( .Not. tarec) Call getridx (procs, nkpt, ik, ikr)
         err = 0
    ! check band range
         If ((i1 .Lt. 1) .Or. (i1 .Gt. nstsv) .Or. (f1 .Lt. 1) .Or. (f1 &
        & .Gt. nstsv) .Or. (i2 .Lt. 1) .Or. (i2 .Gt. nstsv) .Or. (f2 &
        & .Lt. 1) .Or. (f2 .Gt. nstsv) .Or. (i1 .Gt. f1) .Or. (i2 .Gt. &
        & f2)) Then
            Write (unitout,*)
            Write (unitout, '("Error(", a, "): inconsistent limits for &
           &states:")') thisnam
            Write (unitout, '(" limits (lo/hi) : ", 2(2i6, 2x))') i1, &
           & f1, i2, f2
            Write (unitout, '(" maximum value  : ", i6)') nstsv
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If ((size(pm, 2) .Ne. (f1-i1+1)) .Or. (size(pm, 3) .Ne. &
        & (f2-i2+1))) Then
            Write (unitout,*)
            Write (unitout, '("Error(", a, "): output array does not ma&
           &tch for states:")') thisnam
            Write (unitout, '(" limits		     : ", 2(2i6, 2x))') i1, f1, &
           & i2, f2
            Write (unitout, '(" requested number of states : ", 2i6)') &
           & f1 - i1 + 1, f2 - i2 + 1
            Write (unitout, '(" array sizes		     : ", 2i6)') size (pm, &
           & 2), size (pm, 3)
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If (err .Gt. 0) Call terminate
    !------------------------!
    !     get parameters     !
    !------------------------!
         Inquire (IoLength=Recl) vkl_, nstsv_
         Call getunit (un)
         Open (Unit=un, File=trim(filnam), Form='unformatted', Action='&
        &read', Access='direct', Recl=Recl)
         Read (un, Rec=1) vkl_, nstsv_
         Close (un)
         err = 0
    ! check if all states can be read from file
         If ((f1 .Gt. nstsv_) .Or. (f2 .Gt. nstsv_)) Then
            Write (unitout,*)
            Write (unitout, '("Error(", a, "): requested states out of &
           &range for k-point ", I8)') thisnam, ik
            Write (unitout, '(" limits	  : ", 2(2I6, 2x))') i1, f1, i2, &
           & f2
            Write (unitout, '(" cutoff from file: ", I8)') nstsv_
            Write (unitout, '(" filename	  : ", a )') trim (filnam)
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If (err .Gt. 0) Call terminate
    !------------------!
    !     get data     !
    !------------------!
    ! allocate local arrays
         Allocate (pmt(3, nstsv_, nstsv_))
    ! I/O record length
         Inquire (IoLength=Recl) vkl_, nstsv_, pmt
         Call getunit (un)
         Open (Unit=un, File=trim(filnam), Form='unformatted', Action='&
        &read', Access='direct', Recl=Recl)
    ! read from file
         Read (un, Rec=ikr) vkl_, nstsv_, pmt
         Close (un)
    ! check k-point
         If (r3dist(vkl_, vklt(1, ik)) .Gt. input%structure%epslat) &
        & Then
            Write (unitout,*)
            Write (unitout, '(a)') 'Error(' // thisnam // '): differrin&
           &g parameters for matrix elements (current/file): '
            Write (unitout, '(a, i6)') ' k-point index  :', ik
            Write (unitout, '(a, i6)') ' record position:', ikr
            Write (unitout, '(a, 3f12.6, a, 3f12.6)') ' vkl		 :', vklt &
           & (:, ik), ', ', vkl_
            Write (unitout, '(" filename	  : ", a )') trim (filnam)
            Write (unitout,*)
            Call flushifc (unitout)
            Call terminate
         End If
    ! retrieve data within cutoffs
         pm (:, :, :) = pmt (:, i1:f1, i2:f2)
         Deallocate (pmt)
      End Subroutine getpmat
!
End Module m_getpmat
