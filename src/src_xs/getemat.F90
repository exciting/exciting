!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_getemat
      Implicit None
Contains
!
!
      Subroutine getemat (iq, ik, tarec, filnam, ngp, l1, h1, l2, h2, &
     & x12, l3, h3, l4, h4, x34)
         Use modmain
         Use modinput
         Use modxs
         Use modmpi
         Use m_getunit
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq, ik, ngp
         Logical, intent(in) :: tarec
         Character (*), intent(in) :: filnam
         Integer, Intent (In) :: l1, h1, l2, h2
         Complex (8), Intent (Out) :: x12 (:, :, :)
         Integer, Optional, Intent (In) :: l3, h3, l4, h4
         Complex (8), Optional, Intent (Out) :: x34 (:, :, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'getemat'
         Integer :: recl, un, ikr, n1, n2, n3, n4, nstsv_, ngq_, err
         Integer :: l1_, h1_, l2_, h2_, l3_, h3_, l4_, h4_, n1_, n2_, &
        & n3_, n4_
         Real (8) :: vql_ (3), vkl_ (3)
         Logical :: existent, lerr
         Complex (8), Allocatable :: x12t (:, :, :), x34t (:, :, :)
    ! functions
         Real (8) :: r3dist
         External :: r3dist
    ! check if all optional variables are present if any is present
         If ((present(l3) .Or. present(h3) .Or. present(l4) .Or. &
        & present(h4) .Or. present(x34)) .And. ( .Not. (present(l3) &
        & .And. present(h3) .And. present(l4) .And. present(h4) .And. &
        & present(x34)))) Then
            Write (*,*)
            Write (*, '("Error(getemat): optional parameters not comple&
           &te - check calling routines")')
            Write (*,*)
            Call terminate
         End If
    ! check if file exists
         Inquire (File=trim(filnam), Exist=existent)
         If ( .Not. existent) Then
            Write (unitout, '(a)') 'Error(' // thisnam // '): file does&
           & not exist: ' // trim (filnam)
            Call terminate
         End If
    ! record position for k-point
         ikr = ik
         If ( .Not. tarec) Call getridx(nkpt, ik, ikr)
    ! check limits for states
         lerr = (l1 .Lt. 1) .Or. (l1 .Gt. nstsv) .Or. (h1 .Lt. 1) .Or. &
        & (h1 .Gt. nstsv) .Or. (l2 .Lt. 1) .Or. (l2 .Gt. nstsv) .Or. &
        & (h2 .Lt. 1) .Or. (h2 .Gt. nstsv) .Or. (l1 .Gt. h1) .Or. (l2 &
        & .Gt. h2)
         If (present(x34)) lerr = lerr .Or. (l3 .Lt. 1) .Or. (l3 .Gt. &
        & nstsv) .Or. (h3 .Lt. 1) .Or. (h3 .Gt. nstsv) .Or. (l4 .Lt. 1) &
        & .Or. (l4 .Gt. nstsv) .Or. (h4 .Lt. 1) .Or. (h4 .Gt. nstsv) &
        & .Or. (l3 .Gt. h3) .Or. (l4 .Gt. h4)
         err = 0
         If (lerr) Then
            Write (unitout,*)
            Write (unitout, '("Error(", a, "): inconsistent requested l&
           &imits for states:")') thisnam
            If (present(x34)) Then
               Write (unitout, '(" requested state limits (lo, hi): ", &
              &4(2i6, 2x))') l1, h1, l2, h2, l3, h3, l4, h4
            Else
               Write (unitout, '(" requested state limits (lo, hi): ", &
              &2(2i6, 2x))') l1, h1, l2, h2
            End If
            Write (unitout, '(" maximum value 		: ", i6)') nstsv
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         n1 = h1 - l1 + 1
         n2 = h2 - l2 + 1
         If (present(x34)) Then
            n3 = h3 - l3 + 1
            n4 = h4 - l4 + 1
         End If
    ! check block sizes against array
         lerr = (size(x12, 1) .Ne. n1) .Or. (size(x12, 2) .Ne. n2) .Or. &
        & (ngp .Gt. ngq(iq)) .Or. (size(x12, 3) .Ne. ngp)
         If (present(x34)) lerr = lerr .Or. (size(x34, 1) .Ne. n3) .Or. &
        & (size(x34, 2) .Ne. n4) .Or. (size(x34, 3) .Ne. ngp)
         If (lerr) Then
            Write (unitout,*)
            Write (unitout, '("Error(", a, "): output array does not ma&
           &tch for states:")') thisnam
            Write (unitout, '(" requested number of G+q vectors : ", i6&
           &)') ngp
            Write (unitout, '(" current number of G+q vectors   : ", i6&
           &)') ngq (iq)
            If (present(x34)) Then
               Write (unitout, '(" array sizes for G+q vectors     : ",&
              & 2i6)') size (x12, 3), size (x34, 3)
               Write (unitout, '(" block sizes : ", 4i6)') n1, n2, n3, &
              & n4
               Write (unitout, '(" array sizes : ", 4i6)') size (x12, &
              & 1), size (x12, 2), size (x34, 1), size (x34, 2)
            Else
               Write (unitout, '(" array size for G+q vectors      : ",&
              & i6)') size (x12, 3)
               Write (unitout, '(" block sizes : ", 2i6)') n1, n2
               Write (unitout, '(" array sizes : ", 2i6)') size (x12, &
              & 1), size (x12, 2)
            End If
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If (err .Gt. 0) Call terminate
    !------------------------!
    !     get parameters     !
    !------------------------!
         Call getunit (un)
         If (present(x34)) Then
       ! I/O record length
            Inquire (IoLength=Recl) vql_, vkl_, nstsv_, ngq_, l1_, h1_, &
           & l2_, h2_, l3_, h3_, l4_, h4_
            Open (Unit=un, File=trim(filnam), Form='unformatted', &
           & Action='read', Access='direct', Recl=Recl)
            Read (un, Rec=1) vql_, vkl_, nstsv_, ngq_, l1_, h1_, l2_, &
           & h2_, l3_, h3_, l4_, h4_
            Close (un)
         Else
       ! I/O record length
            Inquire (IoLength=Recl) vql_, vkl_, nstsv_, ngq_, l1_, h1_, &
           & l2_, h2_
            Open (Unit=un, File=trim(filnam), Form='unformatted', &
           & Action='read', Access='direct', Recl=Recl)
            Read (un, Rec=1) vql_, vkl_, nstsv_, ngq_, l1_, h1_, l2_, &
           & h2_
            Close (un)
         End If
         err = 0
    ! check block sizes
         lerr = (l1 .Lt. l1_) .Or. (h1 .Gt. h1_) .Or. (l2 .Lt. l2_) &
        & .Or. (h2 .Gt. h2_)
         If (present(x34)) lerr = lerr .Or. (l3 .Lt. l3_) .Or. (h3 .Gt. &
        & h3_) .Or. (l4 .Lt. l4_) .Or. (h4 .Gt. h4_)
         If (lerr) Then
            Write (unitout,*)
            Write (unitout, '("Error(", a, "): limits for states out of&
           & range in file:")') thisnam
            If (present(x34)) Then
               Write (unitout, '(" requested state limits (lo, hi): ", &
              &4(2i6, 2x))') l1, h1, l2, h2, l3, h3, l4, h4
               Write (unitout, '(" state limits from file (lo, hi): ", &
              &4(2i6, 2x))') l1_, h1_, l2_, h2_, l3_, h3_, l4_, h4_
            Else
               Write (unitout, '(" requested state limits (lo, hi): ", &
              &2(2i6, 2x))') l1, h1, l2, h2
               Write (unitout, '(" state limits from file (lo, hi): ", &
              &2(2i6, 2x))') l1_, h1_, l2_, h2_
               Write (unitout, '(" file			   : ", a)') trim (filnam)
            End If
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If (err .Gt. 0) Call terminate
    !------------------!
    !     get data     !
    !------------------!
         n1_ = h1_ - l1_ + 1
         n2_ = h2_ - l2_ + 1
         If (present(x34)) Then
            n3_ = h3_ - l3_ + 1
            n4_ = h4_ - l4_ + 1
         End If
         Allocate (x12t(n1_, n2_, ngq_))
         If (present(x34)) allocate (x34t(n3_, n4_, ngq_))
         Call getunit (un)
         If (present(x34)) Then
       ! I/O record length
            Inquire (IoLength=Recl) vql_, vkl_, nstsv_, ngq_, l1_, h1_, &
           & l2_, h2_, l3_, h3_, l4_, h4_, x12t, x34t
            Open (Unit=un, File=trim(filnam), Form='unformatted', &
           & Action='read', Access='direct', Recl=Recl)
            Read (un, Rec=ikr) vql_, vkl_, nstsv_, ngq_, l1_, h1_, l2_, &
           & h2_, l3_, h3_, l4_, h4_, x12t, x34t
         Else
       ! I/O record length
            Inquire (IoLength=Recl) vql_, vkl_, nstsv_, ngq_, l1_, h1_, &
           & l2_, h2_, x12t
            Open (Unit=un, File=trim(filnam), Form='unformatted', &
           & Action='read', Access='direct', Recl=Recl)
            Read (un, Rec=ikr) vql_, vkl_, nstsv_, ngq_, l1_, h1_, l2_, &
           & h2_, x12t
         End If
         Close (un)
    ! check q-point and k-point
         If ((r3dist(vql_, vql(1, iq)) .Gt. input%structure%epslat) &
        & .Or. (r3dist(vkl_, vkl(1, ik)) .Gt. input%structure%epslat)) &
        & Then
            Write (unitout,*)
            Write (unitout, '(a)') 'Error(' // thisnam // '): differrin&
           &g parameters for matrix elements (current/file): '
            Write (unitout, '(a, 3f12.6, a, 3f12.6)') ' vql :', vql (:, &
           & iq), ', ', vql_
            Write (unitout, '(a, 3f12.6, a, 3f12.6)') ' vkl :', vkl (:, &
           & ik), ', ', vkl_
            Write (unitout, '(a)') ' file: ', trim (filnam)
            Write (unitout,*)
            Call flushifc (unitout)
            Call terminate
         End If
    ! retrieve data within cutoff
         x12 (:, :, :) = x12t (l1-l1_+1:h1-l1_+1, l2-l2_+1:h2-l2_+1, &
        & :ngp)
         Deallocate (x12t)
         If (present(x34)) x34 (:, :, :) = x34t (l3-l3_+1:h3-l3_+1, &
        & l4-l4_+1:h4-l4_+1, :ngp)
         If (present(x34)) deallocate (x34t)
      End Subroutine getemat
!
End Module m_getemat
