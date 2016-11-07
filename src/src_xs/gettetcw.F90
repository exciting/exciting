!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_gettetcw
  use modmpi
      Implicit None
Contains
!
!
      Subroutine gettetcw (iq, ik, i1, i2, n1, n2, nw, fnam, cw, cwa, &
     & cwsurf)
         Use modmain
         Use modinput
         Use modxs
         Use m_getunit
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq, ik, i1, i2, n1, n2, nw
         Character (*), Intent (In) :: fnam
         Real (8), Intent (Out) :: cw (nw), cwa (nw), cwsurf (nw)
    ! local variables
         Integer :: un, irec, recl, nstsv_, n1_, n2_, err
         Real (8) :: vql_ (3), vkl_ (3), vqlt (3), vklt (3)
         Real (8), External :: r3dist
         err = 0
         If ((i1 .Lt. 1) .Or. (i1 .Gt. nstsv) .Or. (i2 .Lt. 1) .Or. (i2 &
        & .Gt. nstsv)) Then
            Write (unitout,*)
            Write (unitout, '("Error(gettetcw): inconsistent band combi&
           &nation:")')
            Write (unitout, '(" bands 	: ", 2i6)') i1, i2
            Write (unitout, '(" maximum value : ", i6)') nstfv
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If (err .Gt. 0) Call terminate
    !------------------------!
    !     get parameters     !
    !------------------------!
         Inquire (IoLength=Recl) vql_, vkl_, nstsv_, n1_, n2_
         Call getunit (un)
         Open (un, File=trim(fnam), Action='read', Form='unformatted', &
        & Status='old', Access='direct', Recl=Recl)
         Read (un, Rec=1) vql_, vkl_, nstsv_, n1_, n2_
         Close (un)
         err = 0
    ! check number of bands
         If (nstsv .Gt. nstsv_) Then
            Write (unitout,*)
            Write (unitout, '("Error(gettetcw): invalid nstsv for k-poi&
           &nt ", I8)') ik
            Write (unitout, '(" q-point    : ", I8)') iq
            Write (unitout, '(" current    : ", I8)') nstsv
            Write (unitout, '(" FILE	     : ", I8)') nstsv_
            Write (unitout, '(" filename   : ", a )') trim (fnam)
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If ((n1 .Ne. n1_) .Or. (n2 .Gt. n2_)) Then
            Write (unitout,*)
            Write (unitout, '("Error(gettetcw): invalid band ranges for&
           & k-point ", I8)') ik
            Write (unitout, '(" q-point    : ", I8)') iq
            Write (unitout, '(" current    : ", 2I8)') n1, n2
            Write (unitout, '(" FILE	     : ", 2I8)') n1_, n2_
            Write (unitout, '(" filename   : ", a )') trim (fnam)
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If (err .Gt. 0) Call terminate
    !------------------!
    !     get data     !
    !------------------!
    ! record position with proper n1 and n2 values
         irec = (ik-1) * n1_ * n2_ + (i1-1) * n2_ + i2
    ! read from file
         Call getunit (un)
         Inquire (IoLength=Recl) vql_, vkl_, nstsv_, n1_, n2_, cw, cwa, &
        & cwsurf
         Open (un, File=trim(fnam), Form='unformatted', Action='read', &
        & Status='old', Access='direct', Recl=Recl)
         Read (un, Rec=irec) vql_, vkl_, nstsv_, n1_, n2_, cw, cwa, &
        & cwsurf
         Close (un)
    ! check q-point and k-point
         If (iq .Eq. 0) Then
       ! Gamma Q-point
            vklt (:) = vkl0 (:, ik)
            vqlt (:) = 0.d0
         Else
            vklt (:) = vkl (:, ik)
            vqlt (:) = vql (:, iq)
         End If
         If ((r3dist(vkl_, vklt) .Gt. input%structure%epslat) .Or. &
        & (r3dist(vql_, vqlt) .Gt. input%structure%epslat)) Then
            Write (unitout,*)
            Write (unitout, '(a)') 'Error(gettetcw): differring paramet&
           &ers for tetrahedron convolution weights (current/file): '
            Write (unitout, '(a, i6)') ' q-point index  :', iq
            Write (unitout, '(a, i6)') ' k-point index  :', ik
            Write (unitout, '(a, 3f12.6, a, 3f12.6)') ' vql		 :', vqlt, &
           & ', ', vql_
            Write (unitout, '(a, 3f12.6, a, 3f12.6)') ' vkl		 :', vklt, &
           & ', ', vkl_
            Write (unitout, '(a)') ' file		 : ' // trim (fnam)
            Write (unitout,*)
            Call flushifc (unitout)
            Call terminate
         End If
      End Subroutine gettetcw
!
End Module m_gettetcw
