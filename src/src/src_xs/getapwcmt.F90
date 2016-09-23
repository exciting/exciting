!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_getapwcmt
      Implicit None
Contains
!
  ! APW functions
!
!
      Subroutine getapwcmt (iq, ik, isti, istf, lmax, apwlm)
         Use modmain
         Use modinput
         Use modxs
         Use m_getunit
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq, ik, isti, istf, lmax
         Complex (8), Intent (Out) :: apwlm (:, :, :, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'getapwcmt'
         Character (256) :: filextt
         Integer :: un, recl, err, nstfv_, apwordmax_, lmaxapw_
         Real (8) :: vql_ (3), vkl_ (3), vklt (3), vqlt (3)
         Complex (8), Allocatable :: apwlmt (:, :, :, :)
         Real (8), External :: r3dist
         err = 0
    ! check band range
         If ((isti .Lt. 1) .Or. (istf .Gt. nstfv) .Or. (istf .Le. &
        & isti)) Then
            Write (unitout,*)
            Write (unitout, '("Error(getapwcmt): inconsistent limits fo&
           &r bands:")')
            Write (unitout, '(" band limits	: ", 2i6)') isti, istf
            Write (unitout, '(" maximum value : ", i6)') nstfv
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If (size(apwlm, 1) .Ne. (istf-isti+1)) Then
            Write (unitout,*)
            Write (unitout, '("Error(getapwcmt): output array does not &
           &match for bands:")')
            Write (unitout, '(" band limits		    : ", 2i6)') isti, istf
            Write (unitout, '(" requested number of bands : ", i6)') &
           & istf - isti + 1
            Write (unitout, '(" array size		    : ", i6)') size (apwlm, &
           & 1)
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
    ! check lmax value
         If ((lmax .Gt. input%groundstate%lmaxapw) .Or. (lmax .Lt. 0)) &
        & Then
            Write (unitout,*)
            Write (unitout, '(a, i8)') 'Error(' // thisnam // '): lmax &
           &> input%groundstate%lmaxapw or < 0:', lmax
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If (err .Gt. 0) Call terminate
    ! set file extension
         filextt = filext
         If (iq .Eq. 0) Call genfilextread (task)
    !------------------------!
    !     get parameters     !
    !------------------------!
         Inquire (IoLength=Recl) vql_, vkl_, nstfv_, apwordmax_, &
        & lmaxapw_
         Call getunit (un)
         Open (un, File='APWCMT'//trim(filext), Action='read', Form='un&
        &formatted', Status='old', Access='direct', Recl=Recl)
         Read (un, Rec=1) vql_, vkl_, nstfv_, apwordmax_, lmaxapw_
         Close (un)
         err = 0
    ! check number of bands
         If (nstfv .Gt. nstfv_) Then
            Write (unitout,*)
            Write (unitout, '("Error(", a, "): invalid nstfv for k-poin&
           &t ", I8)') thisnam, ik
            Write (unitout, '(" q-point    : ", I8)') iq
            Write (unitout, '(" current    : ", I8)') nstfv
            Write (unitout, '(" FILE	     : ", I8)') nstfv_
            Write (unitout, '(" filename   : ", a )') 'APWCMT' // trim &
           & (filext)
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
    ! check APW matching order
         If (apwordmax .Ne. apwordmax_) Then
            Write (unitout,*)
            Write (unitout, '("Error(", a, "): invalid apwordmax for k-&
           &point ", I8)') thisnam, ik
            Write (unitout, '(" q-point    : ", I8)') iq
            Write (unitout, '(" current    : ", I8)') apwordmax
            Write (unitout, '(" FILE	     : ", I8)') apwordmax_
            Write (unitout, '(" filename   : ", a )') 'APWCMT' // trim &
           & (filext)
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
    ! check lmax
         If (input%groundstate%lmaxapw .Gt. lmaxapw_) Then
            Write (unitout,*)
            Write (unitout, '("Error(", a, "): invalid lmaxapw for k-po&
           &int ", I8)') thisnam, ik
            Write (unitout, '(" q-point    : ", I8)') iq
            Write (unitout, '(" current    : ", I8)') &
           & input%groundstate%lmaxapw
            Write (unitout, '(" FILE	     : ", I8)') lmaxapw_
            Write (unitout, '(" filename   : ", a )') 'APWCMT' // trim &
           & (filext)
            Call flushifc (unitout)
            Write (unitout,*)
            err = err + 1
         End If
         If (err .Gt. 0) Call terminate
    !------------------!
    !     get data     !
    !------------------!
    ! assign to output array and apply cutoff
         Allocate (apwlmt(nstfv_, apwordmax, (lmaxapw_+1)**2, natmtot))
    ! read data from file
         Inquire (IoLength=Recl) vql_, vkl_, nstfv_, apwordmax_, &
        & lmaxapw_, apwlmt
         Call getunit (un)
         Open (un, File='APWCMT'//trim(filext), Action='read', Form='un&
        &formatted', Status='old', Access='direct', Recl=Recl)
         Read (un, Rec=ik) vql_, vkl_, nstfv_, apwordmax_, lmaxapw_, &
        & apwlmt
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
        & ((r3dist(vql_, vqlt) .Gt. input%structure%epslat) .And. ( &
        & .Not. tscreen))) Then
            Write (unitout,*)
            Write (unitout, '(a)') 'Error(' // thisnam // '): differrin&
           &g parameters for APW MT coefficients (current/file): '
            Write (unitout, '(a, i6)') ' q-point index  :', iq
            Write (unitout, '(a, i6)') ' k-point index  :', ik
            Write (unitout, '(a, 3f12.6, a, 3f12.6)') ' vql		 :', vqlt, &
           & ', ', vql_
            Write (unitout, '(a, 3f12.6, a, 3f12.6)') ' vkl		 :', vklt, &
           & ', ', vkl_
            Write (unitout, '(a)') ' file		 : APWCMT' // trim (filext)
            Write (unitout,*)
            Call flushifc (unitout)
            Call terminate
         End If
    ! retrieve data within cutoffs
         apwlm (:, :, :, :) = apwlmt (isti:istf, :, 1:(lmax+1)**2, :)
         Deallocate (apwlmt)
    ! restore file extension
         filext = filextt
      End Subroutine getapwcmt
!
End Module m_getapwcmt
