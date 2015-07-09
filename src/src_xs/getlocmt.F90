!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_getlocmt
      Implicit None
Contains
!
  ! local orbitals functions
!
!
      Subroutine getlocmt (iq, ik, isti, istf, lolm)
         Use modmain
         Use modinput
         Use modxs
         Use m_getunit
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq, ik, isti, istf
         Complex (8), Intent (Out) :: lolm (:, :, :, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'getlocmt'
         Character (256) :: filextt
         Integer :: un, recl, err, nstfv_, nlomax_, lolmax_
         Real (8) :: vql_ (3), vkl_ (3), vklt (3), vqlt (3)
         Complex (8), Allocatable :: lolmt (:, :, :, :)
         Real (8), External :: r3dist
         err = 0
    ! check band range
         If ((isti .Lt. 1) .Or. (istf .Gt. nstfv) .Or. (istf .Le. &
        & isti)) Then
            Write (unitout,*)
            Write (unitout, '("Error(getlocmt): inconsistent limits for&
           & bands:")')
            Write (unitout, '(" band limits  : ", 2i6)') isti, istf
            Write (unitout, '(" maximum value: ", i6)') nstfv
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If (size(lolm, 1) .Ne. (istf-isti+1)) Then
            Write (unitout,*)
            Write (unitout, '("Error(getlocmt): output array does not m&
           &atch for bands:")')
            Write (unitout, '(" band limits		   : ", 2i6)') isti, istf
            Write (unitout, '(" requested number of bands: ", i6)') &
           & istf - isti + 1
            Write (unitout, '(" array size		   : ", i6)') size (lolm, &
           & 1)
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
         Inquire (IoLength=Recl) vql_, vkl_, nstfv_, nlomax_, lolmax_
         Call getunit (un)
!write(*,*) 'LOCMT'//trim(filext)
         Open (un, File='LOCMT'//trim(filext), Action='read', Form='unf&
        &ormatted', Status='old', Access='direct', Recl=Recl)
         Read (un, Rec=1) vql_, vkl_, nstfv_, nlomax_, lolmax
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
            Write (unitout, '(" filename   : ", a )') 'LOCMT' // trim &
           & (filext)
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
    ! check number of local orbitals
         If (nlomax .Ne. nlomax_) Then
            Write (unitout,*)
            Write (unitout, '("Error(", a, "): invalid nlomax for k-poi&
           &nt ", I8)') thisnam, ik
            Write (unitout, '(" q-point    : ", I8)') iq
            Write (unitout, '(" current    : ", I8)') nlomax
            Write (unitout, '(" FILE	     : ", I8)') nlomax_
            Write (unitout, '(" filename   : ", a )') 'LOCMT' // trim &
           & (filext)
            Write (unitout,*)
            Call flushifc (unitout)
            err = err + 1
         End If
         If (err .Gt. 0) Call terminate
    !------------------!
    !     get data     !
    !------------------!
    ! assign to output array and apply cutoff
         Allocate (lolmt(nstfv_, nlomax,-lolmax:lolmax, natmtot))
    ! read data from file
         Inquire (IoLength=Recl) vql_, vkl_, nstfv_, nlomax_, lolmax_, &
        & lolmt
         Call getunit (un)
         Open (un, File='LOCMT'//trim(filext), Action='read', Form='unf&
        &ormatted', Status='old', Access='direct', Recl=Recl)
         Read (un, Rec=ik) vql_, vkl_, nstfv_, nlomax_, lolmax_, lolmt
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
           &g parameters for LO MT coefficients (current/file): '
            Write (unitout, '(a, i6)') ' q-point index  :', iq
            Write (unitout, '(a, i6)') ' k-point index  :', ik
            Write (unitout, '(a, 3f12.6, a, 3f12.6)') ' vql		 :', vqlt, &
           & ', ', vql_
            Write (unitout, '(a, 3f12.6, a, 3f12.6)') ' vkl		 :', vklt, &
           & ', ', vkl_
            Write (unitout, '(a)') ' file		 : LOCMT' // trim (filext)
            Write (unitout,*)
            Call flushifc (unitout)
            Call terminate
         End If
    ! retrieve data within cutoffs
         lolm (:, :, :, :) = lolmt (isti:istf, :, :, :)
         Deallocate (lolmt)
    ! restore file extension
         filext = filextt
      End Subroutine getlocmt
!
End Module m_getlocmt
