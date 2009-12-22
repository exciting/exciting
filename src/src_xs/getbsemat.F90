!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine getbsemat (fname, ikkp, n1, n2, zmat)
      Use m_getunit
      Implicit None
  ! arguments
      Character (*), Intent (In) :: fname
      Integer, Intent (In) :: ikkp, n1, n2
      Complex (8), Intent (Out) :: zmat (n1, n2, n1, n2)
  ! local variables
      Integer :: un, recl, ikkp_, iknr_, jknr_, iq_, iqr_, n1_, n2_, &
     & n3_, n4_
      Complex (8), Allocatable :: zm (:, :, :, :)
      Call getunit (un)
  ! get sizes of stored matrix
      Inquire (IoLength=Recl) ikkp_, iknr_, jknr_, iq_, iqr_, n1_, n2_, &
     & n3_, n4_
      Open (Unit=un, File=trim(fname), Form='unformatted', Action='read&
     &', Access='direct', Recl=Recl)
      Read (un, Rec=1) ikkp_, iknr_, jknr_, iq_, iqr_, n1_, n2_, n3_, &
     & n4_
      Close (un)
  ! check if requested size can be retrieved
      If ((n1 .Gt. n1_) .Or. (n2 .Gt. n2_)) Then
         Write (*,*)
         Write (*, '("Error(getbsemat): requested matrix size out of ra&
        &nge")')
         Write (*, '(" requested size : ", 2i8)') n1, n2
         Write (*, '(" stored size	 : ", 2i8)') n1_, n2_
         Write (*,*)
         Call terminate
      End If
  ! read matrix
      Allocate (zm(n1_, n2_, n3_, n4_))
      Inquire (IoLength=Recl) ikkp_, iknr_, jknr_, iq_, iqr_, n1_, n2_, &
     & n3_, n4_, zm
      Open (Unit=un, File=trim(fname), Form='unformatted', Action='read&
     &', Access='direct', Recl=Recl)
      Read (un, Rec=ikkp) ikkp_, iknr_, jknr_, iq_, iqr_, n1_, n2_, &
     & n3_, n4_, zm
      Close (un)
  ! check kkp-index
      If (ikkp .Ne. ikkp_) Then
         Write (*,*)
         Write (*, '("Error(getbsemat): inconsistent (k, kp)-index")')
         Write (*, '(" requested		       : ", i8)') ikkp
         Write (*, '(" stored at requested position : ", i8)') ikkp_
         Write (*,*)
         Call terminate
      End If
  ! cut matrix
      zmat (:, :, :, :) = zm (n1_-n1+1:, :n2, n1_-n1+1:, :n2)
      Deallocate (zm)
End Subroutine getbsemat
