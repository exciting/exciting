!
!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_fxc_bse_ma03
      Implicit None
!
Contains
!
!BOP
! !ROUTINE: fxc_bse_ma03
! !INTERFACE:
!
!
      Subroutine fxc_bse_ma03 (msiz, oct, sw, iw, fxc)
! !USES:
         Use modinput
         Use mod_constants, Only: zzero
         Use modmpi, Only:
         Use modxs, Only: unitout, bzsampl
         Use invert
         Use m_xsgauntgen
         Use m_findgntn0
         Use m_writegqpts
         Use m_genfilname
         Use m_getunit
! !INPUT/OUTPUT PARAMETERS:
!   msiz  : matrix size of local field effects (in,integer)
!   sw    : true for inclusion of local field effects (in,logical)
!   alpha : real constant (in,real)
!   fxc   : xc-kernel Fourier coefficients (out,complex(:,:))
! !DESCRIPTION:
!   BSE-kernel of A. Marini, Phys. Rev. Lett. 91, 256402 (2003).
!   Interface function.
!
! !REVISION HISTORY:
!   Created March 2008 (Sagmeister)
!EOP
!BOC
         Implicit None
    ! arguments
         Integer, Intent (In) :: msiz, oct
    ! true if all G-components of fxc are to be considered
         Logical, Intent (In) :: sw
         Integer, Intent (In) :: iw
         Complex (8), Intent (Out) :: fxc (:, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'fxc_bse_ma03'
         Character (256) :: filnam
         Complex (8), Allocatable :: fxch (:, :), fxcw1 (:, :), fxcw2 &
        & (:, :)
         Complex (8) :: zt1
         Integer :: n, n2, un, recl, n_
         n = size (fxc, 1)
         n2 = size (fxc, 2)
         If ((n .Lt. msiz) .Or. (n .Ne. n2)) Then
            Write (unitout, '(a, 2i9, a, i9, a)') 'Error(' // trim &
           & (thisnam) // '): size of fxc is inconsistent (required)', &
           & n, n2, '(', msiz, ')'
            Call terminate
         End If
         Allocate (fxch(-3:-1,-3:-1), fxcw1(-3:-1, n), fxcw2(n,-3:-1))
    ! filename for BSE-xc-kernel
         Call getunit (un)
    ! filename for xc-kernel
         Call genfilname (basename='FXC_BSE', asc=.False., &
        & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .Not. &
        & input%xs%tddft%aresfxc, tord=input%xs%tddft%tordfxc, iqmt=1, &
        & filnam=filnam)
	  ! get LFE size
         Inquire (IoLength=Recl) n_
         Open (un, File=trim(filnam), Form='unformatted', Action='read',&
        &  Status='old', Access='direct', Recl=Recl)
         Read (un, Rec=1) n_
         Close (un)
    ! check if data from file can be stored in local array
         If (n .Lt. n_) Then
            Write (unitout,*)
            Write (unitout, '("Error(", a, "): LFE size of file too lar&
           &ge")') trim (thisnam)
            Write (unitout, '(" LFE size (file)    : ", i8)') n_
            Write (unitout, '(" LFE size (current) : ", i8)') n
            Write (unitout,*)
            Call terminate
         End If
	  ! get data
         fxc (:, :) = zzero
         Inquire (IoLength=Recl) n_, fxch, fxcw1 (:, :n_), fxcw2 (:n_, &
        & :), fxc (:n_, :n_)
         Open (un, File=trim(filnam), Form='unformatted', Action='read',&
        &  Status='old', Access='direct', Recl=Recl)
         Read (un, Rec=iw) n_, fxch, fxcw1 (:, :n_), fxcw2 (:n_, :), &
        & fxc (:n_, :n_)
    ! assign head
         fxc (1, 1) = fxch (-oct,-oct)
    ! assign wings
         If (msiz .Gt. 1) Then
            fxc (1, 2:n_) = fxcw1 (-oct, 2:n_)
            fxc (2:n_, 1) = fxcw2 (2:n_,-oct)
         End If
    ! no LFE at all
         If ( .Not. sw) Then
            zt1 = fxc (1, 1)
            fxc (:, :) = zzero
            fxc (1, 1) = zt1
         End If
         Close (un)
         Deallocate (fxch, fxcw1, fxcw2)
      End Subroutine fxc_bse_ma03
!EOC
!
End Module m_fxc_bse_ma03
