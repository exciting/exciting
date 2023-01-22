!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: readstate
! !INTERFACE:
!
!
Subroutine readstate
! !USES:
      Use modinput
      Use modmain
#ifdef XS
      Use modxs, Only: isreadstate0
#endif
! !DESCRIPTION:
!   Reads in the charge density and other relevant variables from the file
!   {\tt STATE.OUT}. Checks for version and parameter compatibility.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Logical :: spinpol_
      Integer :: iostat
      Integer :: is, ia, ias, lmmax, lm, ir, jr
      Integer :: idm, ngm, i1, i2, i3, j1, j2, j3
      Integer :: version_ (3), nspecies_, lmmaxvr_
      Integer :: natoms_, nrmt_ (maxspecies), nrmtmax_
      Integer :: ngrid_ (3), ngrtot_, ngvec_, ndmag_
      Integer :: nspinor_, ldapu_, lmmaxlu_
      Character(40) :: githash_
      character(1024) :: message
      Real (8) :: t1
      logical, external :: versions_gt
! allocatable arrays
      Integer, Allocatable :: mapir (:)
      Real (8), Allocatable :: spr_ (:, :)
      Real (8), Allocatable :: rhomt_ (:, :, :)
      Real (8), Allocatable :: rhoir_ (:)
      Real (8), Allocatable :: vclmt_ (:, :, :)
      Real (8), Allocatable :: vclir_ (:)
      Real (8), Allocatable :: vxcmt_ (:, :, :)
      Real (8), Allocatable :: vxcir_ (:)
      Real (8), Allocatable :: veffmt_ (:, :, :)
      Real (8), Allocatable :: veffir_ (:)
      Real (8), Allocatable :: magmt_ (:, :, :, :)
      Real (8), Allocatable :: magir_ (:, :)
      Real (8), Allocatable :: bxcmt_ (:, :, :, :)
      Real (8), Allocatable :: bxcir_ (:, :)
      Complex (8), Allocatable :: veffig_ (:)
      Complex (8), Allocatable :: vmatlu_ (:, :, :, :, :)
#ifdef XS
      If (isreadstate0) Then
         Open (50, File='STATE.OUT', Action='READ', Form='UNFORMATTED', &
        & Status='OLD', IoStat=IoStat)
      Else
#endif
         Open (50, File='STATE'//trim(filext), Action='READ', Form='UNF&
        &ORMATTED', Status='OLD', IoStat=IoStat)
#ifdef XS
      End If
#endif
      If (iostat .Ne. 0) Then
         Write (*,*)
         Write (*, '("Error(readstate): error opening ", A)') 'STATE' &
        & // trim (filext)
         Write (*,*)
         Stop
      End If
      Read (50) version_
      If ((version(1) .Ne. version_(1)) .Or. (version(2) .Ne. &
     & version_(2)) .Or. (version(3) .Ne. version_(3))) Then
         call warning('Warning(readstate):')
         Write(message, '(" Different versions")')
         call warning(message)
         Write(message, '(" current   : ", I3.3, ".", I3.3, ".", I3.3)') version
         call warning(message)
         Write(message, '(" STATE.OUT : ", I3.3, ".", I3.3, ".", I3.3)') version_
        call warning(message)
      End If
! versions > 10.04.14 (April 14, 2010)
      If (versions_gt(version_,refversion_gitstate)) Then
         backspace(50)
         Read (50) version_, githash_
      Else
         githash_ = ''
      End If
      if (githash_ .ne. githash) then
        call warning('Warning(readstate):')
        write(message,'(" Different version hashes")')
        call warning(message)
        Write(message, '(" current   : ", a40)') githash
        call warning(message)
        Write(message, '(" STATE.OUT : ", a40)') githash_
        call warning(message)
      end if
      Read (50) spinpol_
      Read (50) nspecies_
      If (nspecies .Ne. nspecies_) Then
         Write (*,*)
         Write (*, '("Error(readstate): differing nspecies")')
         Write (*, '(" current	 : ", I4)') nspecies
         Write (*, '(" STATE.OUT : ", I4)') nspecies_
         Write (*,*)
         Stop
      End If
      Read (50) lmmaxvr_
#ifdef XS
      If (lmmaxvr .Ne. lmmaxvr_) Then
         call warning('Warning(readstate):')
         Write (message, '(" Differing lmmaxvr values ")')
         call warning(message)
         Write (message, '(" current   : ", I4)') lmmaxvr
         call warning(message)
         Write (message, '(" STATE.OUT : ", I4)') lmmaxvr_
         call warning(message)
      End If
#endif
      Read (50) nrmtmax_
#ifdef XS
      If (nrmtmax .Ne. nrmtmax_) Then
         call warning('Warning(readstate):')
         Write (message, '(" Differing nrmtmax values ")')
         call warning(message)
         Write (message, '(" current   : ", I4)') nrmtmax
         call warning(message)
         Write (message, '(" STATE.OUT : ", I4)') nrmtmax_
         call warning(message)
      End If
#endif
      Allocate (spr_(nrmtmax_, nspecies))
      Do is = 1, nspecies
         Read (50) natoms_
         If (natoms(is) .Ne. natoms_) Then
            Write (*,*)
            Write (*, '("Error(readstate): differing natoms for species&
           & ", I4)') is
            Write (*, '(" current   : ", I4)') natoms (is)
            Write (*, '(" STATE.OUT : ", I4)') natoms_
            Write (*,*)
            Stop
         End If
         Read (50) nrmt_ (is)
#ifdef XS
         If (nrmt(is) .Ne. nrmt_(is)) Then
            call warning('Warning(readstate):')
            Write (message, '(" Differing nrmt for species", I4)') is
            call warning(message)
            Write (message, '(" current   : ", I4)') nrmt (is)
            call warning(message)
            Write (message, '(" STATE.OUT : ", I4)') nrmt_ (is)
            call warning(message)
         End If
#endif
         Read (50) spr_ (1:nrmt_(is), is)
#ifdef XS
         If (nrmt(is) .Eq. nrmt_(is)) Then
            If (any(Abs(spr(1:nrmt(is), is)-spr_(1:nrmt(is), is)) > 1.d-10)) Then
               call warning('Warning(readstate):')
               Write (message, '(" Differing spr for species ", I4)') is
               call warning(message)
               Write (message, '(" average RMS of difference   : ", G18.10)') &
              & Sqrt (sum((spr(1:nrmt(is), is)-spr_(1:nrmt(is), &
              & is))**2)/nrmt(is))
              call warning(message)
            End If
         End If
#endif
      End Do
      Read (50) ngrid_
#ifdef XS
      If (any(ngrid .Ne. ngrid_)) Then
         call warning('Warning(readstate):')
         Write (message, '(" Differing ngrid values ")')
         Write (message, '(" current   : ", 3I4)') ngrid
         Write (message, '(" STATE.OUT : ", 3I4)') ngrid_
      End If
#endif
      Read (50) ngvec_
#ifdef XS
      If (ngvec .Ne. ngvec_) Then
         call warning('Warning(readstate):')
         Write (message, '(" Differing ngvec values ")')
         call warning(message)
         Write (message, '(" current   : ", I9)') ngvec
         call warning(message)
         Write (message, '(" STATE.OUT : ", I9)') ngvec_
         call warning(message)
      End If
#endif
      Read (50) ndmag_
#ifdef XS
      If (ndmag .Ne. ndmag_) Then
         call warning('Warning(readstate):')
         Write (message, '(" Differing ndmag values ")')
         call warning(message)
         Write (message, '(" current   : ", I4)') ndmag
         call warning(message)
         Write (message, '(" STATE.OUT : ", I4)') ndmag_
         call warning(message)
      End If
#endif
      If ((spinpol_) .And. (ndmag_ .Ne. 1) .And. (ndmag_ .Ne. 3)) Then
         Write (*,*)
         Write (*, '("Error(readstate): invalid ndmag in STATE.OUT : ",&
        & I8)') ndmag_
         Write (*,*)
         Stop
      End If
! versions > 0.9.131
      If ((version_(1) .Gt. 0) .Or. (version_(2) .Gt. 9) .Or. &
     & (version_(3) .Gt. 131)) Then
         Read (50) nspinor_
         Read (50) ldapu_
         Read (50) lmmaxlu_
      Else
         ldapu_ = 0
      End If
      ngrtot_ = ngrid_ (1) * ngrid_ (2) * ngrid_ (3)
      Allocate (mapir(ngrtot))
      Allocate (rhomt_(lmmaxvr_, nrmtmax_, natmtot))
      Allocate (rhoir_(ngrtot_))
      Allocate (vclmt_(lmmaxvr_, nrmtmax_, natmtot))
      Allocate (vclir_(ngrtot_))
      Allocate (vxcmt_(lmmaxvr_, nrmtmax_, natmtot))
      Allocate (vxcir_(ngrtot_))
      Allocate (veffmt_(lmmaxvr_, nrmtmax_, natmtot))
      Allocate (veffir_(ngrtot_))
      Allocate (veffig_(ngvec_))
! read muffin-tin density
      Read (50) rhomt_, rhoir_
! read Coulomb potential (spin independent)
      Read (50) vclmt_, vclir_
! read exchange-correlation potential
      Read (50) vxcmt_, vxcir_
! read effective potential
      Read (50) veffmt_, veffir_, veffig_
! read magnetisation and effective field
      If (spinpol_) Then
         Allocate (magmt_(lmmaxvr_, nrmtmax_, natmtot, ndmag_))
         Allocate (magir_(ngrtot_, ndmag_))
         Allocate (bxcmt_(lmmaxvr_, nrmtmax_, natmtot, ndmag_))
         Allocate (bxcir_(ngrtot_, ndmag_))
         Read (50) magmt_, magir_
         Read (50) bxcmt_, bxcir_
      End If
! read LDA+U potential matrix elements
      If ((ldapu .Ne. 0) .And. (ldapu_ .Ne. 0)) Then
         Allocate (vmatlu_(lmmaxlu_, lmmaxlu_, nspinor_, nspinor_, &
        & natmtot))
         Read (50) vmatlu_
         lmmax = Min (lmmaxlu, lmmaxlu_)
         vmatlu (:, :, :, :, :) = 0.d0
         If (nspinor .Eq. nspinor_) Then
            vmatlu (1:lmmax, 1:lmmax, :, :, :) = vmatlu_ (1:lmmax, &
           & 1:lmmax, :, :, :)
         Else If ((nspinor .Eq. 1) .And. (nspinor_ .Eq. 2)) Then
            vmatlu (1:lmmax, 1:lmmax, 1, 1, :) = 0.5d0 * &
           & (vmatlu_(1:lmmax, 1:lmmax, 1, 1, :)+vmatlu_(1:lmmax, &
           & 1:lmmax, 2, 2, :))
         Else
            vmatlu (1:lmmax, 1:lmmax, 1, 1, :) = vmatlu_ (1:lmmax, &
           & 1:lmmax, 1, 1, :)
            vmatlu (1:lmmax, 1:lmmax, 2, 2, :) = vmatlu_ (1:lmmax, &
           & 1:lmmax, 1, 1, :)
         End If
         Deallocate (vmatlu_)
      End If
      Close (50)
!---------------------------!
!     muffin-tin arrays     !
!---------------------------!
      rhomt (:, :, :) = 0.d0
      vclmt (:, :, :) = 0.d0
      vxcmt (:, :, :) = 0.d0
      veffmt (:, :, :) = 0.d0
      If (associated(input%groundstate%spin)) Then
         magmt (:, :, :, :) = 0.d0
         bxcmt (:, :, :, :) = 0.d0
      End If
      lmmax = Min (lmmaxvr, lmmaxvr_)
! interpolate the old arrays on the new radial mesh
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do lm = 1, lmmax
               Call rfinterp (nrmt_(is), spr_(:, is), lmmaxvr_, &
              & rhomt_(lm, 1, ias), nrmt(is), spr(:, is), lmmaxvr, &
              & rhomt(lm, 1, ias))
               Call rfinterp (nrmt_(is), spr_(:, is), lmmaxvr_, &
              & vclmt_(lm, 1, ias), nrmt(is), spr(:, is), lmmaxvr, &
              & vclmt(lm, 1, ias))
               Call rfinterp (nrmt_(is), spr_(:, is), lmmaxvr_, &
              & vxcmt_(lm, 1, ias), nrmt(is), spr(:, is), lmmaxvr, &
              & vxcmt(lm, 1, ias))
               Call rfinterp (nrmt_(is), spr_(:, is), lmmaxvr_, &
              & veffmt_(lm, 1, ias), nrmt(is), spr(:, is), lmmaxvr, &
              & veffmt(lm, 1, ias))
            End Do
            If ((associated(input%groundstate%spin)) .And. (spinpol_)) &
           & Then
               If (ndmag .Eq. ndmag_) Then
                  Do idm = 1, ndmag
                     Do lm = 1, lmmax
                        Call rfinterp (nrmt_(is), spr_(:, is), &
                       & lmmaxvr_, magmt_(lm, 1, ias, idm), nrmt(is), &
                       & spr(:, is), lmmaxvr, magmt(lm, 1, ias, idm))
                        Call rfinterp (nrmt_(is), spr_(:, is), &
                       & lmmaxvr_, bxcmt_(lm, 1, ias, idm), nrmt(is), &
                       & spr(:, is), lmmaxvr, bxcmt(lm, 1, ias, idm))
                     End Do
                  End Do
               Else
                  Do lm = 1, lmmax
                     Call rfinterp (nrmt_(is), spr_(:, is), lmmaxvr_, &
                    & magmt_(lm, 1, ias, ndmag_), nrmt(is), spr(:, is), &
                    & lmmaxvr, magmt(lm, 1, ias, ndmag))
                     Call rfinterp (nrmt_(is), spr_(:, is), lmmaxvr_, &
                    & bxcmt_(lm, 1, ias, ndmag_), nrmt(is), spr(:, is), &
                    & lmmaxvr, bxcmt(lm, 1, ias, ndmag))
                  End Do
               End If
            End If
         End Do
      End Do
!-----------------------------!
!     interstitial arrays     !
!-----------------------------!
      rhoir (:) = 0.d0
      vclir (:) = 0.d0
      vxcir (:) = 0.d0
      veffir (:) = 0.d0
      veffig (:) = 0.d0
      If (associated(input%groundstate%spin)) Then
         magir (:, :) = 0.d0
         bxcir (:, :) = 0.d0
      End If
! map from new grid to old
      Do i3 = 0, ngrid (3) - 1
         t1 = dble (i3*ngrid_(3)) / dble (ngrid(3))
         j3 = modulo (Nint(t1), ngrid_(3))
         Do i2 = 0, ngrid (2) - 1
            t1 = dble (i2*ngrid_(2)) / dble (ngrid(2))
            j2 = modulo (Nint(t1), ngrid_(2))
            Do i1 = 0, ngrid (1) - 1
               t1 = dble (i1*ngrid_(1)) / dble (ngrid(1))
               j1 = modulo (Nint(t1), ngrid_(1))
               ir = i3 * ngrid (2) * ngrid (1) + i2 * ngrid (1) + i1 + &
              & 1
               jr = j3 * ngrid_ (2) * ngrid_ (1) + j2 * ngrid_ (1) + j1 &
              & + 1
               mapir (ir) = jr
            End Do
         End Do
      End Do
      Do ir = 1, ngrtot
         jr = mapir (ir)
         rhoir (ir) = rhoir_ (jr)
         vclir (ir) = vclir_ (jr)
         vxcir (ir) = vxcir_ (jr)
         veffir (ir) = veffir_ (jr)
      End Do
      ngm = Min (ngvec, ngvec_)
      veffig (1:ngm) = veffig_ (1:ngm)
      If ((associated(input%groundstate%spin)) .And. (spinpol_)) Then
         Do ir = 1, ngrtot
            jr = mapir (ir)
            If (ndmag .Eq. ndmag_) Then
               magir (ir, :) = magir_ (jr, :)
               bxcir (ir, :) = bxcir_ (jr, :)
            Else
               magir (ir, ndmag) = magir_ (jr, ndmag_)
               bxcir (ir, ndmag) = bxcir_ (jr, ndmag_)
            End If
         End Do
      End If
      Deallocate (mapir, spr_, rhomt_, rhoir_, vclmt_, vclir_)
      Deallocate (vxcmt_, vxcir_, veffmt_, veffir_, veffig_)
      If (spinpol_) deallocate (magmt_, magir_, bxcmt_, bxcir_)
      Return
End Subroutine
!EOC
