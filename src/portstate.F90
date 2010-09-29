
! Copyright (C) 2002-2007 S. Sagmeister J. K. Dewhurst, S. Sharma and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: portstate
! !INTERFACE:
Subroutine portstate (act)
! !USES:
      Use ioarray
      use mod_misc, only: refversion_gitstate
! !DESCRIPTION:
!   Toggle file format of {\tt STATE.OUT}. If tb2a is true an ASCII
!   file with the name {\tt STATE.xml} is generated and the data
!   from {\tt STATE.OUT} is transferred. If tb2a is false the conversion
!   goes in the other direction. Based upon the routines {\tt readstate}
!   and {\tt writestate}.
!
! !REVISION HISTORY:
!   Created 2007 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Integer, Intent (In) :: act
  ! local variables
      Logical :: tb2a, exist, spinpol_
      Integer :: natmtot, is
      Integer :: version_ (3), nspecies_, lmmaxvr_, nrmtmax_
      Integer :: natoms_ (10000), ngrid_ (3)
      Integer :: ngrtot_, ngvec_, ndmag_, nspinor_, ldapu_, lmmaxlu_
      logical, external :: versions_gt
      Character(40) :: githash_
  ! allocatable arrays
      Integer, Allocatable :: nrmt_ (:)
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
      Select Case (act)
      Case (1, 2,-1,-2)
      Case Default
         Write (*,*)
         Write (*, '("Error(portstate): unknown action: ", i6)') act
         Write (*,*)
         Stop
      End Select
      tb2a = (act .Eq. 1) .Or. (act .Eq.-1)
      If (tb2a) Then
         Open (50, File='STATE.OUT', Action='READ', Form='UNFORMATTED', &
        & Status='OLD')
         If (act .Eq.-1) Then
            Write (*,*)
            Write (*, '("Information on STATE.OUT file:")')
            Write (*,*)
         End If
         Inquire (File='STATE.xml', Exist=Exist)
         If (exist .And. (act .Eq. 1)) Then
            Write (*,*)
            Write (*, '("Error(portstate): not overwriting existent STA&
           &TE.xml file")')
            Write (*,*)
            return
         End If
         If (act .Eq. 1) open (51, file='STATE.xml', action='WRITE', &
        & form='FORMATTED', status='replace')
         Read (50) version_
         ! versions >= 10.4 (after April 2010)
         if (versions_gt(version_,refversion_gitstate)) then
            backspace(50)
            read(50) version_, githash_
         else
            githash_=''
         end if
         Read (50) spinpol_
         Read (50) nspecies_
         Read (50) lmmaxvr_
         Read (50) nrmtmax_
         If (act .Eq.-1) Then
            Write (*, '("version     :", 3i8)') version_
            if (versions_gt(version_,refversion_gitstate)) &
              Write (*, '("versionhash : ", a40)') githash_
            Write (*, '("spinpol     :", l8)') spinpol_
            Write (*, '("nspecies    :", i8)') nspecies_
            Write (*, '("lmaxvr      :", i8)') nint(sqrt(dble(lmmaxvr_)))-1
            Write (*, '("nrmtmax     :", 3i8)') nrmtmax_
         End If
         If (act .Eq. 1) then
            Write (51, '(a)') '<?xml version="1.0"?>'
            Write (51, '(a)') '<state>'
            Write (51, '(a)') '<data name = "version" type = "integer" dim&
           &ension = "1" shape = "3">'
            Call ioarr (un=51, ioa='write', arr1di=version_)
            Write (51, '(a)') '</data>'
            if (versions_gt(version_,refversion_gitstate)) then
              Write (51, '(a)') '<data name = "versionhash" type = "character(40)" &
                &dimension = "1" shape = "1">'
              Write (51, *) githash_
              Write (51, '(a)') '</data>'
            end if
            Write (51, '(a)') '<data name = "spinpol" type = "logical" dim&
           &vension = "1" shape = "1">'
            Write (51,*) spinpol_
            Write (51, '(a)') '</data>'
            Write (51, '(a)') '<data name = "nspecies" type = "integer" di&
           &mension = "1" shape = "1">'
            Write (51,*) nspecies_
            Write (51, '(a)') '</data>'
            Write (51, '(a)') '<data name = "lmmaxvr" type = "integer" dim&
           &ension = "1" shape = "1">'
            Write (51,*) lmmaxvr_
            Write (51, '(a)') '</data>'
            Write (51, '(a)') '<data name = "nrmtmax" type = "integer" dim&
           &ension = "1" shape = "1">'
            Write (51,*) nrmtmax_
            Write (51, '(a)') '</data>'
         end if
      Else
         Open (50, File='STATE.xml', Action='READ', Form='FORMATTED', &
        & Status='OLD')
         If (act .Eq.-2) Then
            Write (*,*)
            Write (*, '("Information on STATE.xml file:")')
            Write (*,*)
         End If
         Inquire (File='STATE.OUT', Exist=Exist)
         If (exist .And. (act .Eq. 2)) Then
            Write (*,*)
            Write (*, '("Error(portstate): not overwriting existent STA&
           &TE.OUT file")')
            Write (*,*)
            return
         End If
         If (act .Eq. 2) open (51, file='STATE.OUT', action='WRITE', &
        & form='UNFORMATTED', status='replace')
         Read (50,*)
         Read (50,*)
         Read (50,*)
         Call ioarr (un=50, ioa='read', arr1di=version_)
         if (versions_gt(version_,refversion_gitstate)) then
           read (50, *)
           read (50, *)
           read (50, *) githash_
         end if
         Read (50,*)
         Read (50,*)
         Read (50,*) spinpol_
         Read (50,*)
         Read (50,*)
         Read (50,*) nspecies_
         Read (50,*)
         Read (50,*)
         Read (50,*) lmmaxvr_
         Read (50,*)
         Read (50,*)
         Read (50,*) nrmtmax_
         Read (50,*)
         If (act .Eq.-2) Then
            Write (*, '("version     :", 3i8)') version_
            if (versions_gt(version_,refversion_gitstate)) &
            Write (*, '("versionhash : ", a40)') githash_
            Write (*, '("spinpol     :", l8)') spinpol_
            Write (*, '("nspecies    :", i8)') nspecies_
            Write (*, '("lmaxvr      :", i8)') nint(sqrt(dble(lmmaxvr_)))-1
            Write (*, '("nrmtmax     :", i8)') nrmtmax_
         End If
         if (act .eq. 2) then
            if (versions_gt(version_,refversion_gitstate)) then
              Write (51) version_, githash_
            else
              Write (51) version_
            end if
            Write (51) spinpol_
            Write (51) nspecies_
            Write (51) lmmaxvr_
            Write (51) nrmtmax_
         end if
      End If
      Allocate (spr_(nrmtmax_, nspecies_))
      Allocate (nrmt_(nspecies_))
      If (tb2a) Then
         natmtot = 0
         Do is = 1, nspecies_
            Read (50) natoms_ (is)
            Read (50) nrmt_ (is)
            Read (50) spr_ (1:nrmt_(is), is)
            if (act .eq. 1) then
               Write (51, '(a)') '<data name = "natoms" type = "integer" d&
              &imension = "1" shape = "1" index = "species" indexval = "' &
              & // trim (i2str(is)) // '">'
               Write (51,*) natoms_ (is)
               Write (51, '(a)') '</data>'
               Write (51, '(a)') '<data name = "nrmt" type = "integer" dim&
              &ension = "1" shape = "1" index = "species" indexval = "' // &
              & trim (i2str(is)) // '">'
               Write (51,*) nrmt_ (is)
               Write (51, '(a)') '</data>'
               Write (51, '(a)') '<data name = "spr" type = "real(8)" dime&
              &nsion = "1" shape = "' // trim (i2str(nrmt_(is))) // '" ind&
              &ex = "species" indexval = "' // trim (i2str(is)) // '">'
               Call ioarr (un=51, ioa='write', arr1dr=spr_(1:nrmt_(is), &
              & is))
               Write (51, '(a)') '</data>'
            end if
            natmtot = natmtot + natoms_ (is)
         End Do
         Read (50) ngrid_
         Read (50) ngvec_
         Read (50) ndmag_
     ! versions > 0.9.131
         If ((version_(1) .Gt. 0) .Or. (version_(2) .Gt. 9) .Or. &
        & (version_(3) .Gt. 131)) Then
            Read (50) nspinor_
            Read (50) ldapu_
            Read (50) lmmaxlu_
         Else
            nspinor_ = 1
            If (spinpol_) nspinor_ = 2
            ldapu_ = 0
            lmmaxlu_ = 0
         End If
         If (act .Eq.-1) Then
            Write (*, '("ngrid       :", 3i8)') ngrid_
            Write (*, '("ngvec       :", i8)') ngvec_
            Write (*, '("ndmag       :", i8)') ndmag_
            Write (*, '("nspinor     :", i8)') nspinor_
            Write (*, '("ldapu       :", i8)') ldapu_
            Write (*, '("lmaxlu      :", i8)') nint(sqrt(dble(lmmaxlu_)))-1
            Write (*,*)
            ! return here if in extracting mode (act=-1,-2)
            return
         End If
         Write (51, '(a)') '<data name = "ngrid" type = "integer" dimen&
        &sion = "1" shape = "3">'
         Call ioarr (un=51, ioa='write', arr1di=ngrid_)
         Write (51, '(a)') '</data>'
         Write (51, '(a)') '<data name = "ngvec" type = "integer" dimen&
        &sion = "1" shape = "1">'
         Write (51,*) ngvec_
         Write (51, '(a)') '</data>'
         Write (51, '(a)') '<data name = "ndmag" type = "integer" dimen&
        &sion = "1" shape = "1">'
         Write (51,*) ndmag_
         Write (51, '(a)') '</data>'
         Write (51, '(a)') '<data name = "nspinor" type = "integer" dim&
        &ension = "1" shape = "1">'
         Write (51,*) nspinor_
         Write (51, '(a)') '</data>'
         Write (51, '(a)') '<data name = "ldapu" type = "integer" dimen&
        &sion = "1" shape = "1">'
         Write (51,*) ldapu_
         Write (51, '(a)') '</data>'
         Write (51, '(a)') '<data name = "lmmaxlu" type = "integer" dim&
        &ension = "1" shape = "1">'
         Write (51,*) lmmaxlu_
         Write (51, '(a)') '</data>'
      Else
         natmtot = 0
         Do is = 1, nspecies_
            Read (50,*)
            Read (50,*) natoms_ (is)
            Read (50,*)
            Read (50,*)
            Read (50,*) nrmt_ (is)
            Read (50,*)
            Read (50,*)
            Call ioarr (un=50, ioa='read', arr1dr=spr_(1:nrmt_(is), &
           & is))
            Read (50,*)
            if (act .eq. 2) then
               Write (51) natoms_ (is)
               Write (51) nrmt_ (is)
               Write (51) spr_ (1:nrmt_(is), is)
            end if
            natmtot = natmtot + natoms_ (is)
         End Do
         Read (50,*)
         Call ioarr (un=50, ioa='read', arr1di=ngrid_)
         Read (50,*)
         Read (50,*)
         Read (50,*) ngvec_
         Read (50,*)
         Read (50,*)
         Read (50,*) ndmag_
         Read (50,*)
         Read (50,*)
         Read (50,*) nspinor_
         Read (50,*)
         Read (50,*)
         Read (50,*) ldapu_
         Read (50,*)
         Read (50,*)
         Read (50,*) lmmaxlu_
         Read (50,*)
         If (act .Eq.-2) Then
            Write (*, '("ngrid       :", 3i8)') ngrid_
            Write (*, '("ngvec       :", i8)') ngvec_
            Write (*, '("ndmag       :", i8)') ndmag_
            Write (*, '("nspinor     :", i8)') nspinor_
            Write (*, '("ldapu       :", i8)') ldapu_
            Write (*, '("lmaxlu      :", i8)') nint(sqrt(dble(lmmaxlu_)))-1
            Write (*,*)
            ! return here if in extracting mode (act=-1,-2)
            return
         End If
         Write (51) ngrid_
         Write (51) ngvec_
         Write (51) ndmag_
         Write (51) nspinor_
         Write (51) ldapu_
         Write (51) lmmaxlu_
      End If
      ngrtot_ = ngrid_ (1) * ngrid_ (2) * ngrid_ (3)
      Allocate (rhomt_(lmmaxvr_, nrmtmax_, natmtot))
      Allocate (rhoir_(ngrtot_))
      Allocate (vclmt_(lmmaxvr_, nrmtmax_, natmtot))
      Allocate (vclir_(ngrtot_))
      Allocate (vxcmt_(lmmaxvr_, nrmtmax_, natmtot))
      Allocate (vxcir_(ngrtot_))
      Allocate (veffmt_(lmmaxvr_, nrmtmax_, natmtot))
      Allocate (veffir_(ngrtot_))
      Allocate (veffig_(ngvec_))
      If (spinpol_) Then
         Allocate (magmt_(lmmaxvr_, nrmtmax_, natmtot, ndmag_))
         Allocate (magir_(ngrtot_, ndmag_))
         Allocate (bxcmt_(lmmaxvr_, nrmtmax_, natmtot, ndmag_))
         Allocate (bxcir_(ngrtot_, ndmag_))
      End If
      If (ldapu_ .Ne. 0) Then
         Allocate (vmatlu_(lmmaxlu_, lmmaxlu_, nspinor_, nspinor_, &
        & natmtot))
      End If
      If (tb2a) Then
     ! read muffin-tin density
         Read (50) rhomt_, rhoir_
     ! read Coulomb potential (spin independent)
         Read (50) vclmt_, vclir_
     ! read exchange-correlation potential
         Read (50) vxcmt_, vxcir_
     ! read effective potential
         Read (50) veffmt_, veffir_, veffig_
     ! write the density
         Write (51, '(a)') '<data name = "rhomt" type = "real(8)" dimen&
        &sion = "3" shape = "' // trim (i2str(lmmaxvr_)) // ', ' // &
        & trim (i2str(nrmtmax_)) // ', ' // trim (i2str(natmtot)) // '"&
        &>'
         Call ioarr (un=51, ioa='write', arr3dr=rhomt_)
         Write (51, '(a)') '</data>'
         Write (51, '(a)') '<data name = "rhoir" type = "real(8)" dimen&
        &sion = "1" shape = "' // trim (i2str(ngrtot_)) // '">'
         Call ioarr (un=51, ioa='write', arr1dr=rhoir_)
         Write (51, '(a)') '</data>'
     ! write the Coulomb potential
         Write (51, '(a)') '<data name = "vclmt" type = "real(8)" dimen&
        &sion = "3" shape = "' // trim (i2str(lmmaxvr_)) // ', ' // &
        & trim (i2str(nrmtmax_)) // ', ' // trim (i2str(natmtot)) // '"&
        &>'
         Call ioarr (un=51, ioa='write', arr3dr=vclmt_)
         Write (51, '(a)') '</data>'
         Write (51, '(a)') '<data name = "vclir" type = "real(8)" dimen&
        &sion = "1" shape = "' // trim (i2str(ngrtot_)) // '">'
         Call ioarr (un=51, ioa='write', arr1dr=vclir_)
         Write (51, '(a)') '</data>'
     ! write the exchange-correlation potential
         Write (51, '(a)') '<data name = "vxcmt" type = "real(8)" dimen&
        &sion = "3" shape = "' // trim (i2str(lmmaxvr_)) // ', ' // &
        & trim (i2str(nrmtmax_)) // ', ' // trim (i2str(natmtot)) // '"&
        &>'
         Call ioarr (un=51, ioa='write', arr3dr=vxcmt_)
         Write (51, '(a)') '</data>'
         Write (51, '(a)') '<data name = "vxcir" type = "real(8)" dimen&
        &sion = "1" shape = "' // trim (i2str(ngrtot_)) // '">'
         Call ioarr (un=51, ioa='write', arr1dr=vxcir_)
         Write (51, '(a)') '</data>'
     ! write the effective potential
         Write (51, '(a)') '<data name = "veffmt" type = "real(8)" dime&
        &nsion = "3" shape = "' // trim (i2str(lmmaxvr_)) // ', ' // &
        & trim (i2str(nrmtmax_)) // ', ' // trim (i2str(natmtot)) // '"&
        &>'
         Call ioarr (un=51, ioa='write', arr3dr=veffmt_)
         Write (51, '(a)') '</data>'
         Write (51, '(a)') '<data name = "veffir" type = "real(8)" dime&
        &nsion = "1" shape = "' // trim (i2str(ngrtot_)) // '">'
         Call ioarr (un=51, ioa='write', arr1dr=veffir_)
         Write (51, '(a)') '</data>'
         Write (51, '(a)') '<data name = "veffig" type = "complex(8)" d&
        &imension = "1" shape = "' // trim (i2str(ngvec_)) // '">'
         Call ioarr (un=51, ioa='write', arr1dc=veffig_)
         Write (51, '(a)') '</data>'
         If (spinpol_) Then
        ! read magnetisation and effective field
            Read (50) magmt_, magir_
            Read (50) bxcmt_, bxcir_
        ! write the magnetisation and effective magnetic fields
            Write (51, '(a)') '<data name = "magmt" type = "real(8)" di&
           &mension = "4" shape = "' // trim (i2str(lmmaxvr_)) // ', ' &
           & // trim (i2str(nrmtmax_)) // ', ' // trim (i2str(natmtot)) &
           & // ', ' // trim (i2str(ndmag_)) // '">'
            Call ioarr (un=51, ioa='write', arr4dr=magmt_)
            Write (51, '(a)') '</data>'
            Write (51, '(a)') '<data name = "magir" type = "real(8)" di&
           &mension = "2" shape = "' // trim (i2str(ngrtot_)) // ', ' &
           & // trim (i2str(ndmag_)) // '">'
            Call ioarr (un=51, ioa='write', arr2dr=magir_)
            Write (51, '(a)') '</data>'
            Write (51, '(a)') '<data name = "bxcmt" type = "real(8)" di&
           &mension = "4" shape = "' // trim (i2str(lmmaxvr_)) // ', ' &
           & // trim (i2str(nrmtmax_)) // ', ' // trim (i2str(natmtot)) &
           & // ', ' // trim (i2str(ndmag_)) // '">'
            Call ioarr (un=51, ioa='write', arr4dr=bxcmt_)
            Write (51, '(a)') '</data>'
            Write (51, '(a)') '<data name = "bxcir" type = "real(8)" di&
           &mension = "2" shape = "' // trim (i2str(ngrtot_)) // ', ' &
           & // trim (i2str(ndmag_)) // '">'
            Call ioarr (un=51, ioa='write', arr2dr=bxcir_)
            Write (51, '(a)') '</data>'
         End If
         If (ldapu_ .Ne. 0) Then
        ! read the LDA+U potential matrix elements
            Read (50) vmatlu_
        ! write the LDA+U potential matrix elements
            Write (51, '(a)') '<data name = "vmatlu" type = "complex(8)&
           &" dimension = "5" shape = "' // trim (i2str(lmmaxlu_)) // '&
           &, ' // trim (i2str(lmmaxlu_)) // ', ' // trim &
           & (i2str(nspinor_)) // ', ' // trim (i2str(nspinor_)) // ', &
           &' // trim (i2str(natmtot)) // '">'
            Call ioarr (un=51, ioa='write', arr5dc=vmatlu_)
            Write (51, '(a)') '</data>'
         End If
         Write (51, '(a)') '</state>'
      Else
     ! read muffin-tin density
         Read (50,*)
         Call ioarr (un=50, ioa='read', arr3dr=rhomt_)
         Read (50,*)
         Read (50,*)
         Call ioarr (un=50, ioa='read', arr1dr=rhoir_)
         Read (50,*)
     ! read Coulomb potential (spin independent)
         Read (50,*)
         Call ioarr (un=50, ioa='read', arr3dr=vclmt_)
         Read (50,*)
         Read (50,*)
         Call ioarr (un=50, ioa='read', arr1dr=vclir_)
         Read (50,*)
     ! read exchange-correlation potential
         Read (50,*)
         Call ioarr (un=50, ioa='read', arr3dr=vxcmt_)
         Read (50,*)
         Read (50,*)
         Call ioarr (un=50, ioa='read', arr1dr=vxcir_)
         Read (50,*)
     ! read effective potential
         Read (50,*)
         Call ioarr (un=50, ioa='read', arr3dr=veffmt_)
         Read (50,*)
         Read (50,*)
         Call ioarr (un=50, ioa='read', arr1dr=veffir_)
         Read (50,*)
         Read (50,*)
         Call ioarr (un=50, ioa='read', arr1dc=veffig_)
         Read (50,*)
     ! write the density
         Write (51) rhomt_, rhoir_
     ! write the Coulomb potential
         Write (51) vclmt_, vclir_
     ! write the exchange-correlation potential
         Write (51) vxcmt_, vxcir_
     ! write the effective potential
         Write (51) veffmt_, veffir_, veffig_
         If (spinpol_) Then
        ! read magnetisation and effective field
            Read (50,*)
            Call ioarr (un=50, ioa='read', arr4dr=magmt_)
            Read (50,*)
            Read (50,*)
            Call ioarr (un=50, ioa='read', arr2dr=magir_)
            Read (50,*)
            Read (50,*)
            Call ioarr (un=50, ioa='read', arr4dr=bxcmt_)
            Read (50,*)
            Read (50,*)
            Call ioarr (un=50, ioa='read', arr2dr=bxcir_)
            Read (50,*)
        ! write the magnetisation and effective magnetic fields
            Write (51) magmt_, magir_
            Write (51) bxcmt_, bxcir_
         End If
         If (ldapu_ .Ne. 0) Then
        ! read the LDA+U potential matrix elements
            Read (50,*)
            Call ioarr (un=50, ioa='read', arr5dc=vmatlu_)
            Read (50,*)
        ! write the LDA+U potential matrix elements
            Write (51) vmatlu_
         End If
      End If
      Close (50)
      Close (51)
      Deallocate (nrmt_, spr_, rhomt_, rhoir_, vclmt_, vclir_)
      Deallocate (vxcmt_, vxcir_, veffmt_, veffir_, veffig_)
      If (spinpol_) deallocate (magmt_, magir_, bxcmt_, bxcir_)
      If (tb2a) Then
         Write (*,*)
         Write (*, '("Info(portstate): generated STATE.xml file from STATE.OUT file")')
         Write (*,*)
      Else
         Write (*,*)
         Write (*, '("Info(portstate): generated STATE.OUT file from STATE.xml file.")')
         Write (*, '(" There is currently no strict XML parsing implemented. Therefore,")')
         Write (*, '(" if you convert a modified XML file or an XML file, not created")')
         Write (*, '(" by this utility, beware that currently there is a particular.")')
         Write (*, '(" sequence of XML-elements assumed, as well as no blank lines")')
         Write (*, '(" and comments are allowed.")')
         Write (*, '(" The numbers, however, are parsed in a free (Fortran) format.")')
         Write (*, '(" Use the resulting STATE.OUT file at your own risk.")')
         Write (*,*)
      End If
Contains
      Character (256) Function i2str (i)
    ! arguments
         Integer, Intent (In) :: i
    ! local variables
         Character (1024) :: str
         Write (str,*) i
         i2str = trim (adjustl(str))
      End Function i2str
End Subroutine portstate
!EOC
