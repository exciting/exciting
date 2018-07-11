!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
!BOP
! !ROUTINE: getoccsv
! !INTERFACE:
!
Subroutine getoccsv (vpl, occsvp)
! !USES:
      Use modmain
      Use modinput
      Use modmpi
! !DESCRIPTION:
!   The file where the (second-variational) occupation numbers are stored is
!    {\tt OCCSV.OUT}.
!   The maximum occupancies for spin-unpolarized systems is $2$, whereas for
!   spin-polarized systems it is $1$.
!   It is a direct-access binary file, the record length of which can be
!   determined
!   with the help of the array sizes and data type information.
!   One record of this file has the following structure
!
!   \begin{tabular}{|l|l|l|}
!   \hline
!   $k_{\rm lat}$ & $N_{\rm stsv}$ & $o$ \\
!   \hline
!   \end{tabular}\newline\newline
!   The following table explains the parts of the record in more detail
!
!   \begin{tabular}{|l|l|l|l|}
!   \hline
!   name & type & shape & description\\
!   \hline \hline
!   $k_{\rm lat}$ & real(8) & 3 & k-point in lattice coordinates \\ \hline
!   $N_{\rm stsv}$ & integer & 1 & number of (second-variational) states \\
!    &  &  & (without core states) \\ \hline
!   $o$ & real(8) & $N_{\rm stsv}$ & (second-variational) occupation number
!    array \\
!   \hline
!   \end{tabular}\newline\newline
!
! !REVISION HISTORY:
!   Documentation added, Dec 2009 (S. Sagmeister)
!EOP
!BOC
      Implicit None
  ! arguments
      Real (8), Intent (In) :: vpl (3)
      Real (8), Intent (Out) :: occsvp (nstsv)
  ! local variables
      Logical :: exist
      Integer :: isym, ik, koffset, i
      Integer :: recl, nstsv_
      Real (8) :: vkl_ (3), t1
      Character (256) :: filetag
      Character (256), External :: outfilenamestring
#ifdef XS
  ! added feature to access arrays for only a subset of bands
      Real (8), Allocatable :: occsv_ (:)
#endif
  ! find the k-point number
      Call findkpt (vpl, isym, ik)
  ! find the record length
!
#ifdef XS
      Inquire (IoLength=Recl) vkl_, nstsv_
#endif
#ifndef XS
      Inquire (IoLength=Recl) vkl_, nstsv_, occsvp
#endif
      filetag = trim (filetag_occsv)
!write(*,*) "getoccsv: Reading from ", outfilenamestring(filetag, ik) 

      Do i = 1, 100
         Inquire (File=outfilenamestring(filetag, ik), Exist=Exist)
         If (exist) Then
            Open (70, File=outfilenamestring(filetag, ik), Action='READ&
           &', Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
            Exit
         Else
            Call system ('sync')
            Write (*,*) "Waiting for other process to write" // ":getoc&
           &csv:" // trim (outfilenamestring(filetag, ik))
            Call sleep (5)
         End If
      End Do
      If (splittfile) Then
         koffset = ik - firstk (procofk(ik)) + 1
      Else
         koffset = ik
      End If
#ifdef XS
      Read (70, Rec=1) vkl_, nstsv_
      Close (70)
      If (nstsv .Gt. nstsv_) Then
         Write (*,*)
         Write (*, '("Error(getoccsv): invalid nstsv for k-point ", I8)&
        &') ik
         Write (*, '(" current    : ", I8)') nstsv
         Write (*, '(" OCCSV.OUT  : ", I8)') nstsv_
         Write (*, '(" file	     : ", a)') trim(outfilenamestring(filetag,ik))
         Write (*,*)
         Stop
      End If
      Allocate (occsv_(nstsv_))
      Inquire (IoLength=Recl) vkl_, nstsv_, occsv_
      Open (70, File=outfilenamestring(filetag, ik), Action='READ', &
     & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
      Read (70, Rec=koffset) vkl_, nstsv_, occsv_
  ! retreive subset
      occsvp (:) = occsv_ (:nstsv)
      Deallocate (occsv_)
#endif
#ifndef XS
      Read (70, Rec=koffset) vkl_, nstsv_, occsvp
#endif
      Close (70)
!
      t1 = Abs (vkl(1, ik)-vkl_(1)) + Abs (vkl(2, ik)-vkl_(2)) + Abs &
     & (vkl(3, ik)-vkl_(3))
      If (t1 .Gt. input%structure%epslat) Then
         Write (*,*)
         Write (*, '("Error(getoccsv): differing vectors for k-point ",&
        & I8)') ik
         Write (*, '(" current    : ", 3G18.10)') vkl (:, ik)
         Write (*, '(" OCCSV.OUT  : ", 3G18.10)') vkl_
         Write (*, '(" file	     : ", a)') trim(outfilenamestring(filetag,ik))
         Write (*,*)
         Stop
      End If
#ifndef XS
      If (nstsv .Ne. nstsv_) Then
         Write (*,*)
         Write (*, '("Error(getoccsv): differing nstsv for k-point ", I&
        &8)') ik
         Write (*, '(" current    : ", I8)') nstsv
         Write (*, '(" OCCSV.OUT  : ", I8)') nstsv_
         Write (*, '(" file	     : ", a)') trim(outfilenamestring(filetag,ik))
         Write (*,*)
         Stop
      End If
#endif
      Return
End Subroutine getoccsv
!EOC

Module m_getoccsvr
      Implicit None
Contains
!
!
      Subroutine getoccsvr (fname, isti, istf, vpl, occsvp)
         Use modmain
         Implicit None
    ! arguments
         Character (*), Intent (In) :: fname
         Integer, Intent (In) :: isti, istf
         Real (8), Intent (In) :: vpl (3)
         Real (8), Intent (Out) :: occsvp (:)
    ! local variables
         Integer :: err
         Real (8), Allocatable :: occsvt (:)
         Character (256) :: strtmp
    ! check correct shapes
         err = 0
         If ((isti .Lt. 1) .Or. (istf .Gt. nstsv) .Or. (istf .Le. &
        & isti)) Then
            Write (*,*)
            Write (*, '("Error(getoccsvr): inconsistent limits for band&
           &s:")')
            Write (*, '(" band limits  : ", 2i6)') isti, istf
            Write (*, '(" maximum value: ", i6)') nstsv
            Write (*,*)
            err = err + 1
         End If
         If (size(occsvp, 1) .Ne. (istf-isti+1)) Then
            Write (*,*)
            Write (*, '("Error(getoccsvr): output array does not match &
           &for bands:")')
            Write (*, '(" band limits 	     : ", 2i6)') isti, istf
            Write (*, '(" requested number of bands: ", i6)') istf - &
           & isti + 1
            Write (*, '(" array size		     : ", i6)') size (occsvp, 1)
            Write (*,*)
            err = err + 1
         End If
         If (err .Ne. 0) Stop
         Allocate (occsvt(nstsv))
         filetag_occsv = trim (fname)
         strtmp = trim (filext)
         filext = ''
         Call getoccsv (vpl, occsvt)
         filetag_occsv = 'OCCSV'
         filext = trim (strtmp)
         occsvp (:) = occsvt (isti:istf)
         Deallocate (occsvt)
      End Subroutine getoccsvr
End Module m_getoccsvr
