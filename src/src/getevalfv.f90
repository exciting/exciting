!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
!BOP
! !ROUTINE: getevalfv
! !INTERFACE:
!
Subroutine getevalfv (vpl, evalfv)
! !USES:
      Use modmain
      Use modinput
      Use modmpi
! !DESCRIPTION:
!   The file where the (first-variational) eigenvalues are stored is
!   {\tt EVALFV.OUT}.
!   It is a direct-access binary file, the record length of which can be
!   determined
!   with the help of the array sizes and data type information.
!   One record of this file has the following structure
!
!   \begin{tabular}{|l|l|l|l|}
!   \hline
!   $k_{\rm lat}$ & $N_{\rm stfv}$ & $N_{\rm spfv}$ & $E$ \\
!   \hline
!   \end{tabular}\newline\newline
!   The following table explains the parts of the record in more detail
!
!   \begin{tabular}{|l|l|l|l|}
!   \hline
!   name & type & shape & description\\
!   \hline \hline
!   $k_{\rm lat}$ & real(8) & 3 & k-point in lattice coordinates \\ \hline
!   $N_{\rm stfv}$ & integer & 1 & number of (first-variational) states \\
!    &  &  & (without core states) \\ \hline
!   $N_{\rm spfv}$ & integer & 1 & first-variational spins (always equals 1)
!       \\ \hline
!   $E$ & real(8) & $N_{\rm stfv}\times N_{\rm spfv}$ & (first-variational)
!    eigenvalue array \\
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
      Real (8), Intent (Out) :: evalfv (nstfv, nspnfv)
! local variables
      Integer :: isym, ik, koffset, i
      Logical :: exist
      Integer :: recl, nstfv_, nspnfv_
      Real (8) :: vkl_ (3), t1
!
      Character (256) :: filetag
      Character (256), External :: outfilenamestring
#ifdef XS
  ! added feature to access arrays for only a subset of bands
      Real (8), Allocatable :: evalfv_ (:, :)
#endif
! find the k-point number
      Call findkpt (vpl, isym, ik)
! find the record length
!
#ifdef XS
      Inquire (IoLength=Recl) vkl_, nstfv_, nspnfv_
#endif
#ifndef XS
      Inquire (IoLength=Recl) vkl_, nstfv_, nspnfv_, evalfv
#endif
!
      filetag = trim (filetag_evalfv)
      Do i = 1, 100
         Inquire (File=outfilenamestring(filetag, ik), Exist=Exist)
         If (exist) Then
            Open (70, File=outfilenamestring(filetag, ik), Action='READ&
           &', Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
            Exit
         Else
            Call system ('sync')
            Write (*,*) "Waiting for other process to write" // ":getev&
           &alfv:" // trim (outfilenamestring(filetag, ik))
            Call sleep (5)
         End If
      End Do
!
      If (splittfile) Then
         koffset = ik - firstk (procofk(ik)) + 1
      Else
         koffset = ik
      End If
!
#ifdef XS
      Read (70, Rec=1) vkl_, nstfv_, nspnfv_
      Close (70)
      If (nstfv .Gt. nstfv_) Then
         Write (*,*)
         Write (*, '("Error(getevalfv): invalid nstfv for k-point ", I8&
        &)') ik
         Write (*, '(" current    : ", I8)') nstfv
         Write (*, '(" EVALFV.OUT : ", I8)') nstfv_
         Write (*, '(" file	     : ", a	 )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
      Allocate (evalfv_(nstfv_, nspnfv_))
      Inquire (IoLength=Recl) vkl_, nstfv_, nspnfv_, evalfv_
      Open (70, File=outfilenamestring(filetag, ik), Action='READ', &
     & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
      Read (70, Rec=koffset) vkl_, nstfv_, nspnfv_, evalfv_
  ! retreive subset
      evalfv (:, :) = evalfv_ (:nstfv, :)
      Deallocate (evalfv_)
#endif
!
#ifndef XS
      Read (70, Rec=koffset) vkl_, nstfv_, nspnfv_, evalfv
#endif
      Close (70)
!
!
!
      t1 = Abs (vkl(1, ik)-vkl_(1)) + Abs (vkl(2, ik)-vkl_(2)) + Abs &
     & (vkl(3, ik)-vkl_(3))
      If (t1 .Gt. input%structure%epslat) Then
         Write (*,*)
         Write (*, '("Error(getevalfv): differing vectors for k-point "&
        &, I8)') ik
         Write (*, '(" current	  : ", 3G18.10)') vkl (:, ik)
         Write (*, '(" EVALFV.OUT : ", 3G18.10)') vkl_
         Write (*, '(" file	  : ", a      )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
#ifndef XS
      If (nstfv .Ne. nstfv_) Then
         Write (*,*)
         Write (*, '("Error(getevalfv): differing nstfv for k-point ", &
        &I8)') ik
         Write (*, '(" current	  : ", I8)') nstfv
         Write (*, '(" EVALFV.OUT : ", I8)') nstfv_
         Write (*, '(" file	  : ", a      )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
#endif
      If (nspnfv .Ne. nspnfv_) Then
         Write (*,*)
         Write (*, '("Error(getevalfv): differing nspnfv for k-point ",&
        & I8)') ik
         Write (*, '(" current	  : ", I8)') nspnfv
         Write (*, '(" EVALFV.OUT : ", I8)') nspnfv_
         Write (*, '(" file	  : ", a      )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
      Return
End Subroutine
!EOC

Module m_getevalfvr
      Implicit None
Contains
!
!
      Subroutine getevalfvr (fname, isti, istf, vpl, evalfv)
         Use modmain
         Implicit None
    ! arguments
         Character (*), Intent (In) :: fname
         Integer, Intent (In) :: isti, istf
         Real (8), Intent (In) :: vpl (3)
         Real (8), Intent (Out) :: evalfv (:, :)
    ! local variables
         Integer :: err
         Real (8), Allocatable :: evalfvt (:, :)
         Character (256) :: str1
    ! check correct shapes
         err = 0
         If ((isti .Lt. 1) .Or. (istf .Gt. nstfv) .Or. (istf .Le. &
        & isti)) Then
            Write (*,*)
            Write (*, '("Error(getevalfvr): inconsistent limits for ban&
           &ds:")')
            Write (*, '(" band limits  : ", 2i6)') isti, istf
            Write (*, '(" maximum value: ", i6)') nstfv
            Write (*,*)
            err = err + 1
         End If
         If (size(evalfv, 1) .Ne. (istf-isti+1)) Then
            Write (*,*)
            Write (*, '("Error(getevalfvr): output array does not match&
           & for bands:")')
            Write (*, '(" band limits 	     : ", 2i6)') isti, istf
            Write (*, '(" requested number of bands: ", i6)') istf - &
           & isti + 1
            Write (*, '(" array size		     : ", i6)') size (evalfv, 1)
            Write (*,*)
            err = err + 1
         End If
         If (size(evalfv, 2) .Ne. nspnfv) Then
            Write (*,*)
            Write (*, '("Error(getevalfvr): output array does not match&
           & for nspnfv:")')
            Write (*, '(" nspnfv    : ", i6)') nspnfv
            Write (*, '(" array size: ", i6)') size (evalfv, 2)
            Write (*,*)
            err = err + 1
         End If
         If (err .Ne. 0) Stop
         Allocate (evalfvt(nstfv, nspnfv))
         filetag_evalfv = trim (fname)
         str1 = trim (filext)
         filext = ''
         Call getevalfv (vpl, evalfvt)
         filetag_evalfv = 'EVALFV'
         filext = trim (str1)
         evalfv (:, :) = evalfvt (isti:istf, :)
         Deallocate (evalfvt)
      End Subroutine getevalfvr
End Module m_getevalfvr
