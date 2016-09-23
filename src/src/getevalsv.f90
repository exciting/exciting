!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
!BOP
! !ROUTINE: getevalsv
! !INTERFACE:
!
Subroutine getevalsv (vpl, evalsvp)
! !USES:
      Use modmain
      Use modinput
      Use modmpi
! !DESCRIPTION:
!   The file where the (second-variational) eigenvalues are stored is
!    {\tt EVALSV.OUT}.
!   It is a direct-access binary file, the record length of which can be
!   determined
!   with the help of the array sizes and data type information.
!   One record of this file has the following structure
!
!   \begin{tabular}{|l|l|l|}
!   \hline
!   $k_{\rm lat}$ & $N_{\rm stsv}$ & $E$ \\
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
!   $E$ & real(8) & $N_{\rm stsv}$ & (second-variational) eigenvalue array \\
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
      Real (8), Intent (Out) :: evalsvp (nstsv)
  ! local variables
      Logical :: exist
      Integer :: isym, ik, koffset, i
      Integer :: recl, nstsv_
      Real (8) :: vkl_ (3), t1
      Character (256) :: filetag
      Character (256), External :: outfilenamestring
!
#ifdef XS
  ! added feature to access arrays for only a subset of bands
      Real (8), Allocatable :: evalsv_ (:)
#endif
  ! find the k-point number
      Call findkpt (vpl, isym, ik)
!
  ! find the record length
#ifdef XS
      Inquire (IoLength=Recl) vkl_, nstsv_
#endif
#ifndef XS
      Inquire (IoLength=Recl) vkl_, nstsv_, evalsvp
#endif
      filetag = trim (filetag_evalsv)
      Do i = 1, 100
         Inquire (File=outfilenamestring(filetag, ik), Exist=Exist)
         If (exist) Then
            Open (70, File=outfilenamestring(filetag, ik), Action='READ&
           &', Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
            Exit
         Else
            Call system ('sync')
            Write (*,*) "Waiting for other process to write" // ":getev&
           &alsv:" // trim (outfilenamestring(filetag, ik))
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
         Write (*, '("Error(getevalsv): invalid nstsv for k-point ",I8)&
        &') ik
         Write (*, '(" current    : ",I8)') nstsv
         Write (*, '(" EVALSV.OUT : ",I8)') nstsv_
         Write (*, '(" file       : ",a      )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
      Allocate (evalsv_(nstsv_))
      Inquire (IoLength=Recl) vkl_, nstsv_, evalsv_
      Open (70, File=outfilenamestring(filetag, ik), Action='READ', &
     & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
      Read (70, Rec=koffset) vkl_, nstsv_, evalsv_
  ! retreive subset
      evalsvp (:) = evalsv_ (:nstsv)
      Deallocate (evalsv_)
#endif
#ifndef XS
      Read (70, Rec=koffset) vkl_, nstsv_, evalsvp
#endif
      Close (70)
!
      t1 = Abs (vkl(1, ik)-vkl_(1)) + Abs (vkl(2, ik)-vkl_(2)) + Abs &
     & (vkl(3, ik)-vkl_(3))
      If (t1 .Gt. input%structure%epslat) Then
         Write (*,*)
         Write (*, '("Error(getevalsv): differing vectors for k-point "&
        &,I8)') ik
         Write (*, '(" current    : ",3G18.10)') vkl (:, ik)
         Write (*, '(" EVALSV.OUT : ",3G18.10)') vkl_
         Write (*, '(" file       : ",a      )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
#ifndef XS
      If (nstsv .Ne. nstsv_) Then
         Write (*,*)
         Write (*, '("Error(getevalsv): differing nstsv for k-point ",I&
        &8)') ik
         Write (*, '(" current    : ",I8)') nstsv
         Write (*, '(" EVALSV.OUT : ",I8)') nstsv_
         Write (*, '(" file       : ",a      )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
#endif
      Return
End Subroutine getevalsv
!EOC

Module m_getevalsvr
      Implicit None
Contains
!
!
      Subroutine getevalsvr (fname, isti, istf, vpl, evalsvp)
         Use modmain
         Implicit None
    ! arguments
         Character (*), Intent (In) :: fname
         Integer, Intent (In) :: isti, istf
         Real (8), Intent (In) :: vpl (3)
         Real (8), Intent (Out) :: evalsvp (:)
    ! local variables
         Integer :: err
         Real (8), Allocatable :: evalsvt (:)
         Character (256) :: tmpstr
    ! check correct shapes
         err = 0
         If ((isti .Lt. 1) .Or. (istf .Gt. nstsv) .Or. (istf .Le. &
        & isti)) Then
            Write (*,*)
            Write (*, '("Error(getevalsvr): inconsistent limits for ban&
           &ds:")')
            Write (*, '(" band limits  : ",2i6)') isti, istf
            Write (*, '(" maximum value: ",i6)') nstsv
            Write (*,*)
            err = err + 1
         End If
         If (size(evalsvp, 1) .Ne. (istf-isti+1)) Then
            Write (*,*)
            Write (*, '("Error(getevalsvr): output array does not match&
           & for bands:")')
            Write (*, '(" band limits              : ",2i6)') isti, &
           & istf
            Write (*, '(" requested number of bands: ",i6)') istf - &
           & isti + 1
            Write (*, '(" array size               : ",i6)') size &
           & (evalsvp, 1)
            Write (*,*)
            err = err + 1
         End If
         If (err .Ne. 0) Stop
         Allocate (evalsvt(nstsv))
         filetag_evalsv = trim (fname)
         tmpstr = trim (filext)
         filext = ''
         Call getevalsv (vpl, evalsvt)
         filetag_evalsv = 'EVALSV'
         filext = trim (tmpstr)
         evalsvp (:) = evalsvt (isti:istf)
         Deallocate (evalsvt)
      End Subroutine getevalsvr
End Module m_getevalsvr
