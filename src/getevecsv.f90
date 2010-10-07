!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
!BOP
! !ROUTINE: getevecsv
! !INTERFACE:
!
Subroutine getevecsv (vpl, evecsv)
! !USES:
      Use modmain
      Use modinput
      Use modmpi
! !DESCRIPTION:
!   The file where the (second-variational) eigenvectors are stored is
!   {\tt EVECSV.OUT}.
!   It is a direct-access binary file, the record length of which can be
!   determined
!   with the help of the array sizes and data type information.
!   One record of this file has the following structure
!
!   \begin{tabular}{|l|l|l|}
!   \hline
!   $k_{\rm lat}$ & $N_{\rm stsv}$ & $\Phi$ \\
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
!   $\Phi$ & complex(8) & $N_{\rm stsv}\times N_{\rm stsv}$ &
!         (second-variational) eigenvector array \\
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
      Complex (8), Intent (Out) :: evecsv (nstsv, nstsv)
! local variables
      Integer :: isym, lspn, ik, ist, i
      Integer :: recl, nstsv_
      Real (8) :: vkl_ (3), det, v (3), th, t1
      Complex (8) su2 (2, 2), zt1, zt2
      Character (256) :: filetag
      Character (256), External :: outfilenamestring
!<chm>
      Logical :: exist
      Integer :: koffset
!</chm>
#ifdef XS
  ! added feature to access arrays for only a subset of bands
      Complex (8), Allocatable :: evecsv_ (:, :)
#endif
! find the k-point number
      Call findkpt (vpl, isym, ik)
! index to global spin rotation in lattice point group
      lspn = lspnsymc (isym)
! find the record length
#ifdef XS
      Inquire (IoLength=Recl) vkl_, nstsv_
#endif
#ifndef XS
      Inquire (IoLength=Recl) vkl_, nstsv_, evecsv
#endif
      filetag = trim (filetag_evecsv)
!$OMP CRITICAL
      Do i = 1, 100
         Inquire (File=outfilenamestring(filetag, ik), Exist=Exist)
         If (exist) Then
            Open (70, File=outfilenamestring(filetag, ik), Action='READ&
           &', Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
            Exit
         Else
            Call system ('sync')
            Write (*,*) "Waiting for other process to write" // ":getev&
           &ecsv:" // trim (outfilenamestring(filetag, ik))
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
         Write (*, '("Error(getevecsv): invalid nstsv for k-point ", I8&
        &)') ik
         Write (*, '(" current    : ", I8)') nstsv
         Write (*, '(" EVECSV.OUT : ", I8)') nstsv_
         Write (*, '(" file	     : ", a	 )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
      Allocate (evecsv_(nstsv_, nstsv_))
      Inquire (IoLength=Recl) vkl_, nstsv_, evecsv_
      Open (70, File=outfilenamestring(filetag, ik), Action='READ', &
     & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
      Read (70, Rec=koffset) vkl_, nstsv_, evecsv_
  ! retreive subset
      evecsv (:, :) = evecsv_ (:nstsv, :nstsv)
      Deallocate (evecsv_)
#endif
#ifndef XS
      Read (70, Rec=koffset) vkl_, nstsv_, evecsv
#endif
      Close (70)
!$OMP END CRITICAL
      t1 = Abs (vkl(1, ik)-vkl_(1)) + Abs (vkl(2, ik)-vkl_(2)) + Abs &
     & (vkl(3, ik)-vkl_(3))
      If (t1 .Gt. input%structure%epslat) Then
         Write (*,*)
         Write (*, '("Error(getevecsv): differing vectors for k-point "&
        &, I8)') ik
         Write (*, '(" current	  : ", 3G18.10)') vkl (:, ik)
         Write (*, '(" EVECSV.OUT : ", 3G18.10)') vkl_
         Write (*, '(" file	  : ", a      )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
#ifndef XS
      If (nstsv .Ne. nstsv_) Then
         Write (*,*)
         Write (*, '("Error(getevecsv): differing nstsv for k-point ", &
        &I8)') ik
         Write (*, '(" current	  : ", I8)') nstsv
         Write (*, '(" EVECSV.OUT : ", I8)') nstsv_
         Write (*, '(" file	  : ", a      )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
#endif
! if symmetry element is the identity return
      If (lspn .Eq. 1) Return
! if eigenvectors are spin-unpolarised return
      If ( .Not. associated(input%groundstate%spin)) Return
! find the SU(2) representation of the spin rotation matrix
      Call rotaxang (input%structure%epslat, symlatc(:, :, lspn), det, &
     & v, th)
      Call axangsu2 (v, th, su2)
! apply SU(2) symmetry matrix to second-variational states
      Do i = 1, nstsv
         Do ist = 1, nstfv
            zt1 = evecsv (ist, i)
            zt2 = evecsv (ist+nstfv, i)
            evecsv (ist, i) = su2 (1, 1) * zt1 + su2 (1, 2) * zt2
            evecsv (ist+nstfv, i) = su2 (2, 1) * zt1 + su2 (2, 2) * zt2
         End Do
      End Do
      Return
End Subroutine
!EOC

Module m_getevecsvr
      Implicit None
Contains
!
!
      Subroutine getevecsvr (fname, isti, istf, vpl, evecsv)
         Use modmain
         Implicit None
    ! arguments
         Character (*), Intent (In) :: fname
         Integer, Intent (In) :: isti, istf
         Real (8), Intent (In) :: vpl (3)
         Complex (8), Intent (Out) :: evecsv (:, :)
    ! local variables
         Integer :: err
         Complex (8), Allocatable :: evecsvt (:, :)
         Character (256) :: strtmp
    ! check correct shapes
         err = 0
         If ((isti .Lt. 1) .Or. (istf .Gt. nstsv) .Or. (istf .Le. &
        & isti)) Then
            Write (*,*)
            Write (*, '("Error(getevecsvr): inconsistent limits for ban&
           &ds:")')
            Write (*, '(" band limits  : ", 2i6)') isti, istf
            Write (*, '(" maximum value: ", i6)') nstsv
            Write (*,*)
            err = err + 1
         End If
         If ((size(evecsv, 1) .Ne. size(evecsv, 2)) .Or. (size(evecsv, &
        & 1) .Ne. (istf-isti+1))) Then
            Write (*,*)
            Write (*, '("Error(getevecsvr): output array does not match&
           & for bands:")')
            Write (*, '(" band limits 	     : ", 2i6)') isti, istf
            Write (*, '(" requested number of bands: ", i6)') istf - &
           & isti + 1
            Write (*, '(" array size		     : ", i6)') size (evecsv, 1)
            Write (*,*)
            err = err + 1
         End If
         If (err .Ne. 0) Stop
         Allocate (evecsvt(nstsv, nstsv))
         filetag_evecsv = trim (fname)
         strtmp = trim (filext)
         filext = ''
         Call getevecsv (vpl, evecsvt)
         filetag_evecsv = 'EVECSV'
         filext = trim (strtmp)
         evecsv (:, :) = evecsvt (isti:istf, isti:istf)
         Deallocate (evecsvt)
      End Subroutine getevecsvr
End Module m_getevecsvr
