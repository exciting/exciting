!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: getevecfv
! !INTERFACE:
!
Subroutine getevecfv (vpl, vgpl, evecfv)
! !USES:
      Use modmain
      Use modinput
      Use modmpi
! !DESCRIPTION:
!   The file where the (first-variational) eigenvectors are stored is
!   {\tt EVECFV.OUT}.
!   It is a direct-access binary file, the record length of which can be
!   determined
!   with the help of the array sizes and data type information.
!   One record of this file corresponds to one k-point in the irreducible
!   Brillouin zone and has the following structure
!
!   \begin{tabular}{|l|l|l|l|l|}
!   \hline
!   $k_{\rm lat}$ & $N_{\rm mat}$ & $N_{\rm stfv}$ & $N_{\rm spfv}$ & $\Phi$ \\
!   \hline
!   \end{tabular}\newline\newline
!   The following table explains the parts of the record in more detail
!
!   \begin{tabular}{|l|l|l|l|}
!   \hline
!   name & type & shape & description\\
!   \hline \hline
!   $k_{\rm lat}$ & real(8) & 3 & k-point in lattice coordinates \\ \hline
!   $N_{\rm mat}$ & integer & 1 & (L)APW basis size including local orbitals \\
!    &  &  & (maximum over k-points) \\ \hline
!   $N_{\rm stfv}$ & integer & 1 & number of (first-variational) states \\
!    &  &  & (without core states) \\ \hline
!   $N_{\rm spfv}$ & integer & 1 & first-variational spins (always equals 1)
!         \\ \hline
!   $\Phi$ & complex(8) & $N_{\rm mat}\times N_{\rm stfv}\times N_{\rm spfv}$ &
!        (first-variational) eigenvector array \\
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
      Real (8), Intent (In) :: vgpl (3, ngkmax)
      Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv, nspnfv)
  ! local variables
      Logical :: exist
      Integer :: isym, lspl, ilo, l, m, lm, koffset
      Integer :: ik, igp, igk, ig, i, ilspl
      Integer :: is, ia, ja, ias, jas
      Integer :: recl, nmatmax_, nstfv_, nspnfv_
      Real (8) :: vkl_ (3), v (3), t1
      Real (8) :: si (3, 3), sc (3, 3)
      Complex (8) zt1
  ! allocatable arrays
#ifdef XS
  ! added feature to access arrays for only a subset of bands
      Complex (8), Allocatable :: evecfv_ (:, :, :)
#endif
      Complex (8), Allocatable :: evecfvt (:, :)
      Complex (8), Allocatable :: zflm1 (:, :), zflm2 (:, :)
      Character (256) :: filetag
      Character (256), External :: outfilenamestring
  ! find the equivalent k-point number and crystal symmetry element
      Call findkpt (vpl, isym, ik)
  ! index to spatial rotation in lattice point group
      lspl = lsplsymc (isym)
  ! find the record length
#ifdef XS
      Inquire (IoLength=Recl) vkl_, nmatmax_, nstfv_, nspnfv_
#endif
#ifndef XS
      Inquire (IoLength=Recl) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
#endif
  !$OMP CRITICAL
      filetag = trim (filetag_evecfv)
      Do i = 1, 100
         Inquire (File=outfilenamestring(filetag, ik), Exist=Exist)
         If (exist) Then
            Open (70, File=outfilenamestring(filetag, ik), Action='READ&
           &', Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
            Exit
         Else
            Call system ('sync')
            Call sleep (5)
            Write (*,*) "Waiting for other process to write" // ":getev&
           &ecfv:" // trim (outfilenamestring(filetag, ik))
         End If
      End Do
      If (splittfile) Then
         koffset = ik - firstk (procofk(ik)) + 1
      Else
         koffset = ik
      End If
#ifdef XS
      Read (70, Rec=1) vkl_, nmatmax_, nstfv_, nspnfv_
      Close (70)
      If (nstfv .Gt. nstfv_) Then
         Write (*,*)
         Write (*, '("Error(getevecfv): invalid nstfv for k-point ", I8&
        &)') ik
         Write (*, '(" current    : ", I8)') nstfv
         Write (*, '(" EVECFV.OUT : ", I8)') nstfv_
         Write (*, '(" file	     : ", a	 )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
      Allocate (evecfv_(nmatmax_, nstfv_, nspnfv_))
      Inquire (IoLength=Recl) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv_
      Open (70, File=outfilenamestring(filetag, ik), Action='READ', &
     & Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
      Read (70, Rec=koffset) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv_
  ! retreive subset
      evecfv (:, :, :) = evecfv_ (:, :nstfv, :)
      Deallocate (evecfv_)
#endif
#ifndef XS
      Read (70, Rec=koffset) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
#endif
      Close (70)
  !$OMP END CRITICAL
      t1 = Abs (vkl(1, ik)-vkl_(1)) + Abs (vkl(2, ik)-vkl_(2)) + Abs &
     & (vkl(3, ik)-vkl_(3))
      If (t1 .Gt. input%structure%epslat) Then
         Write (*,*)
         Write (*, '("Error(getevecfv): differing vectors for k-point "&
        &, I8)') ik
         Write (*, '(" current	  : ", 3G18.10)') vkl (:, ik)
         Write (*, '(" EVECFV.OUT : ", 3G18.10)') vkl_
         Write (*, '(" file	  : ", a      )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
      If (nmatmax .Ne. nmatmax_) Then
         Write (*,*)
         Write (*, '("Error(getevecfv): differing nmatmax for k-point "&
        &, I8)') ik
         Write (*, '(" current	  : ", I8)') nmatmax
         Write (*, '(" EVECFV.OUT : ", I8)') nmatmax_
         Write (*, '(" file	  : ", a      )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
#ifndef XS
      If (nstfv .Ne. nstfv_) Then
         Write (*,*)
         Write (*, '("Error(getevecfv): differing nstfv for k-point ", &
        &I8)') ik
         Write (*, '(" current	  : ", I8)') nstfv
         Write (*, '(" EVECFV.OUT : ", I8)') nstfv_
         Write (*, '(" file	  : ", a      )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
#endif
      If (nspnfv .Ne. nspnfv_) Then
         Write (*,*)
         Write (*, '("Error(getevecfv): differing nspnfv for k-point ",&
        & I8)') ik
         Write (*, '(" current	  : ", I8)') nspnfv
         Write (*, '(" EVECFV.OUT : ", I8)') nspnfv_
         Write (*, '(" file	  : ", a      )') trim &
        & (outfilenamestring(filetag, ik))
         Write (*,*)
         Stop
      End If
  ! if p = k then return
      t1 = Abs (vpl(1)-vkl(1, ik)) + Abs (vpl(2)-vkl(2, ik)) + Abs &
     & (vpl(3)-vkl(3, ik))
      If (t1 .Lt. input%structure%epslat) Return
      If (isspinspiral()) Then
         Write (*,*)
         Write (*, '("Error(getevec): code limitation - cannot rotate s&
        &pin-spiral states")')
         Write (*, '(" (first run one self-consistent loop with no k-po&
        &int reduction)")')
         Write (*,*)
         Stop
      End If
! the inverse of the spatial symmetry rotates k into p
      ilspl = isymlat (lspl)
      si (:, :) = symlat (:, :, ilspl)
!-----------------------------------------------!
!     translate and rotate APW coefficients     !
!-----------------------------------------------!
      Allocate (evecfvt(nmatmax, nstfv))
      Do igk = 1, ngk (1, ik)
         ig = igkig (igk, 1, ik)
         v (:) = dble (ivg(:, ig))
         t1 = - twopi * dot_product (v(:), vtlsymc(:, isym))
         zt1 = cmplx (Cos(t1), Sin(t1), 8)
         evecfvt (igk, :) = zt1 * evecfv (igk, :, 1)
      End Do
      Do igk = 1, ngk (1, ik)
         Call r3mtv (si, vgkl(:, igk, 1, ik), v)
         Do igp = 1, ngk (ik, 1)
            t1 = Abs (v(1)-vgpl(1, igp)) + Abs (v(2)-vgpl(2, igp)) + &
           & Abs (v(3)-vgpl(3, igp))
            If (t1 .Lt. input%structure%epslat) Then
               evecfv (igp, :, 1) = evecfvt (igk, :)
               Go To 10
            End If
         End Do
10       Continue
      End Do
!---------------------------------------------------------!
!     translate and rotate local-orbital coefficients     !
!---------------------------------------------------------!
      If (nlotot .Le. 0) Go To 20
      Allocate (zflm1(lolmmax, nstfv), zflm2(lolmmax, nstfv))
! make a copy of the local-orbital coefficients
      Do i = ngk (1, ik) + 1, nmat (1, ik)
         evecfvt (i, :) = evecfv (i, :, 1)
      End Do
! spatial rotation symmetry matrix in Cartesian coordinates
      sc (:, :) = symlatc (:, :, lspl)
! rotate k-point by inverse symmetry matrix
      Call r3mtv (si, vkl(:, ik), v)
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
    ! equivalent atom for this symmetry
            ja = ieqatom (ia, is, isym)
            jas = idxas (ja, is)
    ! phase factor from translation
            t1 = - twopi * dot_product (vkl(:, ik), input%structure%speciesarray(is)%species%atomarray(ja)%atom%coord(:))
            zt1 = cmplx (Cos(t1), Sin(t1), 8)
            t1 = twopi * dot_product (v(:), input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:))
            zt1 = zt1 * cmplx (Cos(t1), Sin(t1), 8)
    ! rotate local orbitals
            Do ilo = 1, nlorb (is)
               l = lorbl (ilo, is)
               zflm1 (:, :) = 0.d0
               Do m = - l, l
                  lm = idxlm (l, m)
                  i = ngk (1, ik) + idxlo (lm, ilo, jas)
                  zflm1 (lm, :) = evecfvt (i, :)
               End Do
               Call rotzflm (sc, l, nstfv, lolmmax, zflm1, zflm2)
               Do m = - l, l
                  lm = idxlm (l, m)
                  i = ngk (1, ik) + idxlo (lm, ilo, ias)
                  evecfv (i, :, 1) = zt1 * zflm2 (lm, :)
               End Do
            End Do
         End Do
      End Do
      Deallocate (zflm1, zflm2)
20    Continue
      Deallocate (evecfvt)
      Return
End Subroutine getevecfv
!EOC


!
Module m_getevecfvr
      Implicit None
Contains
!
!
      Subroutine getevecfvr (fname, isti, istf, vpl, vgpl, evecfv)
         Use modmain
         Implicit None
    ! arguments
         Character (*), Intent (In) :: fname
         Integer, Intent (In) :: isti, istf
         Real (8), Intent (In) :: vpl (3), vgpl (3, ngkmax)
         Complex (8), Intent (Out) :: evecfv (:, :, :)
    ! local variables
         Integer :: err
         Complex (8), Allocatable :: evecfvt (:, :, :)
         Character (256) :: str1
    ! check correct shapes
         err = 0
         If ((isti .Lt. 1) .Or. (istf .Gt. nstfv) .Or. (istf .Le. &
        & isti)) Then
            Write (*,*)
            Write (*, '("Error(getevecfvr): inconsistent limits for ban&
           &ds:")')
            Write (*, '(" band limits  : ", 2i6)') isti, istf
            Write (*, '(" maximum value: ", i6)') nstfv
            Write (*,*)
            err = err + 1
         End If
         If (size(evecfv, 2) .Ne. (istf-isti+1)) Then
            Write (*,*)
            Write (*, '("Error(getevecfvr): output array does not match&
           & for bands:")')
            Write (*, '(" band limits 	     : ", 2i6)') isti, istf
            Write (*, '(" requested number of bands: ", i6)') istf - &
           & isti + 1
            Write (*, '(" array size		     : ", i6)') size (evecfv, 2)
            Write (*,*)
            err = err + 1
         End If
         If (size(evecfv, 1) .Ne. nmatmax) Then
            Write (*,*)
            Write (*, '("Error(getevecfvr): output array does not match&
           & for nmatmax:")')
            Write (*, '(" nmatmax   : ", i6)') nmatmax
            Write (*, '(" array size: ", i6)') size (evecfv, 1)
            Write (*,*)
            err = err + 1
         End If
         If (size(evecfv, 3) .Ne. nspnfv) Then
            Write (*,*)
            Write (*, '("Error(getevecfvr): output array does not match&
           & for nspnfv:")')
            Write (*, '(" nspnfv    : ", i6)') nspnfv
            Write (*, '(" array size: ", i6)') size (evecfv, 3)
            Write (*,*)
            err = err + 1
         End If
         If (err .Ne. 0) Stop
         Allocate (evecfvt(nmatmax, nstfv, nspnfv))
         filetag_evecfv = trim (fname)
         str1 = trim (filext)
         filext = ''
         Call getevecfv (vpl, vgpl, evecfvt)
         filetag_evecfv = 'EVECFV'
         filext = trim (str1)
         evecfv (:, :, :) = evecfvt (:, isti:istf, :)
         Deallocate (evecfvt)
      End Subroutine getevecfvr
End Module m_getevecfvr
