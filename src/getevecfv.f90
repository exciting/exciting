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
  Use modmpi
  Use modinput
  !use mod_kpoint, only: vkl
  !use mod_Gkvector, only: ngkmax, vgkl, ngk
  !use mod_eigensystem, only: nmatmax, nmat
  use mod_kpoint, only: vkl_ptr
  use mod_Gkvector, only: ngkmax_ptr, vgkl_ptr, ngk_ptr
  use mod_eigensystem, only: nmatmax_ptr, nmat_ptr
  use mod_eigensystem, only: idxlo
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_spin, only: nspnfv
  use mod_names, only: filetag_evecfv
  use mod_symmetry, only: lsplsymc, isymlat, symlat, symlatc, vtlsymc, ieqatom
  use mod_constants, only: twopi
  use mod_APW_LO, only: nlotot, nlorb, lorbl
  use mod_muffin_tin, only: idxlm
  use mod_atoms, only: natoms, nspecies, idxas
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
!   $N_{\rm spfv}$ & integer & 1 & first-variational spins \\
!    &  &  & (2 for spin-spirals, 1 otherwise)
!         \\ \hline
!   $\Phi$ & complex(8) & $N_{\rm mat}\times N_{\rm stfv}\times N_{\rm spfv}$ &
!        (first-variational) eigenvector array \\
!   \hline
!   \end{tabular}\newline\newline
!
! !REVISION HISTORY:
!   Created Feburary 2007 (JKD)
!   Fixed transformation error, October 2007 (JKD, Anton Kozhevnikov)
!   Documentation added, Dec 2009 (S. Sagmeister)
!   Fixed l.o. rotation, June 2010 (A. Kozhevnikov)
!EOP
!BOC
      Implicit None
  ! arguments
      Real (8), Intent (In) :: vpl(3)
      Real (8), Intent (In) :: vgpl(3,ngkmax_ptr,nspnfv)
      Complex (8), Intent (Out) :: evecfv(nmatmax_ptr,nstfv,nspnfv)
  ! local variables
      Logical :: exist
      integer isym,lspl,ilspl
      integer ispn,ilo,l,lm,i,j,m,m1,lm1
      integer ik,igp,igk
      integer is,ia,ja,ias,jas
      Integer :: recl, nmatmax_, nstfv_, nspnfv_, koffset
      Real (8) :: vkl_ (3), v(3), v1(3)
      Real (8) :: sl(3,3), sc(3,3), t1
      Complex(8) ::  zt1
  ! allocatable arrays
#ifdef XS
  ! added feature to access arrays for only a subset of bands
      Complex (8), Allocatable :: evecfv_ (:, :, :)
#endif
      Complex (8), Allocatable :: evecfvt (:, :)
      Character (256) :: filetag
      Character (256), External :: outfilenamestring
      complex(8), external :: getdlmm

  ! find the equivalent k-point number and crystal symmetry element
      Call findkpt(vpl, isym, ik)
      
  ! find the record length
#ifdef XS
      Inquire (IoLength=Recl) vkl_, nmatmax_, nstfv_, nspnfv_
#endif
#ifndef XS
      Inquire (IoLength=Recl) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
#endif
  !$OMP CRITICAL
      filetag = trim (filetag_evecfv)
      Do i = 1, 10
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
         Write (*, '(" file	     : ", a)') trim(outfilenamestring(filetag,ik))
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
      t1 = Abs (vkl_ptr(1, ik)-vkl_(1)) + Abs (vkl_ptr(2, ik)-vkl_(2)) + Abs &
     & (vkl_ptr(3, ik)-vkl_(3))
      If (t1 .Gt. input%structure%epslat) Then
         Write (*,*)
         Write (*, '("Error(getevecfv): differing vectors for k-point "&
        &, I8)') ik
         Write (*, '(" current	  : ", 3G18.10)') vkl_ptr (:, ik)
         Write (*, '(" EVECFV.OUT : ", 3G18.10)') vkl_
         Write (*, '(" file	  : ", a)') trim(outfilenamestring(filetag,ik))
         Write (*,*)
         Stop
      End If
      If (nmatmax_ptr .Ne. nmatmax_) Then
         Write (*,*)
         Write (*, '("Error(getevecfv): differing nmatmax for k-point "&
        &, I8)') ik
         Write (*, '(" current	  : ", I8)') nmatmax_ptr
         Write (*, '(" EVECFV.OUT : ", I8)') nmatmax_
         Write (*, '(" file	  : ", a)') trim(outfilenamestring(filetag,ik))
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
         Write (*, '(" file	  : ", a)') trim(outfilenamestring(filetag,ik))
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
         Write (*, '(" file	  : ", a)') trim(outfilenamestring(filetag,ik))
         Write (*,*)
         Stop
      End If
            
      ! if p = k then return
      t1 = Abs (vpl(1)-vkl_ptr(1, ik)) + &
      &    Abs (vpl(2)-vkl_ptr(2, ik)) + &
      &    Abs (vpl(3)-vkl_ptr(3, ik))
      If (t1 .Lt. input%structure%epslat) Return
      
      ! index to spatial rotation in lattice point group
      lspl = lsplsymc(isym)
      
      ! the inverse of the spatial symmetry rotates k into p
      ilspl = isymlat(lspl)
      sl(:,:) = dble(symlat(:,:,ilspl))
      sc(:,:) = symlatc(:,:,ilspl)
            
!-----------------------------------------------!
!     translate and rotate APW coefficients     !
!-----------------------------------------------!
      allocate(evecfvt(nmatmax_ptr,nstfv))
      evecfvt(:,:) = 0.d0
      
      do ispn = 1, nspnfv
        do igk = 1, ngk_ptr(ispn,ik)
          v(:) = dble(vgkl_ptr(:,igk,ispn,ik))
          t1 = -twopi*dot_product(v(:),vtlsymc(:,isym))
          zt1 = cmplx(cos(t1),sin(t1),8)
          evecfvt(igk,:) = zt1*evecfv(igk,:,ispn)
        end do
        do igk = 1, ngk_ptr(ispn,ik)
          call r3mtv(sl,vgkl_ptr(:,igk,ispn,ik),v)
          do igp = 1, ngk_ptr(ispn,ik)
            t1 = abs(v(1)-vgpl(1,igp,ispn)) &
               + abs(v(2)-vgpl(2,igp,ispn)) &
               + abs(v(3)-vgpl(3,igp,ispn))
            if (t1<input%structure%epslat) then
              evecfv(igp,:,ispn) = evecfvt(igk,:)
              exit ! continue with new igp
            end if
          end do
        end do
      end do
      
!---------------------------------------------------------!
!     translate and rotate local-orbital coefficients     !
!---------------------------------------------------------!

      if (nlotot>0) then
      
        call r3mtv(sl,vkl_ptr(:,ik),v)

! loop over the first-variational spins
        do ispn = 1, nspnfv
        
! make a copy of the local-orbital coefficients
          do i = ngk_ptr(ispn,ik)+1, nmat_ptr(ispn,ik)
            evecfvt(i,:) = evecfv(i,:,ispn)
          end do
          
          do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia,is)
            
            ! equivalent atom for this symmetry
            ja = ieqatom(ia,is,isym)
            jas = idxas(ja,is)
            
            !---------------
            ! phase factor
            !---------------
            
            v1(:) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            t1 = -twopi*dot_product(vkl_ptr(:,ik),v1(:)+vtlsymc(:,isym))
            zt1 = cmplx(cos(t1),sin(t1),8)
            t1 = twopi*dot_product(v(:),input%structure%speciesarray(is)%species%atomarray(ja)%atom%coord(:))
            zt1 = zt1*cmplx(cos(t1),sin(t1),8)
            
            do ilo = 1, nlorb(is)
              l = lorbl(ilo,is)
              do m1 = -l, l
                lm1 = idxlm(l,m1)
                i = ngk_ptr(ispn,ik)+idxlo(lm1,ilo,jas)
                evecfv(i,:,ispn) = 0.d0
                do m = -l, l
                  lm = idxlm(l,m)
                  j = ngk_ptr(ispn,ik)+idxlo(lm,ilo,ias)
                  evecfv(i,:,ispn) = evecfv(i,:,ispn)+ &
                  &                  zt1*evecfvt(j,:)*getdlmm(sc,l,m1,m)
                end do
              end do
            end do ! ilo
            
          end do
          end do
        end do
      end if
      
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
         Real (8), Intent (In) :: vpl (3)
         Real (8), Intent (In) :: vgpl(3,ngkmax_ptr,nspnfv)
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
         If (size(evecfv, 1) .Ne. nmatmax_ptr) Then
            Write (*,*)
            Write (*, '("Error(getevecfvr): output array does not match&
           & for nmatmax:")')
            Write (*, '(" nmatmax   : ", i6)') nmatmax_ptr
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
         Allocate (evecfvt(nmatmax_ptr, nstfv, nspnfv))
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
