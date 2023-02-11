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
  use mod_kpoint, only: vkl_ptr, nkpt
  use mod_Gkvector, only: ngkmax_ptr, vgkl_ptr, ngk_ptr
  use mod_eigensystem, only: nmatmax_ptr
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_spin, only: nspnfv
  use mod_names, only: filetag_evecfv
  use constants, only: twopi
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
      integer :: isym, ispn, i, ik
      Integer :: recl, nmatmax_, nstfv_, nspnfv_, koffset
      Real (8) :: vkl_ (3), t1
  ! allocatable arrays
#ifdef XS
  ! added feature to access arrays for only a subset of bands
      Complex (8), Allocatable :: evecfv_ (:, :, :)
#endif
      Character (256) :: filetag
      Character (256), External :: outfilenamestring
      complex(8), external :: getdlmm

  ! find the equivalent k-point number and crystal symmetry element
      Call findkpt(vpl, isym, ik)

  ! find the record length
#ifdef XS
      Inquire (IoLength=Recl) vkl_, nmatmax_, nstfv_, nspnfv_
#else
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
         koffset = ik - firstk (procofk(ik, nkpt), nkpt) + 1
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
#else
      Read (70, Rec=koffset) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
#endif
      Close (70)
  !$OMP END CRITICAL
      t1 = Abs (vkl_ptr(1, ik)-vkl_(1)) + &
           Abs (vkl_ptr(2, ik)-vkl_(2)) + &
           Abs (vkl_ptr(3, ik)-vkl_(3))
      If (t1 > 1.d-6) Then
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
      If (t1 < 1.d-6) Return

      do ispn = 1, nspnfv
        call rotate_evecfv( isym, vkl_ptr(:,ik), vpl, ngk_ptr(ispn,ik), vgkl_ptr(:,:,ispn,ik), vgpl(:,:,ispn), &
               evecfv(:,:,ispn), nmatmax_ptr, nstfv)
      end do
      Return
End Subroutine getevecfv
!EOC

!> Rotates the eigenvectors corresponding to wavefunctions \(\psi_{n{\bf p}}\)
!> to the ones corresponing to \(\psi_{n{\bf p}'}\) using the symmetry operation
!> `isym` which rotates \({\bf p}\) into \({\bf p}'\), i.e., 
!> \({\bf p}' = {\bf S}^\top \cdot {\bf p}\). 
subroutine rotate_evecfv( isym, vpl, vprl, ngp, vgpl, vgprl, evecfv, ld, nst)
  use precision, only: dp
  use modinput
  use constants, only: zzero, twopi
  use mod_symmetry, only: lsplsymc, isymlat, symlat, symlatc, vtlsymc, ieqatom
  use mod_APW_LO, only: nlotot, nlorb, lorbl
  use mod_muffin_tin, only: idxlm
  use mod_eigensystem, only: idxlo
  use mod_atoms, only: nspecies, natoms, idxas

  !> global index of symmetry operation
  integer, intent(in) :: isym
  !> k-point \({\bf p}\) in lattice coordinates
  real(dp), intent(in) :: vpl(3)
  !> rotated k-point \({\bf p}' = {\bf S}^\top \cdot {\bf p}\) in lattice coordinates
  real(dp), intent(in) :: vprl(3)
  !> number of \({\bf G+p}\) vectors
  integer, intent(in) :: ngp
  !> \({\bf G+p}\) vectors at \({\bf p}\)
  real(dp), intent(in) :: vgpl(3,*)
  !> \(\bf G+p}'\) vectors at \({bf p}'\)
  real(dp), intent(in) :: vgprl(3,*)
  !> leading dimension of `evecfv`
  integer, intent(in) :: ld
  !> on input: eigenvectors at \({\bf p}\);
  !> on output: eigenvectors at \({\bf p}'\)
  complex(dp), intent(inout) :: evecfv(ld,*)
  !> number of states
  integer, intent(in) :: nst

  integer :: lspl, ilspl, igp, igpr, i, j, is, ia, ja, ias, jas, l, m, m1, lm, lm1, ilo, ist
  real(dp) :: sl(3,3), sc(3,3), v(3), v1(3), t1, t2
  complex(dp) :: zt1

  integer, allocatable :: idxlom(:)
  real(dp), allocatable :: dotp(:)
  complex(dp), allocatable :: evecfvt(:,:), dlmm(:), ztv(:)

  complex(dp), external :: getdlmm, zdotu

  ! index to spatial rotation in lattice point group
  lspl = lsplsymc( isym)
  ! the inverse of the spatial symmetry rotates k into p
  ilspl = isymlat( lspl)
  sl = dble( symlat(:,:,ilspl))
  sc = symlatc(:,:,ilspl)

!-----------------------------------------------!
!     translate and rotate APW coefficients     !
!-----------------------------------------------!
  allocate( evecfvt( ngp+nlotot, nst))
  evecfvt = zzero

  ! get phase factors exp(-i (G+p).tau) with tau being the translation of the symmetry
  allocate( dotp(ngp), ztv(ngp))
  call dgemv( 't', 3, ngp, -twopi, vgpl, 3, vtlsymc(1,isym), 1, 0.d0, dotp, 1)
  ztv = cmplx( cos(dotp), sin(dotp), dp)
!$omp parallel default(shared) private(ist)
!$omp do
  do ist = 1, nst
    evecfvt(1:ngp,ist) = ztv*evecfv(1:ngp,ist)
  end do
!$omp end do
!$omp end parallel
  deallocate( dotp, ztv)

!$omp parallel default(shared) private(igp,igpr,v,ist)
!$omp do
  do igp = 1, ngp
    call r3mtv( sl, vgpl(:,igp), v)
    do igpr = 1, ngp
      if( sum(abs(v-vgprl(:,igpr))) < 1.d-6) then
        do ist = 1, nst
          evecfv( igpr, ist) = evecfvt( igp, ist)
        end do
        exit ! continue with new igp
      end if
    end do
  end do
!$omp end do
!$omp end parallel

!---------------------------------------------------------!
!     translate and rotate local-orbital coefficients     !
!---------------------------------------------------------!
  if (nlotot>0) then
    call r3mtv( sl, vpl, v)
    ! make a copy of the local-orbital coefficients
    evecfvt(ngp+1:ngp+nlotot,1:nst) = evecfv(ngp+1:ngp+nlotot,1:nst)
    do is = 1, nspecies
      do ja = 1, natoms(is)
        jas = idxas(ja,is)
        ! equivalent atom for this symmetry
        ia = ieqatom(ja,is,isym)
        ias = idxas(ia,is)
        !---------------
        ! phase factor
        !---------------
        v1 = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord + vtlsymc(:,isym)
        t1 = dot_product( vpl, v1)
        t2 = dot_product( v, input%structure%speciesarray(is)%species%atomarray(ja)%atom%coord)
        t1 = twopi*(t2 - t1)
        zt1 = cmplx( cos(t1), sin(t1), dp)

        do ilo = 1, nlorb(is)
          l = lorbl(ilo,is)
          ! get index of first basis coefficients of this LO at atom ias
          i = ngp + idxlo( idxlm(l,-l), ilo, ias)
          ! get index of first basis coefficients of this LO at atom jas
          j = ngp + idxlo( idxlm(l,-l), ilo, jas)
          ! get D(l) matrix that rotates spherical harmonics
          dlmm = [((getdlmm(sc,l,m1,m), m1=-l, l), m=-l, l)]
          ! rotate local orbitals
          call zgemm('n', 'n', 2*l+1, nst, 2*l+1, zt1, &
                 dlmm, 2*l+1, &
                 evecfvt(i,1), ngp+nlotot, zzero, &
                 evecfv(j,1), ld)
        end do
      end do
    end do
  end if

  deallocate( evecfvt)
end subroutine


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
