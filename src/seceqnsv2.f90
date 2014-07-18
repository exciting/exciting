!
!
!
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma, C. Ambrosch-Draxl
! F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine seceqnsv2 (ik, apwalm, evalfv, evecfv, evecsv)
      Use modmain
      Use modinput
      use modfvsystem
      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Real (8), Intent (In) :: evalfv (nstfv)
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv)
      Complex (8), Intent (Out) :: evecsv (nstsv, nstsv)
! local variables
      Integer :: ispn, jspn, ia, is, ias
      Integer :: ist, jst, i, j, k, l, lm, nm
      Integer :: ir, irc, igk, ifg
      Integer :: nsc, lwork, info
! fine structure constant
      Real (8), Parameter :: alpha = 1.d0 / 137.03599911d0
! electron g factor
      Real (8), Parameter :: ge = 2.0023193043718d0
      Real (8), Parameter :: ga4 = ge * alpha / 4.d0
      Real (8), Parameter :: a24 = alpha ** 2 / 4.d0
      Real (8) :: rm, t1
      Real (8) :: ts0, ts1, ta,tb,tc,td
! automatic arrays
      Complex (8) zftp1 (lmmaxvr), zftp2 (lmmaxvr)
      Complex (8) zlflm (lmmaxvr, 3)
! allocatable arrays
      Type (evsystem) :: system
      Logical :: packed

      Real (8), Allocatable :: bmt (:, :, :),vmtbackup(:, :, :),virbackup(:)
      Real (8), Allocatable :: bir (:, :)
      Real (8), Allocatable :: vr (:)
      Real (8), Allocatable :: drv (:)
      Real (8), Allocatable :: cf (:, :)
      Real (8), Allocatable :: sor (:)
      Real (8), Allocatable :: rwork (:)
      Complex (8), Allocatable :: zm (:, :)
      Complex (8), Allocatable :: wfmt1 (:, :, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :)
      Complex (8), Allocatable :: zfft1 (:)
      Complex (8), Allocatable :: zfft2 (:)
      Complex (8), Allocatable :: zv (:, :)
      Complex (8), Allocatable :: work (:)
!     Complex (8), Allocatable :: wfmt3 (:, :)
      Complex (8) :: wfmt3 (lmmaxvr, nrcmtmax),wfmt4 (lmmaxvr, nrcmtmax)
! external functions
      Complex (8) zdotc, zfmtinp
      External zdotc, zfmtinp
! spin-unpolarised case
      If (( .Not. associated(input%groundstate%spin)) .And. (ldapu .Eq. 0)) Then
         Do i = 1, nstsv
            evalsv (i, ik) = evalfv (i)
         End Do
         evecsv (:, :) = 0.d0
         Do i = 1, nstsv
            evecsv (i, i) = 1.d0
         End Do
!         Return
!      End If
      elseif (.not.ncmag.and..not.isspinorb().and.(ldapu.eq.0)) then
        evecsv (:, :) = 0.d0
        allocate(vmtbackup(lmmaxvr, nrmtmax, natmtot)) 
        allocate(virbackup(ngrtot))
        allocate(zm(nmat(1,ik),nstfv))
        Allocate (rwork(3*nstsv))
        lwork = 2 * nstsv
        Allocate (work(lwork))
        vmtbackup=veffmt
        virbackup=veffir

! Spin up
        veffmt(:,:,:)=vmtbackup(:,:,:)+bxcmt(:,:,:,1)
        veffir(:)=virbackup(:)+bxcir(:,1)
        Do is = 1, nspecies
          Do ia = 1, natoms (is)
            ias=idxas(ia,is)
            veffmt (1, :, ias) = veffmt (1, :, ias) + ga4/y00 * &
              & (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(3)+input%groundstate%spin%bfieldc(3))
          End Do
        End Do
        veffir (:) = veffir (:)+ga4*input%groundstate%spin%bfieldc(3) 


        call hmlint
        packed = input%groundstate%solver%packedmatrixstorage
        
        Call newsystem (system, packed, nmat(1,ik))
        h1on=(input%groundstate%ValenceRelativity.eq.'iora*')
        Call hamiltonandoverlapsetup (system, ngk(1,ik), apwalm, igkig, vgkc)
        call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                    nmat(1,ik), &          ! M ... rows of op( A ) = rows of C
                    nstfv, &           ! N ... cols of op( B ) = cols of C
                    nmat(1,ik), &          ! K ... cols of op( A ) = rows of op( B )
                    zone, &          ! alpha
                    system%hamilton%za, &           ! A
                    nmat(1,ik),&           ! LDA ... leading dimension of A
                    evecfv, &           ! B
                    nmatmax, &          ! LDB ... leading dimension of B
                    zzero, &          ! beta
                    zm, &  ! C
                    nmat(1,ik) &      ! LDC ... leading dimension of C
                   )
        call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                    nstfv, &          ! M ... rows of op( A ) = rows of C
                    nstfv, &           ! N ... cols of op( B ) = cols of C
                    nmat(1,ik), &          ! K ... cols of op( A ) = rows of op( B )
                    zone, &          ! alpha
                    evecfv, &           ! A
                    nmatmax,&           ! LDA ... leading dimension of A
                    zm, &           ! B
                    nmat(1,ik), &          ! LDB ... leading dimension of B
                    zzero, &          ! beta
                    evecsv, &  ! C
                    2*nstfv &      ! LDC ... leading dimension of C
                   )                

        Call zheev ('V', 'U', nstfv, evecsv, nstsv, evalsv(1, ik), work, lwork, rwork, info)

        if (info.ne.0) then
          Write (*,*)
          Write (*, '("Error(seceqnsv): diagonalisation of the second-varia&
            &tional Hamiltonian failed")')
          Write (*, '(" for k-point ", I8)') ik
          Write (*, '(" ZHEEV returned INFO = ", I8)') info
          Write (*,*)
          Stop
        endif

        call deleteystem(system)
! Spin down
        veffmt(:,:,:)=vmtbackup(:,:,:)-bxcmt(:,:,:,1)
        veffir(:)=virbackup(:)-bxcir(:,1)

        Do is = 1, nspecies
          Do ia = 1, natoms (is)
            ias=idxas(ia,is)
            veffmt (1, :, ias) = veffmt (1, :, ias) - ga4/y00 * &
              & (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(3)+input%groundstate%spin%bfieldc(3))
          End Do
        End Do
        veffir (:) = veffir (:)-ga4*input%groundstate%spin%bfieldc(3)

        call hmlint
        packed = input%groundstate%solver%packedmatrixstorage

        Call newsystem (system, packed, nmat(1,ik))
        h1on=(input%groundstate%ValenceRelativity.eq.'iora*')
        Call hamiltonandoverlapsetup (system, ngk(1,ik), apwalm, igkig, vgkc)
        call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                    nmat(1,ik), &          ! M ... rows of op( A ) = rows of C
                    nstfv, &           ! N ... cols of op( B ) = cols of C
                    nmat(1,ik), &          ! K ... cols of op( A ) = rows of op( B )
                    zone, &          ! alpha
                    system%hamilton%za, &           ! A
                    nmat(1,ik),&           ! LDA ... leading dimension of A
                    evecfv, &           ! B
                    nmatmax, &          ! LDB ... leading dimension of B
                    zzero, &          ! beta
                    zm, &  ! C
                    nmat(1,ik) &      ! LDC ... leading dimension of C
                   )
        call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                    nstfv, &          ! M ... rows of op( A ) = rows of C
                    nstfv, &           ! N ... cols of op( B ) = cols of C
                    nmat(1,ik), &          ! K ... cols of op( A ) = rows of op( B )
                    zone, &          ! alpha
                    evecfv, &           ! A
                    nmatmax,&           ! LDA ... leading dimension of A
                    zm, &           ! B
                    nmat(1,ik), &          ! LDB ... leading dimension of B
                    zzero, &          ! beta
                    evecsv(nstfv+1,nstfv+1), &  ! C
                    2*nstfv &      ! LDC ... leading dimension of C
                   )

        Call zheev ('V', 'U', nstfv, evecsv(nstfv+1,nstfv+1), nstsv, evalsv(1+nstfv, ik), work, lwork, rwork, info)

        if (info.ne.0) then
          Write (*,*)
          Write (*, '("Error(seceqnsv): diagonalisation of the second-varia&
            &tional Hamiltonian failed")')
          Write (*, '(" for k-point ", I8)') ik
          Write (*, '(" ZHEEV returned INFO = ", I8)') info
          Write (*,*)
          Stop
        endif

        call deleteystem(system)
        evecsv(1:nstfv, nstfv+1:2*nstfv) = 0.d0
        evecsv(nstfv+1:2*nstfv, 1:nstfv) = 0.d0

        veffmt=vmtbackup
        veffir=virbackup

        deallocate(vmtbackup,virbackup,zm,work,rwork)
      
      else
        write(*,*) 'seceqnsv2 does not support LDA+U and non-collinear calculations'
        stop
      endif

      Return
End Subroutine
