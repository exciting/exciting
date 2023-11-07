!
!
!
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma, C. Ambrosch-Draxl
! F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine seceqnsv (ik, apwalm, evalfv, evecfv, evecsv)
      Use modinput, only: input, isspinorb, issvlo
      use constants, only: y00, zone, zzero, zi
      Use mod_LDA_LU, only: ldapu, llu, vmatlu, lmmaxlu
      Use mod_Gvector, only: ngrtot, cfunir, igfft, ngrid
      Use mod_Gkvector, only: ngkmax, ngk, igkig, vgkc
      Use mod_atoms, only: natmtot, nspecies, natoms, idxas, spr
      Use mod_muffin_tin, only: lmmaxvr, nrcmtmax, lmmaxapw, nrmtmax,&
                              & nrmt, nrcmt, idxlm, rcmt
      Use mod_potential_and_density, only: bxcmt, veffmt, bxcir, ex_coef, ec_coef, xctype
      Use mod_SHT, only: rbshtvr, zbshtvr, zfshtvr
      Use mod_eigensystem !, only: nmatmax
      Use mod_spin, only: ncmag, nspinor, ndmag
      Use mod_eigenvalue_occupancy, only: nstfv, nstsv, evalsv
      Use mod_APW_LO
      Use mod_misc, only: task
      Use mod_timing, only: timesv
      Use generation_wavefunction, only: generate_basisfunction_secondvariation_MT
      Use svlo, only: get_num_of_basis_funs_sv, construct_H_and_S_in_evecfv_plus_lo_basis
      Use modfvsystem, only: newsystem, evsystem, deletesystem, solvewithlapack

      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, natmtot)
      Real (8), Intent (In) :: evalfv (nstfv)
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv)
      Complex (8), Intent (Out) :: evecsv (nstsv, nstsv)
! local variables
      Complex (8) :: evecfv_tmp (nmatmax, get_num_of_basis_funs_sv())
      Integer :: ispn, jspn, ia, is, ias
      Integer :: ist, jst, i, j, k, l, lm, nm, m, io
      Integer :: ir, irc, igk, ifg
      Integer :: nsc, lwork, info
      Integer :: if3, offset, svlo_offset, losize, b, e
      logical :: realspace
      Integer :: num_of_basis_funs_sv
! fine structure constant
      Real (8), Parameter :: alpha = 1.d0 / 137.03599911d0
! electron g factor
      Real (8), Parameter :: ge = 2.0023193043718d0
      Real (8), Parameter :: ga4 = ge * alpha / 4.d0
      Real (8), Parameter :: a24 = alpha ** 2 / 4.d0
      Real (8) :: rm, t1
      Real (8) :: ts0, ts1, ta,tb,tc,td
! automatic arrays
      Complex (8) zlflm (lmmaxvr, 3)
! allocatable arrays
      Real (8), Allocatable :: bmt (:, :, :)
      Real (8), Allocatable :: bir (:, :)
      Real (8), Allocatable :: vr (:)
      Real (8), Allocatable :: drv (:)
      Real (8), Allocatable :: cf (:, :)
      Real (8), Allocatable :: sor (:)
      Real (8), Allocatable :: rwork (:)
      Real (8), Allocatable :: veffmt_pbe(:,:,:)
      Complex (8), Allocatable :: wfmt1 (:, :, :)
      Complex (8), Allocatable :: wfmt2 (:, :, :)
      Complex (8), Allocatable :: zfft1 (:)
      Complex (8), Allocatable :: zfft2 (:)
      Complex (8), Allocatable :: zv (:, :)
      Complex (8), Allocatable :: work (:)
!     Complex (8), Allocatable :: wfmt3 (:, :)
      Complex (8) :: wfmt3 (lmmaxvr, nrcmtmax),wfmt4 (lmmaxvr, nrcmtmax)
      Complex (8) , Allocatable :: zwf (:,:),apwi(:,:),zhwf(:,:),zhlo(:,:)
      Type (evsystem) :: systemfv, systemsv
      Logical :: packed
! external functions
      Complex (8) zdotc, zfmtinp
      External zdotc, zfmtinp
      ! Type (MTHamiltonianList) :: mt_h
      ! Type (apw_lo_basis_type) :: mt_basis

      num_of_basis_funs_sv = get_num_of_basis_funs_sv()
      
      If (issvlo()) then
         evecfv_tmp(:,:) = zzero
         evecfv_tmp(1:ngk(1, ik), 1:nstfv) = evecfv(1:ngk(1, ik), 1:nstfv)
         do i = 1, nlotot
            evecfv_tmp(ngk(1, ik) + i, nstfv + i) = zone
         end do
      else
         evecfv_tmp(:,:) = evecfv(:,:)
      end If

      If (task==7) then
         if (allocated(veffmt_pbe)) deallocate(veffmt_pbe)
         allocate(veffmt_pbe(lmmaxvr,nrmtmax,natmtot))
         call poteff_soc(veffmt_pbe)
      endif

! spin-unpolarised case
      If (( .Not. associated(input%groundstate%spin)) .And. (ldapu .Eq. 0)) Then
         Do i = 1, nstsv
            evalsv (i, ik) = evalfv (i)
         End Do
         evecsv (:, :) = 0.d0
         Do i = 1, nstsv
            evecsv (i, i) = 1.d0
         End Do
         Return
      End If

! number of spin combinations after application of Hamiltonian
      If (associated(input%groundstate%spin)) Then
         If ((ncmag) .Or. (isspinorb())) Then
            nsc = 3
         Else
            nsc = 2
         End If
      Else
         nsc = 1
      End If

      Call timesec (ts0)
      call timesec(ta)

      Allocate (bmt(lmmaxvr, nrcmtmax, 3))
      Allocate (bir(ngrtot, 3))
      Allocate (vr(nrmtmax))
      Allocate (drv(nrmtmax))
      Allocate (cf(3, nrmtmax))
      Allocate (sor(nrcmtmax))
      Allocate (rwork(3*nstsv))
      Allocate (wfmt1(lmmaxvr, nrcmtmax, num_of_basis_funs_sv))
      Allocate (wfmt2(lmmaxvr, nrcmtmax, nsc))
      lwork = 2 * nstsv
      Allocate (work(lwork))
! zero the second-variational Hamiltonian (stored in the eigenvector array)
      evecsv (:, :) = 0.d0

! Which algorithm are we using?
! True - transform FV wave functions to the real space
! False - work with the original basis and the transform to the LAPW basis
      if (.not.associated(input%groundstate%spin)) then
         realspace = .True.
      else
         realspace = input%groundstate%spin%realspace
      end if

      if (realspace) then
!-------------------------!
!     muffin-tin part     !
!-------------------------!
      Do is = 1, nspecies
!        allocate(wfmt3(lmmaxvr, nrcmt(is)))
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            If (associated(input%groundstate%spin)) Then
! exchange-correlation magnetic field in spherical coordinates
               call timesec(tc)
               If (ncmag) Then
! non-collinear
                  irc = 0
                  Do ir = 1, nrmt (is), input%groundstate%lradstep
                     irc = irc + 1
                     Do i = 1, 3
                        Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, &
                       & rbshtvr, lmmaxvr, bxcmt(:, ir, ias, i), 1, &
                       & 0.d0, bmt(:, irc, i), 1)
                     End Do
                  End Do
               Else
! collinear
                  irc = 0
                  Do ir = 1, nrmt (is), input%groundstate%lradstep
                     irc = irc + 1
                     bmt (:, irc, 1:2) = 0.d0
                     Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                    & lmmaxvr, bxcmt(:, ir, ias, 1), 1, 0.d0, bmt(:, &
                    & irc, 3), 1)
                  End Do
               End If
               call timesec(td)
!               write(*,*) td-tc
               call timesec(tc)
! external muffin-tin magnetic field
               Do irc = 1, nrcmt (is)
                  Do i = 1, 3
                     bmt (:, irc, i) = bmt (:, irc, i) + ga4 * &
                    & (input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(i)+input%groundstate%spin%bfieldc(i))
                  End Do
               End Do
               call timesec(td)
!               write(*,*) td-tc
! spin-orbit radial function
               If (isspinorb()) Then
                  If (task==7) then
                    vr(1:nrmt(is)) = veffmt_pbe(1, 1:nrmt(is), ias) * y00
                  else
                    vr(1:nrmt(is)) = veffmt(1, 1:nrmt(is), ias) * y00
                  endif
                  Call fderiv (1, nrmt(is), spr(:, is), vr, drv, cf)
! spin-orbit coupling prefactor
                  irc = 0
                  Do ir = 1, nrmt (is), input%groundstate%lradstep
                     irc = irc + 1
                     rm = 1.d0 - 0.5d0 * (alpha**2) * vr (ir)
                     sor (irc) = a24 * drv (ir) / (spr(ir, is)*rm**2)
                  End Do
               End If
            End If
  ! compute the first-variational wavefunctions
            call timesec(tc)
            call generate_basisfunction_secondvariation_MT(input%groundstate%lmaxvr,lmmaxvr, ia, is, ngk(1,ik), apwalm, evecfv, wfmt1)
             call timesec(td)
!            write(*,*) td-tc
            call timesec(tc)
! begin loop over states
!xOMP PARALLEL DEFAULT(SHARED) PRIVATE(wfmt3,wfmt2,i,j,wfmt4)
!xOMP DO
            Do jst = 1, num_of_basis_funs_sv
               If (associated(input%groundstate%spin)) Then
                  wfmt4(:, :)=wfmt1(:, :, jst)
                  Call zgemm ('N', 'N', lmmaxvr, &
                                & nrcmt(is), lmmaxvr, zone, zbshtvr, &
                                & lmmaxvr, wfmt4(:,:), lmmaxvr, zzero, &
                                & wfmt3(:, :), lmmaxvr)
                  wfmt4(:, :)=wfmt3(:, :)*bmt (:, :, 3)
                  Call zgemm ('N', 'N', lmmaxvr, &
                                & nrcmt(is), lmmaxvr, zone, zfshtvr, &
                                & lmmaxvr, wfmt4(:,:), lmmaxvr, zzero, &
                                & wfmt2(:, :,1), lmmaxvr)
                  wfmt2 (:, :, 2) = - wfmt2 (:, :, 1)
                  If (nsc .Eq. 3) Then
                     wfmt4(:, :)=wfmt3(:, :) * cmplx (bmt(:, :, &
                          & 1),-bmt(:, :, 2), 8)
                     Call zgemm ('N', 'N', lmmaxvr, &
                          & nrcmt(is), lmmaxvr, zone, zfshtvr, &
                          & lmmaxvr, wfmt4(:,:), lmmaxvr, zzero, &
                          & wfmt2(:, :,3), lmmaxvr)
                  End If

                  If (isspinorb()) Then
                     Do irc = 1, nrcmt (is)
                        Call lopzflm (input%groundstate%lmaxvr, &
                             & wfmt1(:, irc, jst), lmmaxvr, zlflm)
                        t1 = sor (irc)
                        Do lm = 1, lmmaxvr
                           wfmt2 (lm, irc, 1) = wfmt2 (lm, irc, 1) + t1 &
                                & * zlflm (lm, 3)
                           wfmt2 (lm, irc, 2) = wfmt2 (lm, irc, 2) - t1 &
                                & * zlflm (lm, 3)
                           wfmt2 (lm, irc, 3) = wfmt2 (lm, irc, 3) + t1 &
                                & * (zlflm(lm, 1)-zi*zlflm(lm, 2))
                        End Do
                     End Do
                  End If
               Else
                  wfmt2 (:, :, :) = 0.d0
               End If
! apply LDA+U potential if required
               If ((ldapu .Ne. 0) .And. (llu(is) .Ge. 0)) Then
                  l = llu (is)
                  nm = 2 * l + 1
                  lm = idxlm (l,-l)
                  Do k = 1, nsc
                     If (k .Eq. 1) Then
                        ispn = 1
                        jspn = 1
                     Else If (k .Eq. 2) Then
                        ispn = 2
                        jspn = 2
                     Else
                        ispn = 1
                        jspn = 2
                     End If
                     Call zgemm ('N', 'N', nm, nrcmt(is), nm, zone, &
                    & vmatlu(lm, lm, ispn, jspn, ias), lmmaxlu, &
                    & wfmt1(lm, 1, jst), lmmaxvr, zone, wfmt2(lm, 1, &
                    & k), lmmaxvr)
                  End Do
               End If
! second-variational Hamiltonian matrix
               Do ist = 1, num_of_basis_funs_sv
                  Do k = 1, nsc
                     If (k .Eq. 1) Then
                        i = ist
                        j = jst
                     Else If (k .Eq. 2) Then
                        i = ist + num_of_basis_funs_sv
                        j = jst + num_of_basis_funs_sv
                     Else
                        i = ist
                        j = jst + num_of_basis_funs_sv
                     End If
                        evecsv (i, j) = evecsv (i, j) + zfmtinp &
                       & (.True., input%groundstate%lmaxmat, nrcmt(is), &
                       & rcmt(:, is), lmmaxvr, wfmt1(:, :, ist), &
                       & wfmt2(:, :, k))
                  End Do
               End Do
            End Do
!xOMP END DO
!xOMP END PARALLEL
            call timesec(td)
!            write(*,*) td-tc
! end loops over atoms and species
         End Do
!         deallocate(wfmt3)
      End Do
      call timesec(tb)
!      write(*,*) 'sv / MT part',tb-ta


! New algorithm for the second variation follows
      else


! debugging info
!do ist=1,50
!  write(*,*) mt_h%main%aa(ist,ist,1)
!enddo
!stop



        allocate(zwf(mt_hscf%maxaa,num_of_basis_funs_sv))
        allocate(apwi(mt_hscf%maxaa,ngk(1, ik)))
        allocate(zhwf(mt_hscf%maxaa,num_of_basis_funs_sv))
        allocate(zhlo(mt_hscf%maxnlo,num_of_basis_funs_sv))

        offset=ngk(1,ik)
        Do is = 1, nspecies
          Do ia = 1, natoms (is)
!--Hamiltonian--
! APW-APW part
          Call timesec (ts0)
          ias = idxas (ia, is)
          zhlo=zzero
          zhwf=zzero
          zwf=zzero
          apwi=zzero
          if3=0
          Do l = 0, input%groundstate%lmaxmat
            Do m = - l, l
            lm = idxlm (l, m)
              Do io = 1, apword (l, is)
                if3=if3+1
                apwi(if3,:)=apwalm(1:ngk(1, ik), io, lm, ias)
              End Do
            End Do
          End Do
          zwf=zzero

! express the wave functions in terms of the muffin-tin basis
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%maxaa, &          ! M ... rows of op( A ) = rows of C
                      nstfv, &           ! N ... cols of op( B ) = cols of C
                      ngk(1, ik), &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      apwi, &        ! A
                      mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                      evecfv_tmp, &           ! B
                      nmatmax, &          ! LDB ... leading dimension of B
                      zzero, &          ! beta
                      zwf, &  ! C
                      mt_hscf%maxaa &      ! LDC ... leading dimension of C
                      )
! calculate parts of T_SO
! alpha-alpha block
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%maxaa, &          ! M ... rows of op( A ) = rows of C
                      nstfv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%alpha%aa(:,:,ias), &        ! A
                      mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                      zwf, &           ! B
                      mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                      zzero, &          ! beta
                      zhwf, &  ! C
                      mt_hscf%maxaa &      ! LDC ... leading dimension of C
                      )
if (mt_hscf%losize(is).gt.0) then
   if (issvlo()) then
      svlo_offset = nstfv + offset - ngk(1,ik)
      losize = mt_hscf%losize(is)

      zhwf(1:mt_hscf%maxaa, svlo_offset+1:svlo_offset+losize) = &
      zhwf(1:mt_hscf%maxaa, svlo_offset+1:svlo_offset+losize) + &
           mt_hscf%alpha%alo(1:mt_hscf%maxaa, 1:losize, ias)

      evecsv(svlo_offset+1:svlo_offset+losize, &
             svlo_offset+1:svlo_offset+losize) = &
           mt_hscf%alpha%lolo(1:losize, 1:losize, ias)
   else
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%losize(is), &          ! M ... rows of op( A ) = rows of C
                      nstfv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%alpha%loa(:,:,ias), &        ! A
                      mt_hscf%maxnlo,&           ! LDA ... leading dimension of A
                      zwf, &           ! B
                      mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                      zzero, &          ! beta
                      zhlo, &  ! C
                      mt_hscf%maxnlo &      ! LDC ... leading dimension of C
                      )
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%losize(is), &          ! M ... rows of op( A ) = rows of C
                      nstfv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%alpha%lolo(:,:,ias), &        ! A
                      mt_hscf%maxnlo,&           ! LDA ... leading dimension of A
                      evecfv_tmp(offset+1,1), &           ! B
                      nmatmax, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      zhlo, &  ! C
                      mt_hscf%maxnlo &      ! LDC ... leading dimension of C
                      )
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%maxaa, &          ! M ... rows of op( A ) = rows of C
                      nstfv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%alpha%alo(:,:,ias), &        ! A
                      mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                      evecfv_tmp(offset+1,1), &           ! B
                      nmatmax, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      zhwf, &  ! C
                      mt_hscf%maxaa &      ! LDC ... leading dimension of C
                      )

          call zgemm('C', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      nstfv, &          ! M ... rows of op( A ) = rows of C
                      nstfv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      evecfv_tmp(offset+1,1), &        ! A
                      nmatmax,&           ! LDA ... leading dimension of A
                      zhlo, &           ! B
                      mt_hscf%maxnlo, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      evecsv(1,1), &  ! C
                      nstsv &      ! LDC ... leading dimension of C
                      )
       end if ! issvlo
endif
          call zgemm('C', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      num_of_basis_funs_sv, &          ! M ... rows of op( A ) = rows of C
                      num_of_basis_funs_sv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      zwf, &        ! A
                      mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                      zhwf, &           ! B
                      mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      evecsv(1,1), &  ! C
                      nstsv &      ! LDC ... leading dimension of C
                      )

! beta-beta block
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%maxaa, &          ! M ... rows of op( A ) = rows of C
                      num_of_basis_funs_sv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%beta%aa(:,:,ias), &        ! A
                      mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                      zwf, &           ! B
                      mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                      zzero, &          ! beta
                      zhwf, &  ! C
                      mt_hscf%maxaa &      ! LDC ... leading dimension of C
                      )
if (mt_hscf%losize(is).gt.0) then
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%losize(is), &          ! M ... rows of op( A ) = rows of C
                      num_of_basis_funs_sv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%beta%loa(:,:,ias), &        ! A
                      mt_hscf%maxnlo,&           ! LDA ... leading dimension of A
                      zwf, &           ! B
                      mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                      zzero, &          ! beta
                      zhlo, &  ! C
                      mt_hscf%maxnlo &      ! LDC ... leading dimension of C
                      )
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%losize(is), &          ! M ... rows of op( A ) = rows of C
                      num_of_basis_funs_sv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%beta%lolo(:,:,ias), &        ! A
                      mt_hscf%maxnlo,&           ! LDA ... leading dimension of A
                      evecfv_tmp(offset+1,1), &           ! B
                      nmatmax, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      zhlo, &  ! C
                      mt_hscf%maxnlo &      ! LDC ... leading dimension of C
                      )
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%maxaa, &          ! M ... rows of op( A ) = rows of C
                      num_of_basis_funs_sv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%beta%alo(:,:,ias), &        ! A
                      mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                      evecfv_tmp(offset+1,1), &           ! B
                      nmatmax, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      zhwf, &  ! C
                      mt_hscf%maxaa &      ! LDC ... leading dimension of C
                      )

          call zgemm('C', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      num_of_basis_funs_sv, &          ! M ... rows of op( A ) = rows of C
                      num_of_basis_funs_sv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      evecfv_tmp(offset+1,1), &        ! A
                      nmatmax,&           ! LDA ... leading dimension of A
                      zhlo, &           ! B
                      mt_hscf%maxnlo, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      evecsv(num_of_basis_funs_sv+1,num_of_basis_funs_sv+1), &  ! C
                      nstsv &      ! LDC ... leading dimension of C
                      )

endif

          call zgemm('C', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      num_of_basis_funs_sv, &          ! M ... rows of op( A ) = rows of C
                      num_of_basis_funs_sv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      zwf, &        ! A
                      mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                      zhwf, &           ! B
                      mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      evecsv(num_of_basis_funs_sv+1,num_of_basis_funs_sv+1), &  ! C
                      nstsv &      ! LDC ... leading dimension of C
                      )
! alpha-beta block
if (ncmag) then
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%maxaa, &          ! M ... rows of op( A ) = rows of C
                      num_of_basis_funs_sv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%ab%aa(:,:,ias), &        ! A
                      mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                      zwf, &           ! B
                      mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                      zzero, &          ! beta
                      zhwf, &  ! C
                      mt_hscf%maxaa &      ! LDC ... leading dimension of C
                      )
if (mt_hscf%losize(is).gt.0) then
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%losize(is), &          ! M ... rows of op( A ) = rows of C
                      num_of_basis_funs_sv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%ab%loa(:,:,ias), &        ! A
                      mt_hscf%maxnlo,&           ! LDA ... leading dimension of A
                      zwf, &           ! B
                      mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                      zzero, &          ! beta
                      zhlo, &  ! C
                      mt_hscf%maxnlo &      ! LDC ... leading dimension of C
                      )
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%losize(is), &          ! M ... rows of op( A ) = rows of C
                      num_of_basis_funs_sv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%ab%lolo(:,:,ias), &        ! A
                      mt_hscf%maxnlo,&           ! LDA ... leading dimension of A
                      evecfv_tmp(offset+1,1), &           ! B
                      nmatmax, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      zhlo, &  ! C
                      mt_hscf%maxnlo &      ! LDC ... leading dimension of C
                      )
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%maxaa, &          ! M ... rows of op( A ) = rows of C
                      num_of_basis_funs_sv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%ab%alo(:,:,ias), &        ! A
                      mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                      evecfv_tmp(offset+1,1), &           ! B
                      nmatmax, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      zhwf, &  ! C
                      mt_hscf%maxaa &      ! LDC ... leading dimension of C
                      )

          call zgemm('C', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      num_of_basis_funs_sv, &          ! M ... rows of op( A ) = rows of C
                      num_of_basis_funs_sv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      evecfv_tmp(offset+1,1), &        ! A
                      nmatmax,&           ! LDA ... leading dimension of A
                      zhlo, &           ! B
                      mt_hscf%maxnlo, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      evecsv(1,num_of_basis_funs_sv+1), &  ! C
                      nstsv &      ! LDC ... leading dimension of C
                      )

endif

          call zgemm('C', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      num_of_basis_funs_sv, &          ! M ... rows of op( A ) = rows of C
                      num_of_basis_funs_sv, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      zwf, &        ! A
                      mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                      zhwf, &           ! B
                      mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      evecsv(1,num_of_basis_funs_sv+1), &  ! C
                      nstsv &      ! LDC ... leading dimension of C
                      )
endif

        offset=offset+mt_hscf%losize(is)
      enddo
      enddo
!do i=1,mt_h%maxnlo
!  write(*,*) dble(mt_h%alpha%lolo(i,2,1))
!enddo

!write(*,*)

        deallocate(zhwf,zhlo)
        deallocate(apwi)
        deallocate(zwf)

     endif

! Debugging info
!do ist=1,nstsv
!  write(*,*) evecsv(2,ist)
!enddo
!stop
If (.not. issvlo()) then
!---------------------------!
!     interstitial part     !
!---------------------------!
      call timesec(ta)
      If (associated(input%groundstate%spin)) Then
         If (ncmag) Then
! non-collinear
            Do ir = 1, ngrtot
               bir (ir, :) = &
               &  (bxcir(ir,:)+ga4*input%groundstate%spin%bfieldc(:))*cfunir(ir)
            End Do
         Else
! collinear
            Do ir = 1, ngrtot
               bir (ir, 1:2) = 0.d0
               bir (ir, 3) = (bxcir(ir, &
              & 1)+ga4*input%groundstate%spin%bfieldc(3)) * cfunir (ir)
            End Do
         End If
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) SHARED(nstfv,ngk,igfft,igkig,ngrid,ik,evecfv_tmp,evecsv,nsc,bir,ngkmax,ngrtot) PRIVATE(zfft1,zfft2,jst,ist,igk,ifg,zv,k,i,j)
#endif
      Allocate (zfft1(ngrtot))
      Allocate (zfft2(ngrtot))
      Allocate (zv(ngkmax, nsc))
#ifdef USEOMP
!$OMP DO
#endif

         Do jst = 1, nstfv
            zfft1 (:) = 0.d0
            Do igk = 1, ngk (1, ik)
               ifg = igfft (igkig(igk, 1, ik))
               zfft1 (ifg) = evecfv_tmp (igk, jst)
            End Do
! Fourier transform wavefunction to real-space
            Call zfftifc (3, ngrid, 1, zfft1)
! multiply with magnetic field and transform to G-space
            zfft2 (:) = zfft1 (:) * bir (:, 3)
            Call zfftifc (3, ngrid,-1, zfft2)
            Do igk = 1, ngk (1, ik)
               ifg = igfft (igkig(igk, 1, ik))
               zv (igk, 1) = zfft2 (ifg)
               zv (igk, 2) = - zfft2 (ifg)
            End Do
            If (nsc .Eq. 3) Then
               zfft2 (:) = zfft1 (:) * cmplx (bir(:, 1),-bir(:, 2), 8)
               Call zfftifc (3, ngrid,-1, zfft2)
               Do igk = 1, ngk (1, ik)
                  ifg = igfft (igkig(igk, 1, ik))
                  zv (igk, 3) = zfft2 (ifg)
               End Do
            End If
! add to Hamiltonian matrix
            Do ist = 1, nstfv
               Do k = 1, nsc
                  If (k .Eq. 1) Then
                     i = ist
                     j = jst
                  Else If (k .Eq. 2) Then
                     i = ist + nstfv
                     j = jst + nstfv
                  Else
                     i = ist
                     j = jst + nstfv
                  End If
                  If (i .Le. j) Then
                     evecsv (i, j) = evecsv (i, j) + zdotc (ngk(1, ik), &
                    & evecfv_tmp(:, ist), 1, zv(:, k), 1)
                  End If
               End Do
            End Do
         End Do
#ifdef USEOMP
!$OMP END DO
#endif
      deallocate (zfft1,zfft2,zv)
#ifdef USEOMP
!$OMP END PARALLEL
#endif
      End If
   End If

!do i=1,nstsv
!  write(*,*) dble(evecsv(i,nstfv+4))
!enddo
!stop

   If (issvlo()) then
      !--------------------------------------------------
      ! add the first-variational part for svlo and solve
      !--------------------------------------------------
      Call construct_H_and_S_in_evecfv_plus_lo_basis(ik, apwalm, evecfv, igkig, vgkc, ngk, nmat, systemfv)
      packed = .false.
      Call newsystem (systemsv, packed, nstsv)
      systemsv%hamilton%za(:, :) = evecsv(:, :)
      systemsv%overlap%za(:, :) = zzero
      Do ispn = 1, nspinor
         b = (ispn-1)*num_of_basis_funs_sv+1
         e = (ispn-1)*num_of_basis_funs_sv+num_of_basis_funs_sv
         systemsv%hamilton%za(b:e, b:e) = systemsv%hamilton%za(b:e, b:e) + systemfv%hamilton%za(:,:)
         systemsv%overlap%za(b:e, b:e) = systemsv%overlap%za(b:e, b:e) + systemfv%overlap%za(:,:)
      End Do
      Call deletesystem (systemfv)

      ! diagonalise second-variational Hamiltonian
      call timesec(ta)
      Call solvewithlapack(systemsv,nstsv,evecsv,evalsv(1,ik))

   else
      !--------------------------------------------------------------------------------
      ! add the diagonal first-variational part for standard second variation and solve
      !--------------------------------------------------------------------------------
      i = 0
      Do ispn = 1, nspinor
         Do ist = 1, nstfv
            i = i + 1
            evecsv (i, i) = evecsv (i, i) + evalfv (ist)
         End Do
      End Do
! diagonalise second-variational Hamiltonian
      call timesec(ta)
      If (ndmag .Eq. 1) Then
! collinear: block diagonalise H
         Call zheev ('V', 'U', nstfv, evecsv, nstsv, evalsv(:, ik), &
        & work, lwork, rwork, info)
         If (info .Ne. 0) Go To 20
         i = nstfv + 1
         Call zheev ('V', 'U', nstfv, evecsv(i, i), nstsv, evalsv(i, &
        & ik), work, lwork, rwork, info)
         If (info .Ne. 0) Go To 20
         Do i = 1, nstfv
            Do j = 1, nstfv
               evecsv (i, j+nstfv) = 0.d0
               evecsv (i+nstfv, j) = 0.d0
            End Do
         End Do
      Else
! non-collinear or spin-unpolarised: full diagonalisation
         Call zheev ('V', 'U', nstsv, evecsv, nstsv, evalsv(:, ik), &
        & work, lwork, rwork, info)
         If (info .Ne. 0) Go To 20
      End If
   End If
      Deallocate (bmt, bir, vr, drv, cf, sor, rwork)
      Deallocate (wfmt1, wfmt2, work)
      Call timesec (ts1)
      call timesec(tb)
!      write(*,*) 'sv / diagonalization', tb-ta

      timesv = timesv + ts1 - ts0
!$OMP CRITICAL
!!      timesv = timesv + ts1 - ts0
!$OMP END CRITICAL
      Return
20    Continue
      Write (*,*)
      Write (*, '("Error(seceqnsv):&
     & diagonalisation of the second-variational Hamiltonian failed")')
      Write (*, '(" for k-point ", I8)') ik
      Write (*, '(" ZHEEV returned INFO = ", I8)') info
      Write (*,*)
      Stop
End Subroutine
