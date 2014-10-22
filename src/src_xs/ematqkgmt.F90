!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematqkgmt (iq, ik, igq,integrals)
      Use modmain
      Use modinput
      Use modxs
!      Use m_zaxpyc
!      Use m_xszoutpr
!      Use m_xszoutpr3
#ifdef USEOMP
      use omp_lib
#endif
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq, ik, igq
      type (mtints_type) :: integrals
  ! local variables
      Character (*), Parameter :: thisnam = 'ematqkgmt'
      Integer :: is, ia, ias, l1, m1, lm1, l3, m3, lm3, io, io1, io2, &
     & ilo, ilo1, ilo2
      Integer :: lmax1, lmax3, ikt, i, j,zmsize, whichthread
      Complex (8), Allocatable :: zm(:,:)
      Complex (8) :: prefactor
      Real (8) :: cmt0, cmt1, cmt2, cmt3, cmt4

#ifdef USEOMP
      whichthread=omp_get_thread_num()
#endif


      ikt = ik
      lmax1 = input%xs%lmaxapwwf
      lmax3 = lmax1
      zmsize=max(apwmaxsize,lomaxsize)
      allocate(zm(1:istu2-istl2+1,zmsize))

!      xih (:, :) = zzero
!      xiuhloa (:, :) = zzero
!      xiohalo (:, :) = zzero
      xiou (:, :, igq) = zzero

  ! loop over species and atoms
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Call timesec (cmt0)
        !---------------------------!
        !     APW-APW contribution  !
        !---------------------------!
          prefactor=fourpi*conjg(sfacgq(igq, ias, iq))
          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      nst2, &          ! M ... rows of op( A ) = rows of C
                      apwmaxsize, &           ! N ... cols of op( B ) = cols of C
                      apwmaxsize, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      apwcmtfun(1,1,ias), &           ! B
                      nst2, &          ! LDB ... leading dimension of B
                      integrals%aa(1,1,ias), &           ! A
                      apwmaxsize,&           ! LDA ... leading dimension of A
                      zzero, &          ! beta
                      zm, &  ! C
                      nst2 & ! LDC ... leading dimension of C
                     )
          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'T', &           ! TRANSB = 'N'  op( B ) = B.
                      nst1, &    ! M ... rows of op( A ) = rows of C
                      nst2, &           ! N ... cols of op( B ) = cols of C
                      apwmaxsize, &          ! K ... cols of op( A ) = rows of op( B )
                      prefactor, &          ! alpha
                      apwcmtfun0(1,1,ias), &           ! B
                      nst1, &          ! LDB ... leading dimension of B
                      zm, &           ! A
                      nst2,&           ! LDA ... leading dimension of A
                      zone, &          ! beta
                      xiou(1,1,igq), &  ! C
                      nst1 &      ! LDC ... leading dimension of C
                     )
 

            Call timesec (cmt1)
           !--------------------------------------!
           !     local-orbital-APW contribution   !
           !--------------------------------------!
           ! loop over local orbitals
          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      nst2, &          ! M ... rows of op( A ) = rows of C
                      losize(is), &           ! N ... cols of op( B ) = cols of C
                      apwmaxsize, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      apwcmtfun(1,1,ias), &           ! B
                      nst2, &          ! LDB ... leading dimension of B
                      integrals%loa(1,1,ias), &           ! A
                      apwmaxsize,&           ! LDA ... leading dimension of A
                      zzero, &          ! beta
                      zm, &  ! C
                      nst2 & ! LDC ... leading dimension of C
                     )
          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'T', &           ! TRANSB = 'N'  op( B ) = B.
                      nst1, &    ! M ... rows of op( A ) = rows of C
                      nst2, &           ! N ... cols of op( B ) = cols of C
                      losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                      prefactor, &          ! alpha
                      locmtfun0(1,1,ias), &           ! B
                      nst1, &          ! LDB ... leading dimension of B
                      zm, &           ! A
                      nst2,&           ! LDA ... leading dimension of A
                      zone, &          ! beta
                      xiou(1,1,igq), &  ! C
                      nst1 &      ! LDC ... leading dimension of C
                     )

               Call timesec (cmt2)
           !--------------------------------------!
           !     APW-local-orbital contribution   !
           !--------------------------------------!
          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      nst2, &          ! M ... rows of op( A ) = rows of C
                      apwmaxsize, &           ! N ... cols of op( B ) = cols of C
                      losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      locmtfun(1,1,ias), &           ! B
                      nst2, &          ! LDB ... leading dimension of B
                      integrals%alo(1,1,ias), &           ! A
                      lomaxsize,&           ! LDA ... leading dimension of A
                      zzero, &          ! beta
                      zm, &  ! C
                      nst2 & ! LDC ... leading dimension of C
                     )
          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'T', &           ! TRANSB = 'N'  op( B ) = B.
                      nst1, &    ! M ... rows of op( A ) = rows of C
                      nst2, &           ! N ... cols of op( B ) = cols of C
                      apwmaxsize, &          ! K ... cols of op( A ) = rows of op( B )
                      prefactor, &          ! alpha
                      apwcmtfun0(1,1,ias), &           ! B
                      nst1, &          ! LDB ... leading dimension of B
                      zm, &           ! A
                      nst2,&           ! LDA ... leading dimension of A
                      zone, &          ! beta
                      xiou(1,1,igq), &  ! C
                      nst1 &      ! LDC ... leading dimension of C
                     )
 
               Call timesec (cmt3)
           !------------------------------------------------!
           !     local-orbital-local-orbital contribution   !
           !------------------------------------------------!
          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      nst2, &          ! M ... rows of op( A ) = rows of C
                      losize(is), &           ! N ... cols of op( B ) = cols of C
                      losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      locmtfun(1,1,ias), &           ! B
                      nst2, &          ! LDB ... leading dimension of B
                      integrals%lolo(1,1,ias), &           ! A
                      lomaxsize,&           ! LDA ... leading dimension of A
                      zzero, &          ! beta
                      zm, &  ! C
                      nst2 & ! LDC ... leading dimension of C
                     )
          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'T', &           ! TRANSB = 'N'  op( B ) = B.
                      nst1, &    ! M ... rows of op( A ) = rows of C
                      nst2, &           ! N ... cols of op( B ) = cols of C
                      losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                      prefactor, &          ! alpha
                      locmtfun0(1,1,ias), &           ! B
                      nst1, &          ! LDB ... leading dimension of B
                      zm, &           ! A
                      nst2,&           ! LDA ... leading dimension of A
                      zone, &          ! beta
                      xiou(1,1,igq), &  ! C
                      nst1 &      ! LDC ... leading dimension of C
                     )

            Call timesec (cmt4)
           if (whichthread.eq.0) then
            cpumtaa = cpumtaa + cmt1 - cmt0
            cpumtloa = cpumtloa + cmt2 - cmt1
            cpumtalo = cpumtalo + cmt3 - cmt2
            cpumtlolo = cpumtlolo + cmt4 - cmt3
           endif
        ! end loop over species and atoms
         End Do ! ia
      End Do ! is
   deallocate(zm)
End Subroutine ematqkgmt
