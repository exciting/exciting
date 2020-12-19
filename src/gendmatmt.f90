!
!BOP
! !ROUTINE: gendmatmt
! !INTERFACE:
!
!
Subroutine gendmatmt (ik, evecfv, evecsv)
! !USES:
      Use modinput
      Use modmain
      use constants, only: zzero
      use modmpi
! !INPUT/OUTPUT PARAMETERS:
!   ik     : k-point number (in,integer)
!   evecfv : first-variational eigenvectors (in,complex(nmatmax,nstfv,nspnfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
! !DESCRIPTION:
!   Constructs a contribution to the density matrix  due to the $\mathbf{k}$-point ik.
!
! !REVISION HISTORY:
!   Created May 2014 (Andris)
!   Revised Aug 2020 (Ronaldo)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv, nspnfv)
      Complex (8), Intent (In) :: evecsv (nstsv, nstsv)
! local variables
      Integer :: ispn, is, ia, ias
      Integer :: i, n, maxaa, maxnlo, wfsize
      Integer :: l, m, lm, io
      integer , pointer :: losize(:)
      Real (8) :: t1, t2
      Real (8) :: ts0, ts1
      Complex (8) weight
! allocatable arrays
      Complex (8), Allocatable :: apwalm (:, :, :, :, :),apwi(:,:)
      Complex (8), pointer :: wf1(:,:), wf2prime(:,:), wfalpha(:,:),wfbeta(:,:)
      Complex (8), Allocatable :: dm2(:,:)
      integer :: l3,lm3,if3,ngp,l1,lm1,j1,j3

      Call timesec (ts0)

!      ist=1
!      jst=nstsv
!      do while (jst-ist.gt.1)
!        kst=(ist+jst)/2
!        if (occnum(kst).gt.input%groundstate%epsocc) then
!          ist=kst
!        else
!          jst=kst
!        endif
!      enddo
!      nst=ist ! last properly occupied level
!      nst=nstsv

      weight=dcmplx(wkpt(ik),0d0)
      ngp=ngk(1, ik)
      ispn=1
      maxaa=mt_dm%maxaa
      maxnlo=mt_dm%maxnlo
      losize => mt_dm%losize

      wfsize=maxaa+maxnlo
      allocate(dm2(wfsize,wfsize))

      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))

      Do ispn = 1, nspnfv
         Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, ispn, ik), sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, ispn))
      End Do

     allocate(wf1(wfsize,nstfv))
     allocate(wf2prime(wfsize,nstsv))
     allocate(wfalpha(wfsize,nstsv))
     allocate(wfbeta(wfsize,nstsv))
!     wf1=zzero

     allocate(apwi(maxaa,ngk(1,ik)))
!---------------------------------!
!     muffin-tin density matrix   !
!---------------------------------!
      Do is = 1, nspecies
        n = lmmaxvr * nrcmt (is)
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
wf1=zzero

call timesec(t1)
          if3=0
          Do l = 0, input%groundstate%lmaxmat
            Do io = 1, apword (l, is)
              Do m = - l, l
                lm = idxlm (l, m)
                if3=if3+1
                apwi(if3,:)=apwalm(1:ngk(1, ik), io, lm, ias,1)
              End Do
            End Do
          End Do

call timesec(t2)

!     APW coefficients
call timesec(t1)

          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      maxaa, &          ! M ... rows of op( A ) = rows of C
                      nstfv, &           ! N ... cols of op( B ) = cols of C
                      ngk(1,ik), &        ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      apwi, &        ! A
                      maxaa,&           ! LDA ... leading dimension of A
                      evecfv, &           ! B
                      nmatmax, &          ! LDB ... leading dimension of B
                      zzero, &          ! beta
                      wf1, &  ! C
                      wfsize &      ! LDC ... leading dimension of C
                      )
!    LO coefficients
       if (losize(is).gt.0) then
         l1 = lorbl (1, is)
         l3 = lorbl (nlorb(is), is)
         lm1=idxlm (l1,-l1)
         lm3=idxlm (l3,l3)
         j1= ngp + idxlo (lm1, 1, ias)
         j3= ngp + idxlo (lm3, nlorb(is), ias)
         wf1(maxaa+1:maxaa+losize(is),1:nstfv)=evecfv(j1:j3,1:nstfv,1)
       endif


          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      wfsize, &          ! M ... rows of op( A ) = rows of C
                      nstsv, &           ! N ... cols of op( B ) = cols of C
                      nstfv, &        ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      wf1, &        ! A
                      wfsize,&           ! LDA ... leading dimension of A
                      evecsv(1,1), &           ! B
                      nstsv, &          ! LDB ... leading dimension of B
                      zzero, &          ! beta
                      wfalpha, &  ! C
                      wfsize &      ! LDC ... leading dimension of C
                      )

call timesec(t2)

          do i=1,nstsv
             wf2prime(:,i)=wfalpha(:,i)*occsv(i,ik)
          enddo


call timesec(t1)

! Density matrix
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'C', &           ! TRANSB = 'N'  op( B ) = B.
                      wfsize, &          ! M ... rows of op( A ) = rows of C
                      wfsize, &    ! N ... cols of op( B ) = cols of C
                      nstsv, &        ! K ... cols of op( A ) = rows of op( B )
                      weight, &          ! alpha
                      wfalpha, &        ! A
                      wfsize,&           ! LDA ... leading dimension of A
                      wf2prime, &           ! B
                      wfsize, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      mt_dm%alpha%ff(1,1,ias), &  ! C
                      wfsize &      ! LDC ... leading dimension of C
                      )

         if (associated(input%groundstate%spin)) then
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      wfsize, &          ! M ... rows of op( A ) = rows of C
                      nstsv, &           ! N ... cols of op( B ) = cols of C
                      nstfv, &        ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      wf1, &        ! A
                      wfsize,&           ! LDA ... leading dimension of A
                      evecsv(nstfv+1,1), &           ! B
                      nstsv, &          ! LDB ... leading dimension of B
                      zzero, &          ! beta
                      wfbeta, &  ! C
                      wfsize &      ! LDC ... leading dimension of C
                      )

          do i=1,nstsv
             wf2prime(:,i)=wfbeta(:,i)*occsv(i,ik)
          enddo


          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'C', &           ! TRANSB = 'N'  op( B ) = B.
                      wfsize, &          ! M ... rows of op( A ) = rows of C
                      wfsize, &    ! N ... cols of op( B ) = cols of C
                      nstsv, &        ! K ... cols of op( A ) = rows of op( B )
                      weight, &          ! alpha
                      wfbeta, &        ! A
                      wfsize,&           ! LDA ... leading dimension of A
                      wf2prime, &           ! B
                      wfsize, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      mt_dm%beta%ff(1,1,ias), &  ! C
                      wfsize &      ! LDC ... leading dimension of C
                      )

! What we really need are not the densities for spin up and down individually,
! but rather the total density rho=rho_up+rho_down and the magnetisation m=rho_up-rho_down.
! This is the place where we do the transformation.
!          dm2(:,:)=mt_dm%alpha%ff(:,:,ias)-mt_dm%beta%ff(:,:,ias)
!          mt_dm%alpha%ff(:,:,ias)=mt_dm%alpha%ff(:,:,ias)+mt_dm%beta%ff(:,:,ias)
!          mt_dm%beta%ff(:,:,ias)=dm2(:,:)
!          mt_dm%alpha%ff(:,:,ias)=mt_dm%beta%ff(:,:,ias)
!          mt_dm%beta%ff(:,:,ias)=0d0

          if (ncmag) then
           call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                      'C', &           ! TRANSB = 'N'  op( B ) = B.
                       wfsize, &          ! M ... rows of op( A ) = rows of C
                       wfsize, &    ! N ... cols of op( B ) = cols of C
                       nstfv, &        ! K ... cols of op( A ) = rows of op( B )
                       weight, &          ! alpha
                       wfalpha, &        ! A
                       wfsize,&           ! LDA ... leading dimension of A
                       wf2prime, &           ! B
                       wfsize, &          ! LDB ... leading dimension of B
                       zone, &          ! beta
                       mt_dm%ab%ff(1,1,ias), &  ! C
                       wfsize &      ! LDC ... leading dimension of C
                       )
          endif
         endif
call timesec(t2)



         End Do
      End Do


      deallocate(wf1,wf2prime,wfalpha,wfbeta)
      deallocate(apwi,apwalm)
      deallocate(dm2)
      Call timesec (ts1)
      timerho = timerho + ts1 - ts0

      Return
End Subroutine
!EOC
