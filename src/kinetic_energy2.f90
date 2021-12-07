!
!
!
! Copyright (C) 2014-2021 A. Gulans
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: kinetic_energy
! !INTERFACE:
!
!
Subroutine kinetic_energy2(ik,evecfv,apwalm,ngp,vgpc,igpig)
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   An alternative to the kinetic_energy subroutine. Created for debugging purposes. 
!
! !REVISION HISTORY:
!   Created December 2021 (A. Gulans)
!EOP
!BOC
      Implicit None
      integer, intent (in) :: ngp,ik
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, natmtot)
      Complex (8), Intent (in) :: evecfv (nmatmax, nstfv)
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)

! local variables

      Integer :: is, ia, ias, nr, ir, if1,if3,inonz,ireset1,ireset3,maxi,i,ist
      Integer :: l1, l2, l3, m2, lm2, m1, m3, lm1, lm3
      Integer :: ilo, ilo1, ilo2, io, io1, io2, nalo1, maxnlo
      Real (8) :: t1,t2,angular
      real (8), allocatable :: t_aa(:,:,:),t_alo(:,:),t_lolo(:,:)
      Real (8) :: rmtable(nrmtmax),r2inv(nrmtmax)
      complex(8) :: zsum
      Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax),a,rm,alpha
      parameter (alpha=1d0 / 137.03599911d0)
      complex(8), allocatable :: zm(:,:),zvec(:),zveclo(:),zfft(:),zfft2(:),apwi(:,:),zwf(:,:),zt(:,:),zkn(:,:)
      complex(8) :: zdotc
      external zdotc
      logical :: applyiora

      integer :: ifg,ix,igk,l,m,lm,LOoffset,nfnmax, j1,j3,nbf
      real(8) :: Eiora
      Type(MTHamiltonianList) :: mt_kn
      Type (apw_lo_basis_type) :: mt_basis

write(*,*) 'kinetic starts'
      engyknst(:,ik) = 0.d0

! Initialisation of some variables that exist just for the sake of convenience

      if (input%groundstate%ValenceRelativity.ne.'none') then
        a=0.5d0*alpha**2
      else
        a=0d0
      endif
      applyiora=((input%groundstate%ValenceRelativity.eq.'iora*').or.(input%groundstate%ValenceRelativity.eq.'iora*'))
      if (applyiora) then
        write(*,*) 'Direct kinetic energy calculations with IORA are not implemented yet... '
        stop
! evalfv needs to be incorporated
      endif

! MT part

      mt_basis%lofr=>lofr
      mt_basis%apwfr=>apwfr

      call MTNullify(mt_kn)
      call MTInitAll(mt_kn)
      call MTRedirect(mt_kn%main,mt_kn%spinless)
!write(*,*) '/////////'
!write(*,*) mt_kn%main%lolo
      call mt_kin(veffmt,mt_basis,mt_kn)
!write(*,*) '/////////'
!write(*,*) mt_kn%main%lolo

! APW-APW storage initialisation
!write(*,*) 'maxnlo',mt_kn%maxnlo
!write(*,*) 'maxa',mt_kn%maxaa
      nfnmax=mt_kn%maxaa+mt_kn%maxnlo
      allocate(zwf(nfnmax,nstfv))
      allocate(zm(nfnmax,nstfv))
      allocate(zkn(nstfv,nstfv))
      allocate(zt(nfnmax,nfnmax))
      allocate(apwi(mt_kn%maxaa,ngp))
!write(*,*) 'passed'

! begin loops over atoms and species
      LOoffset=ngp
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
           ias = idxas (ia, is)
          apwi=dcmplx(0d0,0d0)
          if3=0
          Do l3 = 0, input%groundstate%lmaxmat
            Do m3 = - l3, l3
            lm3 = idxlm (l3, m3)
              Do io2 = 1, apword (l3, is)
                if3=if3+1
                apwi(if3,:)=apwalm(1:ngp, io2, lm3, ias)
              End Do
            End Do
          End Do
!---------------------------!
!     APW-APW integrals     !
!---------------------------!

          nbf=if3
          zwf=0d0
          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                     nbf, &          ! M ... rows of op( A ) = rows of C
                     nstfv, &           ! N ... cols of op( B ) = cols of C
                     ngp, &          ! K ... cols of op( A ) = rows of op( B )
                     zone, &          ! alpha
                     apwi, &           ! A
                     mt_kn%maxaa,&           ! LDA ... leading dimension of A
                     evecfv, &           ! B
                     nmatmax, &          ! LDB ... leading dimension of B
                     zzero, &          ! beta
                     zwf, &  ! C
                     nfnmax &      ! LDC ... leading dimension of C
                    )

zt=0d0
          zt(1:nbf,1:nbf)=mt_kn%main%aa(1:nbf,1:nbf,ias)
          if (nlorb(is).ne.0) then
           

           l1=lorbl (1, is)
           lm1=idxlm (l1,-l1)
           j1= ngp + idxlo (lm1, 1, ias)

           l3=lorbl (nlorb(is),is)
           lm3=idxlm (l3,l3)
           j3= ngp + idxlo (lm3,nlorb(is),ias)

!write(*,*) 'j1,j3',j1,j3

           zwf(nbf+1:nbf+j3-j1+1,:)=evecfv(j1:j3,:)           

           zt(nbf+1:nbf+j3-j1+1 , 1:if3)=mt_kn%main%loa(1:j3-j1+1,1:nbf,ias)
           do ilo=1,j3-j1+1
            zt(1:nbf , nbf+ilo)=conjg(mt_kn%main%loa(ilo,1:nbf,ias))
!            zt(1:nbf , nbf+1:nbf+j3-j1+1)=mt_kn%main%alo(1:nbf,1:j3-j1+1,ias)
           enddo
           zt(nbf+1:nbf+j3-j1+1 , nbf+1:nbf+j3-j1+1)=mt_kn%main%lolo(1:j3-j1+1,1:j3-j1+1,ias)
           nbf=nbf+j3-j1+1
          endif

!do ilo=1,nbf
! write(*,*) zwf(ilo,1)
!enddo
!write(*,*) evecfv(ngp+1,1) 
!write(*,*) '/////////'
!write(*,*) mt_kn%main%lolo
!write(*,*) '/////////'
!stop
!          write(*,*) 'nbf,loindex',nbf,j3-j1+1
!          write(*,*) zt(nbf,nbf),mt_kn%main%lolo(1,1,ias)

          call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                     nbf, &          ! M ... rows of op( A ) = rows of C
                     nstfv, &           ! N ... cols of op( B ) = cols of C
                     nbf, &          ! K ... cols of op( A ) = rows of op( B )
                     zone, &          ! alpha
                     zt, &           ! A
                     nfnmax,&           ! LDA ... leading dimension of A
                     zwf, &           ! B
                     nfnmax, &          ! LDB ... leading dimension of B
                     zzero, &          ! beta
                     zm, &  ! C
                     nfnmax &      ! LDC ... leading dimension of C
                    )
!do ilo=1,nbf
! write(*,*) zt(ilo,nbf)
!enddo
!write(*,*) '/////////'
!do ilo=1,nbf
! write(*,*) zm(ilo,1)
!enddo
!stop
          call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                     nstfv, &          ! M ... rows of op( A ) = rows of C
                     nstfv, &           ! N ... cols of op( B ) = cols of C
                     nbf, &          ! K ... cols of op( A ) = rows of op( B )
                     zone, &          ! alpha
                     zwf, &           ! A
                     nfnmax,&           ! LDA ... leading dimension of A
                     zm, &           ! B
                     nfnmax, &          ! LDB ... leading dimension of B
                     zzero, &          ! beta
                     zkn, &  ! C
                     nstfv &      ! LDC ... leading dimension of C
                    )

          do ist=1,nstfv
            engyknst(ist,ik)=engyknst(ist,ik)+zkn(ist,ist)
          enddo
         End Do
      End Do
      deallocate(zm,zwf,zkn,zt,apwi)

      call mt_kn%release()




! Interstitial contribution
      allocate(zfft(ngrtot),zfft2(ngrtot))
      allocate(zvec(ngp))
      do ist=1,nstfv
        do ix=1,3

          zfft=zzero
          Do igk = 1, ngp
            ifg = igfft (igpig(igk))
            zfft (ifg) = evecfv (igk, ist) *vgpc(ix, igk)
          End Do
          Call zfftifc (3, ngrid, 1, zfft)
if  (input%groundstate%ValenceRelativity.eq.'none') then
          Do ir = 1, ngrtot
            zfft (ir)=zfft (ir)*cfunir(ir)
          End Do
else
          if (applyiora) zfft2=zfft
          Do ir = 1, ngrtot
            zfft (ir)=zfft (ir)*cfunir(ir)/(1d0-0.5d0*alpha*alpha*veffir(ir))
          End Do

endif
          Call zfftifc (3, ngrid,-1, zfft)
          zvec=zzero
          do igk=1,ngp
            zvec(igk)=zfft(igfft(igpig(igk)))*vgpc(ix, igk)
          enddo
          engyknst(ist,ik)=engyknst(ist,ik)+0.5d0*zdotc(ngp,evecfv(1,ist),1,zvec,1)

if (applyiora) then
          zfft=zfft2
          Do ir = 1, ngrtot
            zfft (ir)=zfft (ir)*cfunir(ir)/(1d0-0.5d0*alpha*alpha*veffir(ir))**2
          End Do
          Call zfftifc (3, ngrid,-1, zfft)
          zvec=zzero
          do igk=1,ngp
            zvec(igk)=zfft(igfft(igpig(igk)))*vgpc(ix, igk)
          enddo
          engyknst(ist,ik)=engyknst(ist,ik)-0.25d0*a*zdotc(ngp,evecfv(1,ist),1,zvec,1) ! evalfv is missing
endif

        enddo
      enddo

      deallocate(zvec)
      deallocate(zfft,zfft2)

      Return
End Subroutine
!EOC
