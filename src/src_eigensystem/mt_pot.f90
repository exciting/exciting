!
!
!
!
! Copyright (C) 2013 exciting team
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: mt_pot
! !INTERFACE:
!
!
Subroutine mt_pot(pot,basis,mt_h)
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Calculates the potential energy contribution to the muffin-tin Hamiltonian.
!
!EOP
!BOC
      Implicit None
      Real(8), intent(in) :: pot(lmmaxvr,nrmtmax,natmtot)
      type(apw_lo_basis_type) :: basis
      Type (MTHamiltonianList) :: mt_h

! local variables
      Integer :: is, ia, ias, nr, ir, if1,if3,inonz,ireset1,ireset3
      Integer :: l1, l2, l3, m2, lm2, m1, m3, lm1, lm3
      Integer :: ilo, ilo1, ilo2, io, io1, io2, maxnlo, maxaa
      Real (8) :: t1,t2,angular
      Real (8), allocatable :: haaintegrals(:,:,:,:,:),halointegrals(:,:,:,:),hlolointegrals(:,:,:)
      Real (8) :: rmtable(nrmtmax),r2inv(nrmtmax)
      complex(8) :: zsum
! automatic arrays
      Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax),a,rm,alpha
      parameter (alpha=1d0 / 137.03599911d0)
      integer, allocatable :: lfromlm(:),mfromlm(:)

! Initialisation of some variables that exist just for the sake of convenience    

      allocate (lfromlm(lmmaxvr))
      allocate (mfromlm(lmmaxvr))
      Do l1 = 0, input%groundstate%lmaxvr
        Do m1 = - l1, l1
          lm1 = idxlm (l1, m1)
          lfromlm(lm1)=l1
          mfromlm(lm1)=m1
        End Do
      End Do


      if (input%groundstate%ValenceRelativity.ne.'none') then
        a=0.5d0*alpha**2
      else
        a=0d0
      endif

! APW-APW storage initialisation
      haaijSize=0
      Do is = 1, nspecies
        if1=0
        Do l1 = 0, input%groundstate%lmaxmat
          Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do io1 = 1, apword (l1, is)
              if1=if1+1
            End Do
          End Do
        End Do
        if (if1.gt.haaijSize) haaijSize=if1
      Enddo
!      if (allocated(haaij)) deallocate(haaij)
!      allocate(haaij(haaijSize,haaijSize,natmtot))
!      haaij=dcmplx(0d0,0d0)
      Allocate (haaintegrals(lmmaxvr, apwordmax, 0:input%groundstate%lmaxapw, apwordmax, 0:input%groundstate%lmaxmat))
      haaintegrals (:, :, :, :, :)=1d100
! APW-LO storage initialisation
!      if (allocated(haloij)) deallocate(haloij)
!      if (allocated(haloijSize)) deallocate(haloijSize)
!      allocate(haloijSize(nspecies))
!      haloijSize=0
!      maxnlo=0
!      Do is = 1, nspecies
!        ias=idxas (1, is)
!        ilo=nlorb (is)
!        if (ilo.gt.0) then
!          l1 = lorbl (ilo, is)
!          lm1 = idxlm (l1, l1)
!          l3 = lorbl (1, is)
!          lm3 = idxlm (l3, -l3)
!          haloijSize(is)=idxlo (lm1, ilo, ias)- idxlo (lm3, 1, ias)+1
!          if (maxnlo.lt.haloijSize(is)) maxnlo=haloijSize(is)
!        endif
!      Enddo
      maxnlo=mt_h%maxnlo
      if (maxnlo.gt.0) then 
!        allocate(haloij(maxnlo,haaijSize,natmtot))
!        haloij=dcmplx(0d0,0d0)
        Allocate (halointegrals(lmmaxvr, apwordmax, 0:input%groundstate%lmaxmat, nlomax))

! LO-LO storage initialisation
!        if (allocated(hloloij)) deallocate(hloloij)
        allocate(hlolointegrals(lmmaxvr,nlomax,nlomax))
!        allocate(hloloij(maxnlo,maxnlo,natmtot))
!        hloloij=dcmplx(0d0,0d0)
      endif

! begin loops over atoms and species
      Do is = 1, nspecies
         nr = nrmt (is)
         Do ir = 1, nr
            r2 (ir) = spr (ir, is) ** 2
            r2inv(ir)=1d0/r2(ir)
         End Do
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, nr
              rmtable (ir) = 1d0/(1d0-a*pot (1, ir, ias)*y00)
            End Do
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
! Radial integrals first
#ifdef USEOMP
!xOMP PARALLEL DEFAULT(NONE) SHARED(input,apword,lmmaxvr,mfromlm,lfromlm,apwfr,r2,pot,spr,nr,haaintegrals,is,ias,rmtable,r2inv) PRIVATE(lm2,m2,l2,ir,t1,fr,gr,cf,l1,l3,t2,angular,io1,io2)
!$OMP PARALLEL DEFAULT(NONE) SHARED(lorbl,nlorb,input,apword,lmmaxvr,mfromlm,lfromlm,apwfr,lofr,r2,pot,spr,nr,haaintegrals,hlolointegrals,halointegrals,is,ias,rmtable,r2inv) PRIVATE(lm2,m2,l2,ir,t1,fr,gr,cf,l1,l3,t2,angular,io1,io2,ilo1,ilo2,io,ilo)
#endif
            Do l1 = 0, input%groundstate%lmaxmat
               Do io1 = 1, apword (l1, is)
                  Do l3 = 0, input%groundstate%lmaxmat
                     Do io2 = 1, apword (l3, is)
#ifdef USEOMP
!$OMP DO
#endif
                        Do lm2 = 1, lmmaxvr
                          m2 = mfromlm(lm2)
                          l2 = lfromlm(lm2)
                          Do ir = 1, nr
                            t1=apwfr(ir,1,io1,l1,ias)*apwfr(ir,1,io2,l3,ias)*r2(ir)
                            fr (ir) = t1 * pot (lm2, ir, ias)
                          End Do
                            Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                            haaintegrals (lm2, io2, l3, io1, l1)=gr (nr)
                        End Do
#ifdef USEOMP
!$OMP END DO NOWAIT
#endif

                     End Do
                  End Do
               End Do
            End Do
#ifdef USEOMP
!xOMP END PARALLEL
#endif

!--------------------------------------!
!     local-orbital-APW integtrals     !
!--------------------------------------!
#ifdef USEOMP
!xOMP PARALLEL DEFAULT(NONE) SHARED(apword,nlorb,lorbl,rmtable,lmmaxvr,mfromlm,lfromlm,apwfr,lofr,r2,pot,spr,nr,halointegrals,is,ias,input,r2inv) PRIVATE(lm2,m2,l2,ir,t1,t2,fr,gr,cf,ilo,io,l1,l3,angular)
#endif
            Do ilo = 1, nlorb (is)
               l1 = lorbl (ilo, is)
               Do l3 = 0, input%groundstate%lmaxmat
                  Do io = 1, apword (l3, is)
#ifdef USEOMP
!$OMP DO
#endif
                     Do lm2 = 1, lmmaxvr
                       m2 = mfromlm(lm2)
                       l2 = lfromlm(lm2)
                       Do ir = 1, nr
                         t1 = lofr (ir, 1, ilo, ias) * apwfr (ir, 1, io, l3, ias) * r2 (ir)
                         fr (ir) = t1 * pot (lm2, ir, ias)
                       End Do
                       Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                       halointegrals (lm2, io, l3, ilo) = gr (nr)  
                     End Do
#ifdef USEOMP
!$OMP END DO NOWAIT
#endif
                  End Do
               End Do
            End Do
#ifdef USEOMP
!xOMP END PARALLEL
#endif
 
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
            Do ilo1 = 1, nlorb (is)
               l1 = lorbl (ilo1, is)
               Do ilo2 = 1, nlorb (is)
                  l3 = lorbl (ilo2, is)
#ifdef USEOMP
!xOMP PARALLEL DEFAULT(NONE) SHARED(lmmaxvr,mfromlm,lfromlm,lofr,r2,pot,spr,nr,hlolointegrals,ilo1,ilo2,is,ias) PRIVATE(lm2,m2,l2,ir,t1,fr,gr,cf)
!$OMP DO
#endif
                  Do lm2 = 1, lmmaxvr
                    m2 = mfromlm(lm2)
                    l2 = lfromlm(lm2)
                    Do ir = 1, nr
                      t1 = lofr (ir, 1, ilo1, ias) * lofr (ir, 1, ilo2, ias) * r2 (ir)
                      fr (ir) = t1 * pot (lm2, ir, ias)
                    End Do
                    Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                    hlolointegrals (lm2, ilo1, ilo2) = gr (nr)
                  End Do
#ifdef USEOMP
!$OMP END DO
#endif
               End Do
            End Do

#ifdef USEOMP
!$OMP END PARALLEL
#endif
!write(*,*) 'gntnonz',gntnonz(1),y00
!read(*,*)

! Now the angular integrals
      t1 = 0.5d0 * rmt (is) ** 2
      if1=0
      inonz=1
      Do l1 = 0, input%groundstate%lmaxmat
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            ireset1=inonz
            Do io1 = 1, apword (l1, is)
               if1=if1+1
               if3=0
               Do l3 = 0, input%groundstate%lmaxmat
                  Do m3 = - l3, l3
                     lm3 = idxlm (l3, m3)
                        ireset3=inonz
                        Do io2 = 1, apword (l3, is)
                          if3=if3+1
                          zsum = 0.d0
                          do while ((gntnonzlm3(inonz).eq.lm3).and.(gntnonzlm1(inonz).eq.lm1))
                            zsum=zsum+gntnonz(inonz)*haaintegrals (gntnonzlm2(inonz), io2, l3, io1, l1)
                            inonz=inonz+1
                          enddo

                          mt_h%main%aa(if1,if3,ias)=mt_h%main%aa(if1,if3,ias)+zsum

                          if (io2.ne.apword (l3, is)) inonz=ireset3
                        End Do
                  End Do
               End Do
              if (io1.ne.apword (l1, is)) inonz=ireset1
            End Do
         End Do
      End Do

     if1=0
      Do ilo = 1, nlorb (is)
         l1 = lorbl (ilo, is)
         inonz=gntnonzlindex(l1)
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            if1=if1+1
            if3=0
            Do l3 = 0, input%groundstate%lmaxmat
               Do m3 = - l3, l3
                  lm3 = idxlm (l3, m3)
                  ireset3=inonz
                  Do io = 1, apword (l3, is)
                     if3=if3+1
                     zsum = 0.d0
                     do while ((gntnonzlm3(inonz).eq.lm3).and.(gntnonzlm1(inonz).eq.lm1))
                       zsum=zsum+gntnonz(inonz)*halointegrals(gntnonzlm2(inonz),io, l3, ilo)
                       inonz=inonz+1
                     enddo
                     mt_h%main%loa(if1,if3,ias)=mt_h%main%loa(if1,if3,ias)+zsum
                     mt_h%main%alo(if3,if1,ias)=mt_h%main%loa(if1,if3,ias)

                     if (io.ne.apword(l3,is)) inonz=ireset3
                  End Do
               End Do
            End Do
         End Do
      End Do

            if1=0
            Do ilo1 = 1, nlorb (is)
              l1 = lorbl (ilo1, is)
              Do m1 = - l1, l1
                lm1 = idxlm (l1, m1)
                if1=if1+1
                if3=0
                Do ilo2 = 1, nlorb (is)
                  l3 = lorbl (ilo2, is)
                  Do m3 = - l3, l3
                    lm3 = idxlm (l3, m3)
                    if3=if3+1
                    zsum = 0.d0
                    inonz=gntnonzl2index(lm1,lm3)
                    do while ((gntnonzlm3(inonz).eq.lm3).and.(gntnonzlm1(inonz).eq.lm1))
                      zsum=zsum+gntnonz(inonz)*dcmplx(hlolointegrals(gntnonzlm2(inonz),ilo1,ilo2),0d0)
                      inonz=inonz+1
                    enddo
                    mt_h%main%lolo(if1,if3,ias)=mt_h%main%lolo(if1,if3,ias)+zsum
                  End Do
                End Do
              End Do
            End Do

! end loops over atoms and species
         End Do
      End Do
! cleaning up 
      deallocate(haaintegrals)
      deallocate(lfromlm,mfromlm)  
      if (allocated(halointegrals)) deallocate(halointegrals)
      if (allocated(hlolointegrals)) deallocate(hlolointegrals)

      Return
End Subroutine
!EOC
