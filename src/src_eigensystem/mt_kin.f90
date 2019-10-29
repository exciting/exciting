!
!
!
!
! Copyright (C) 2013 exciting team
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: mt_kin
! !INTERFACE:
!
!
Subroutine mt_kin(pot,basis,mt_h)
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!  Calculates the potential energy contribution to the muffin-tin Hamiltonian. 
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
      Integer :: ilo, ilo1, ilo2, io, io1, io2, maxnlo
      Real (8) :: t1,t2,angular
      Real (8), allocatable :: haaintegrals(:,:,:,:,:),halointegrals(:,:,:,:),hlolointegrals(:,:,:)
      Real (8) :: rmtable(nrmtmax),r2inv(nrmtmax)
      complex(8) :: zsum
! automatic arrays
      Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, nrmtmax),a,rm,alpha
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
              rmtable (ir) = 1d0/(1d0-a*veffmt (1, ir, ias)*y00)
            End Do
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
! Radial integrals first
!#ifdef USEOMP
!xOMP PARALLEL DEFAULT(NONE) SHARED(lorbl,nlorb,input,apword,lmmaxvr,mfromlm,lfromlm,apwfr,lofr,r2,veffmt,spr,nr,haaintegrals,hlolointegrals,halointegrals,is,ias,rmtable,r2inv) PRIVATE(lm2,m2,l2,ir,t1,fr,gr,cf,l1,l3,t2,angular,io1,io2,ilo1,ilo2,io,ilo)
!#endif
            Do l1 = 0, input%groundstate%lmaxmat
              Do io1 = 1, apword (l1, is)
                Do io2 = 1, apword (l1, is)
                  angular=dble(l1*(l1+1))
                  Do ir = 1, nr
                    t1=basis%apwfr(ir, 1, io1, l1, ias)*basis%apwfr(ir, 1, io2, l1, ias)
                    t2=basis%apwfr(ir, 2, io1, l1, ias)*basis%apwfr(ir, 2, io2, l1, ias)
                    fr (ir) = (0.5d0*t2*rmtable(ir) + 0.5d0*angular*t1*rmtable(ir)*r2inv(ir) )*r2 (ir)
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  haaintegrals (1, io2, l1, io1, l1)= gr (nr) / y00
                End Do
              End Do
            End Do
!--------------------------------------!
!     local-orbital-APW integtrals     !
!--------------------------------------!
            Do ilo = 1, nlorb (is)
              l1 = lorbl (ilo, is)
              Do l3 = 0, input%groundstate%lmaxmat
                If (l1 .Eq. l3) Then
                  Do io = 1, apword (l3, is)
                    angular=dble(l1*(l1+1))
                    Do ir = 1, nr
                       t1=basis%apwfr(ir, 1, io, l1, ias)*basis%lofr(ir, 1, ilo, ias)
                       t2=basis%apwfr(ir, 2, io, l1, ias)*basis%lofr(ir, 2, ilo, ias)
                       fr (ir) = (0.5d0*t2*rmtable(ir) + 0.5d0*angular*t1*rmtable(ir)*r2inv(ir))*r2 (ir)
                    End Do
                    Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                    halointegrals (1, io, l3, ilo) = gr (nr) / y00
                  End Do
                End If
              End Do
            End Do
 
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
            Do ilo1 = 1, nlorb (is)
              l1 = lorbl (ilo1, is)
              Do ilo2 = 1, nlorb (is)
                 l3 = lorbl (ilo2, is)
                 If (l1 .Eq. l3) Then
                   angular=dble(l1*(l1+1))
                   Do ir = 1, nr
                     t1=basis%lofr(ir, 1, ilo1, ias)*basis%lofr(ir, 1, ilo2, ias)
                     t2=basis%lofr(ir, 2, ilo1, ias)*basis%lofr(ir, 2, ilo2, ias)
                     fr (ir) = (0.5d0*t2*rmtable(ir) + 0.5d0*angular*t1*rmtable(ir)*r2inv(ir))*r2 (ir)
                   End Do
                   Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                   hlolointegrals (1, ilo1, ilo2) = gr (nr) / y00
                 End If
              End Do
            End Do

! Now the angular integrals
            t1 = 0.5d0 * rmt (is) ** 2
            if1=0
            inonz=1
            Do l1 = 0, input%groundstate%lmaxmat
              Do m1 = - l1, l1
                lm1 = idxlm (l1, m1)
                Do io1 = 1, apword (l1, is)
                  Do io2 = 1, apword (l1, is)
                    mt_h%main%aa(if1+io1,if1+io2,ias)=mt_h%main%aa(if1+io1,if1+io2,ias)+haaintegrals (1, io2, l1, io1, l1)*y00
                  End Do
                End Do
                if1=if1+apword (l1, is)
              End Do
            End Do


            if1=0
            Do ilo = 1, nlorb (is)
              l1 = lorbl (ilo, is)
              Do m1 = - l1, l1
                lm1 = idxlm (l1, m1)
                if1=if1+1
                if3=0
                Do l3 = 0, input%groundstate%lmaxmat
                  Do m3 = - l3, l3
                    lm3 = idxlm (l3, m3)
                    Do io = 1, apword (l3, is)
                      if3=if3+1
                      if (lm1.eq.lm3) mt_h%main%loa(if1,if3,ias)=mt_h%main%loa(if1,if3,ias)+halointegrals(1,io, l3, ilo)*y00
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
                    if (lm1.eq.lm3) then
                      mt_h%main%lolo(if1,if3,ias)=mt_h%main%lolo(if1,if3,ias)+hlolointegrals(1,ilo1,ilo2)*y00
                    endif
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
End Subroutine mt_kin
!EOC
