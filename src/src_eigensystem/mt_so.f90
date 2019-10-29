!
!
!
!
!BOP
! !ROUTINE: mt_nc
! !INTERFACE:
!
!
Subroutine mt_so (pot,basis_alpha,basis_beta,mt_list,rel_level)

! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
! Calculates matrix elements of spin-orbit hamiltonian (and eventually non-collinear magnetic field).
!
! !REVISION HISTORY:
!   Created April 2015 (Andris)
!EOP
!BOC
      Implicit None
      Real(8), intent(in) :: pot(lmmaxvr,nrmtmax,natmtot)
      type(apw_lo_basis_type) :: basis_alpha,basis_beta
      type(MTHamiltonianList) :: mt_list      
      integer, intent(in) :: rel_level


! local variables
      Integer :: is, ia, ias, nr, ir, if1,if3,inonz,ireset1,ireset3,maxi,i
      Integer :: l1, l2, l3, m2, lm2, m1, m3, lm1, lm3
      Integer :: ilo, ilo1, ilo2, io, io1, io2, nalo1, maxnlo
      Real (8) :: t1,t2,angular
!     Real (8), allocatable :: haaintegrals(:,:,:,:,:),halointegrals(:,:,:,:),hlolointegrals(:,:,:)
      Real (8), allocatable :: prefsminus(:,:),prefsplus(:,:)
      Real (8) :: rmtable(nrmtmax),r2inv(nrmtmax),radint(2,apwordmax,apwordmax),rmtable2(nrmtmax),radint2
      complex(8) :: zsum
      Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, nrmtmax),a,alpha
      parameter (alpha=1d0 / 137.03599911d0)
      logical :: Tsymmetric,InterstitialSO
      integer, allocatable :: lfromlm(:),mfromlm(:)
      real(8), allocatable :: rm(:,:),rmlm(:,:)
      logical :: Tfp


InterstitialSO=.false.
if (.true.) then

      l3=input%groundstate%lmaxmat
      allocate(prefsplus(1:l3,-l3:l3))
      allocate(prefsminus(1:l3,-l3:l3))
      prefsplus=0d0
      prefsminus=0d0
      do l1=1,l3
        do m1=-l1,l1-1
          prefsminus(l1,m1)=                  sqrt(dble((l1+1)**2-m1**2)/dble((2*l1+1)*(2*l1+3)))*    sqrt(dble((l1-m1)*(l1-m1+1))/dble(2*(2*l1+1)*(2*l1+3)))
          prefsminus(l1,m1)=prefsminus(l1,m1)+sqrt(dble((l1+1)**2-(m1+1)**2)/dble((2*l1+1)*(2*l1+3)))*sqrt(dble((l1+m1+1)*(l1+m1+2))/dble(2*(2*l1+1)*(2*l1+3)))
          prefsminus(l1,m1)=prefsminus(l1,m1)/sqrt(2d0)
          prefsplus(l1,m1)=                 sqrt(dble(l1**2-m1**2)/dble((2*l1-1)*(2*l1+1)))*    sqrt(dble((l1+m1)*(l1+m1+1))/dble(2*(2*l1-1)*(2*l1+1)))
          prefsplus(l1,m1)=prefsplus(l1,m1)+sqrt(dble(l1**2-(m1+1)**2)/dble((2*l1-1)*(2*l1+1)))*sqrt(dble((l1-m1-1)*(l1-m1))/dble(2*(2*l1-1)*(2*l1+1)))
          prefsplus(l1,m1)=prefsplus(l1,m1)/sqrt(2d0)
        enddo
      enddo

      a=0.5d0*alpha**2


      maxnlo=mt_list%maxnlo

! begin loops over atoms and species
      Do is = 1, nspecies
         nr = nrmt (is)
         Do ir = 1, nr
            r2 (ir) = spr (ir, is) ** 2
            r2inv(ir)=1d0/r2(ir)
         End Do
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)

!if (InterstitialSO) then
            if (rel_level.eq.level_nr) then
              Do ir = 1, nr
                rmtable (ir) = a*pot (1, ir, ias)*y00
              End Do
            elseif (rel_level.eq.level_zora) then
              Do ir = 1, nr
                rmtable (ir) = 1d0/(1d0-a*pot (1, ir, ias)*y00)
              End Do
            else
              Do ir = 1, nr
                rmtable (ir) = a/((1d0-a*pot (1, ir, ias)*y00)**2)
              End Do
            endif
!else
            fr(1:nr)=pot (1, :, ias)
            Call fderiv (1, nr, spr(:, is), fr, gr, cf)
            fr(1:nr)=rmtable (1:nr)
            Call fderiv (1, nr, spr(:, is), fr, gr, cf)
            Do ir = 1, nr
              rmtable2(ir)=gr(ir)/spr(ir,is)
            enddo

           
!endif


!---------------------------!
!     APW-APW integrals     !
!---------------------------!
! Radial integrals first
            if1=apword (0, is)
            Do l1 = 1, input%groundstate%lmaxmat
              Do io1 = 1, apword (l1, is)
                Do io2 = 1, apword (l1, is)
if (InterstitialSO) then
                  Do ir = 1, nr
                    t1=basis_alpha%apwfr(ir, 2, io1, l1, ias)-dble(l1)/spr(ir,is)*basis_alpha%apwfr(ir, 1, io1, l1, ias)
                    t2=basis_beta%apwfr(ir, 2, io2, l1, ias)-dble(l1)/spr(ir,is)*basis_beta%apwfr(ir, 1, io2, l1, ias)
                    fr (ir) =t1*t2*rmtable(ir)*spr(ir,is)**2 
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  radint(1,io1,io2)=gr (nr)

                  Do ir = 1, nr
                    t1=basis_alpha%apwfr(ir, 2, io1, l1, ias)+dble(l1+1)/spr(ir,is)*basis_alpha%apwfr(ir, 1, io1, l1, ias)
                    t2=basis_beta%apwfr(ir, 2, io2, l1, ias)+dble(l1+1)/spr(ir,is)*basis_beta%apwfr(ir, 1, io2, l1, ias)
                    fr (ir) =t1*t2*rmtable(ir)*spr(ir,is)**2
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  radint(2,io1,io2)=gr (nr)
else
                  Do ir = 1, nr
                    fr(ir)=basis_alpha%apwfr(ir, 1, io1, l1, ias)*basis_beta%apwfr(ir, 1, io2, l1, ias)*rmtable2(ir)*spr(ir,is)**2
                  enddo
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  radint(1,io1,io2)=gr(nr)
endif
                End Do
              End Do

              do m1=-l1,l1
                do io1 = 1, apword (l1, is)
                  do io2 = 1, apword (l1, is)
#ifndef SOaa
                    if (m1.ne.l1) then
if (InterstitialSO) then
                       mt_list%ab%aa(if1+(m1+l1)*apword(l1,is)+io1,if1+(m1+l1+1)*apword(l1,is)+io2,ias)=mt_list%ab%aa(if1+(m1+l1)*apword(l1,is)+io1,if1+(m1+l1+1)*apword(l1,is)+io2,ias)+prefsminus(l1,m1)*radint(1,io1,io2)-prefsplus(l1,m1)*radint(2,io1,io2)
else
                       mt_list%ab%aa(if1+(m1+l1)*apword(l1,is)+io1,if1+(m1+l1+1)*apword(l1,is)+io2,ias)=mt_list%ab%aa(if1+(m1+l1)*apword(l1,is)+io1,if1+(m1+l1+1)*apword(l1,is)+io2,ias)+0.5d0*sqrt(dble((l1*(l1+1)-m1*(m1+1))))*radint(1,io1,io2)
endif
                    endif

if (InterstitialSO) then
                    mt_list%alpha%aa(if1+(m1+l1)*apword(l1,is)+io1,if1+(m1+l1)*apword(l1,is)+io2,ias)=mt_list%alpha%aa(if1+(m1+l1)*apword(l1,is)+io1,if1+(m1+l1)*apword(l1,is)+io2,ias)- &
                               0.5d0*(radint(1,io1,io2)*(dble((l1-m1+1)*(l1-m1+2))/dble(2*(2*l1+1)*(2*l1+3))-dble((l1+m1+1)*(l1+m1+2))/dble(2*(2*l1+1)*(2*l1+3)))+ &
                                      radint(2,io1,io2)*(dble((l1+m1-1)*(l1+m1  ))/dble(2*(2*l1-1)*(2*l1+1))-dble((l1-m1-1)*(l1-m1  ))/dble(2*(2*l1-1)*(2*l1+1))))
                    mt_list%beta%aa (if1+(m1+l1)*apword(l1,is)+io1,if1+(m1+l1)*apword(l1,is)+io2,ias)=mt_list%beta%aa (if1+(m1+l1)*apword(l1,is)+io1,if1+(m1+l1)*apword(l1,is)+io2,ias)+ &
                               0.5d0*(radint(1,io1,io2)*(dble((l1-m1+1)*(l1-m1+2))/dble(2*(2*l1+1)*(2*l1+3))-dble((l1+m1+1)*(l1+m1+2))/dble(2*(2*l1+1)*(2*l1+3)))+ &
                                      radint(2,io1,io2)*(dble((l1+m1-1)*(l1+m1  ))/dble(2*(2*l1-1)*(2*l1+1))-dble((l1-m1-1)*(l1-m1  ))/dble(2*(2*l1-1)*(2*l1+1))))
else
                    mt_list%alpha%aa(if1+(m1+l1)*apword(l1,is)+io1,if1+(m1+l1)*apword(l1,is)+io2,ias)= mt_list%alpha%aa(if1+(m1+l1)*apword(l1,is)+io1,if1+(m1+l1)*apword(l1,is)+io2,ias) + 0.5d0*dble(m1)*radint(1,io1,io2)
                    mt_list%beta%aa (if1+(m1+l1)*apword(l1,is)+io1,if1+(m1+l1)*apword(l1,is)+io2,ias)= mt_list%beta%aa (if1+(m1+l1)*apword(l1,is)+io1,if1+(m1+l1)*apword(l1,is)+io2,ias) - 0.5d0*dble(m1)*radint(1,io1,io2)
endif
#endif

                   enddo
                enddo
              enddo

              if1=if1+apword (l1,is)*(2*l1+1)
            End Do

#ifndef SOalo
!--------------------------------------!
!     local-orbital-APW integtrals     !
!--------------------------------------!
            if3=2
            Do ilo = 1, nlorb (is)
              l1 = lorbl (ilo, is)
              if (l1.ne.0) then
                Do io = 1, apword (l1, is)
if (InterstitialSO) then
                  Do ir = 1, nr
                    t1=basis_alpha%apwfr(ir, 2, io, l1, ias)-dble(l1)/spr(ir,is)*basis_alpha%apwfr(ir, 1, io, l1, ias)
                    t2=basis_beta%lofr(ir, 2, ilo, ias)-dble(l1)/spr(ir,is)*basis_beta%lofr(ir, 1, ilo, ias)
                    fr (ir) =t1*t2*rmtable(ir)*spr(ir,is)**2
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  radint(1,io,1)=gr (nr)

                  Do ir = 1, nr
                    t1=basis_alpha%apwfr(ir, 2, io, l1, ias)+dble(l1+1)/spr(ir,is)*basis_alpha%apwfr(ir, 1, io, l1, ias)
                    t2=basis_beta%lofr(ir, 2, ilo, ias)+dble(l1+1)/spr(ir,is)*basis_beta%lofr(ir, 1, ilo, ias)
                    fr (ir) =t1*t2*rmtable(ir)*spr(ir,is)**2
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  radint(2,io,1)=gr (nr)
else
                  Do ir = 1, nr
                    fr(ir)=basis_alpha%apwfr(ir, 1, io, l1, ias)*basis_beta%lofr(ir, 1, ilo, ias)*rmtable2(ir)*spr(ir,is)**2
                  enddo
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  radint(1,io,1)=gr(nr)
endif
                End Do



                if1=0
                do l3=0,l1-1
                  if1=if1+apword(l3,is)*(2*l3+1)
                enddo

                do m1=-l1,l1
                  do io1 = 1, apword (l1, is)
                    if (m1.ne.l1) then
if (InterstitialSO) then
                      mt_list%ab%alo(if1+(m1+l1)*apword(l1,is)+io1,if3+m1+l1,ias)=mt_list%ab%alo(if1+(m1+l1)*apword(l1,is)+io1,if3+m1+l1,ias)+prefsminus(l1,m1)*radint(1,io1,1)-prefsplus(l1,m1)*radint(2,io1,1)
else
                      mt_list%ab%alo(if1+(m1+l1)*apword(l1,is)+io1,if3+m1+l1,ias)=mt_list%ab%alo(if1+(m1+l1)*apword(l1,is)+io1,if3+m1+l1,ias)+0.5d0*sqrt(dble((l1*(l1+1)-m1*(m1+1))))*radint(1,io1,1)
endif
                    endif

                  enddo
                enddo

              endif

              if3=if3+2*l1+1
            End Do

            if3=0
            Do ilo = 1, nlorb (is)
              l1 = lorbl (ilo, is)
              if (l1.ne.0) then
                Do io = 1, apword (l1, is)
if (InterstitialSO) then
                  Do ir = 1, nr
                    t1=basis_alpha%lofr(ir, 2, ilo, ias)-dble(l1)/spr(ir,is)*basis_alpha%lofr(ir, 1, ilo, ias)
                    t2=basis_beta%apwfr(ir, 2, io, l1, ias)-dble(l1)/spr(ir,is)*basis_beta%apwfr(ir, 1, io, l1, ias)
                    fr (ir) =t1*t2*rmtable(ir)*spr(ir,is)**2
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  radint(1,io,1)=gr (nr)

                  Do ir = 1, nr
                    t1=basis_alpha%lofr(ir, 2, ilo, ias)+dble(l1+1)/spr(ir,is)*basis_alpha%lofr(ir, 1, ilo, ias)
                    t2=basis_beta%apwfr(ir, 2, io, l1, ias)+dble(l1+1)/spr(ir,is)*basis_beta%apwfr(ir, 1, io, l1, ias)
                    fr (ir) =t1*t2*rmtable(ir)*spr(ir,is)**2
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  radint(2,io,1)=gr (nr)
else                
                  Do ir = 1, nr
                    fr(ir)=basis_alpha%lofr(ir, 1, ilo, ias)*basis_beta%apwfr(ir, 1, io, l1, ias)*rmtable2(ir)*spr(ir,is)**2
                  enddo
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  radint(1,io,1)=gr(nr)
endif
                End Do

                
                if1=0
                do l3=0,l1-1
                  if1=if1+apword(l3,is)*(2*l3+1)
                enddo

                do m1=-l1,l1
                  do io1 = 1, apword (l1, is)
                    if (m1.ne.l1) then
if (InterstitialSO) then
                      mt_list%ab%loa(if3+m1+l1+1,if1+(m1+l1+1)*apword(l1,is)+io1,ias)=mt_list%ab%loa(if3+m1+l1+1,if1+(m1+l1+1)*apword(l1,is)+io1,ias)+prefsminus(l1,m1)*radint(1,io1,1)-prefsplus(l1,m1)*radint(2,io1,1)
else
                      mt_list%ab%loa(if3+m1+l1+1,if1+(m1+l1+1)*apword(l1,is)+io1,ias)=mt_list%ab%loa(if3+m1+l1+1,if1+(m1+l1+1)*apword(l1,is)+io1,ias)+0.5d0*sqrt(dble((l1*(l1+1)-m1*(m1+1))))*radint(1,io1,1)
endif
                    endif

#ifndef SOloa
if (InterstitialSO) then
                    mt_list%alpha%loa(if3+m1+l1+1,if1+(m1+l1)*apword(l1,is)+io1,ias)=mt_list%alpha%loa(if3+m1+l1+1,if1+(m1+l1)*apword(l1,is)+io1,ias)- &
                               0.5d0*(radint(1,io1,1)*(dble((l1-m1+1)*(l1-m1+2))/dble(2*(2*l1+1)*(2*l1+3))-dble((l1+m1+1)*(l1+m1+2))/dble(2*(2*l1+1)*(2*l1+3)))+ &
                                      radint(2,io1,1)*(dble((l1+m1-1)*(l1+m1  ))/dble(2*(2*l1-1)*(2*l1+1))-dble((l1-m1-1)*(l1-m1  ))/dble(2*(2*l1-1)*(2*l1+1))))
                    mt_list%beta%loa(if3+m1+l1+1,if1+(m1+l1)*apword(l1,is)+io1,ias)=mt_list%beta%loa(if3+m1+l1+1,if1+(m1+l1)*apword(l1,is)+io1,ias)+ &
                               0.5d0*(radint(1,io1,1)*(dble((l1-m1+1)*(l1-m1+2))/dble(2*(2*l1+1)*(2*l1+3))-dble((l1+m1+1)*(l1+m1+2))/dble(2*(2*l1+1)*(2*l1+3)))+ &
                                      radint(2,io1,1)*(dble((l1+m1-1)*(l1+m1  ))/dble(2*(2*l1-1)*(2*l1+1))-dble((l1-m1-1)*(l1-m1  ))/dble(2*(2*l1-1)*(2*l1+1))))
else
                    mt_list%alpha%loa(if3+m1+l1+1,if1+(m1+l1)*apword(l1,is)+io1,ias)= mt_list%alpha%loa(if3+m1+l1+1,if1+(m1+l1)*apword(l1,is)+io1,ias) + 0.5d0*dble(m1)*radint(1,io1,1)
                    mt_list%beta%loa(if3+m1+l1+1,if1+(m1+l1)*apword(l1,is)+io1,ias) = mt_list%beta%loa(if3+m1+l1+1,if1+(m1+l1)*apword(l1,is)+io1,ias)  - 0.5d0*dble(m1)*radint(1,io1,1)

endif

#endif

                  enddo
                enddo

              endif

              if3=if3+2*l1+1
            End Do
#endif


#ifndef SOlolo
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
!            write(*,*) 'lolo'
            if1=0
            Do ilo1 = 1, nlorb (is)
              l1 = lorbl (ilo1, is)
              if3=0
              Do ilo2 = 1, nlorb (is)
                l3 = lorbl (ilo2, is)
                If ((l1 .Eq. l3).and.(l1.ne.0)) Then
if (InterstitialSO) then
                  Do ir = 1, nr
                    t1=basis_alpha%lofr(ir, 2, ilo1, ias)-dble(l1)/spr(ir,is)*basis_alpha%lofr(ir, 1, ilo1, ias)
                    t2=basis_beta%lofr (ir, 2, ilo2, ias)-dble(l1)/spr(ir,is)*basis_beta%lofr (ir, 1, ilo2, ias)
                    fr (ir) =t1*t2*rmtable(ir)*spr(ir,is)**2
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  radint(1,1,1)=gr (nr)
                  Do ir = 1, nr
                    t1=basis_alpha%lofr(ir, 2, ilo1, ias)+dble(l1+1)/spr(ir,is)*basis_alpha%lofr(ir, 1, ilo1, ias)
                    t2=basis_beta%lofr (ir, 2, ilo2, ias)+dble(l1+1)/spr(ir,is)*basis_beta%lofr (ir, 1, ilo2, ias)
                    fr (ir) =t1*t2*rmtable(ir)*spr(ir,is)**2
                  End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  radint(2,1,1)=gr (nr)
else
                  Do ir = 1, nr
                    fr(ir)=basis_alpha%lofr(ir, 1, ilo1, ias)*basis_beta%lofr (ir, 1, ilo2, ias)*rmtable2(ir)*spr(ir,is)**2
                  enddo
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  radint2=gr(nr)
endif

                  do m1=-l1,l1
if (InterstitialSO) then
                    if (m1.ne.l1) then
                      mt_list%ab%lolo(if1+m1+l1+1,if3+m1+l1+2,ias)=mt_list%ab%lolo(if1+m1+l1+1,if3+m1+l1+2,ias)+prefsminus(l1,m1)*radint(1,1,1)-prefsplus(l1,m1)*radint(2,1,1)
                    endif
else
                    if (m1.ne.l1) then
                      mt_list%ab%lolo(if1+m1+l1+1,if3+m1+l1+2,ias)=mt_list%ab%lolo(if1+m1+l1+1,if3+m1+l1+2,ias)+0.5d0*sqrt(dble((l1*(l1+1)-m1*(m1+1))))*radint2 !(1,1,1)
                    endif
endif

#ifndef SOlolo
if (InterstitialSO) then
                    mt_list%alpha%lolo(if1+m1+l1+1,if3+m1+l1+1,ias)=mt_list%alpha%lolo(if1+m1+l1+1,if3+m1+l1+1,ias)- &
                                0.5d0*(radint(1,1,1)*(dble((l1-m1+1)*(l1-m1+2))/dble(2*(2*l1+1)*(2*l1+3))-dble((l1+m1+1)*(l1+m1+2))/dble(2*(2*l1+1)*(2*l1+3)))+ &
                                       radint(2,1,1)*(dble((l1+m1-1)*(l1+m1  ))/dble(2*(2*l1-1)*(2*l1+1))-dble((l1-m1-1)*(l1-m1  ))/dble(2*(2*l1-1)*(2*l1+1))))
                    mt_list%beta%lolo (if1+m1+l1+1,if3+m1+l1+1,ias)=mt_list%beta%lolo (if1+m1+l1+1,if3+m1+l1+1,ias)+ &
                                0.5d0*(radint(1,1,1)*(dble((l1-m1+1)*(l1-m1+2))/dble(2*(2*l1+1)*(2*l1+3))-dble((l1+m1+1)*(l1+m1+2))/dble(2*(2*l1+1)*(2*l1+3)))+ &
                                       radint(2,1,1)*(dble((l1+m1-1)*(l1+m1  ))/dble(2*(2*l1-1)*(2*l1+1))-dble((l1-m1-1)*(l1-m1  ))/dble(2*(2*l1-1)*(2*l1+1))))

else
                    mt_list%alpha%lolo(if1+m1+l1+1,if3+m1+l1+1,ias)=mt_list%alpha%lolo(if1+m1+l1+1,if3+m1+l1+1,ias) + 0.5d0*dble(m1)*radint2 !(1,1,1)
                    mt_list%beta%lolo(if1+m1+l1+1,if3+m1+l1+1,ias) =mt_list%beta%lolo(if1+m1+l1+1,if3+m1+l1+1,ias)  - 0.5d0*dble(m1)*radint2 !(1,1,1)
endif

#endif

                  enddo

                End If
                if3=if3+2*l3+1
              End Do
              if1=if1+2*l1+1
            End Do

#endif

            do if1=1,mt_list%maxaa
              do if3=1,mt_list%maxnlo
                mt_list%alpha%alo(if1,if3,ias)=conjg(mt_list%alpha%loa(if3,if1,ias))
                mt_list%beta%alo(if1,if3,ias)=conjg(mt_list%beta%loa(if3,if1,ias))
              enddo 
            enddo


! end loops over atoms and species
         End Do
      End Do
! cleaning up 



!write (*,*) 'mt_kin ends'

endif
      Return
End Subroutine
!EOC

