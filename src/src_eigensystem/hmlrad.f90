!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: hmlrad
! !INTERFACE:
!
!
Subroutine hmlrad
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Calculates the radial Hamiltonian integrals of the APW and local-orbital
!   basis functions. In other words, for spin $\sigma$ and atom $j$ of species
!   $i$, it computes integrals of the form
!   $$ h^{\sigma;ij}_{qq';ll'l''m''}=\begin{cases}
!    \int_0^{R_i}u^{\sigma;ij}_{q;l}(r)Hu^{\sigma;ij}_{q';l'}(r)r^2dr & l''=0 \\
!    \int_0^{R_i}u^{\sigma;ij}_{q;l}(r)V^{\sigma;ij}_{l''m''}(r)
!    u^{\sigma;ij}_{q';l'}(r)r^2dr & l''>0 \end{cases}, $$
!   where $u^{\sigma;ij}_{q;l}$ is the $q$th APW radial function for angular
!   momentum $l$; $H$ is the Hamiltonian of the radial Schr\"{o}dinger equation;
!   and $V^{\sigma;ij}_{l''m''}$ is the effective muffin-tin potential. Similar
!   integrals are calculated for APW-local-orbital and
!   local-orbital-local-orbital contributions.
!
! !REVISION HISTORY:
!   Created December 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, nr, ir, if1,if3,inonz,ireset1,ireset3
      Integer :: l1, l2, l3, m2, lm2, m1, m3, lm1, lm3
      Integer :: ilo, ilo1, ilo2, io, io1, io2, nalo1, maxnlo
      Real (8) :: t1,t2,angular
      Real (8), allocatable :: hintegrals(:,:,:,:,:),halointegrals(:,:,:,:)
      complex(8) :: zsum
! automatic arrays
      Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax),a,rm,energyref,alpha
      parameter (alpha=1d0 / 137.03599911d0)
! begin loops over atoms and species
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
      if (allocated(haaij)) deallocate(haaij)
      allocate(haaij(haaijSize,haaijSize,natmtot))
      Allocate (hintegrals(lmmaxvr, apwordmax, 0:input%groundstate%lmaxapw, apwordmax, 0:input%groundstate%lmaxmat))
            energyref=input%groundstate%energyref
            if (input%groundstate%ValenceRelativity.ne.'none') then
              a=0.5d0*alpha**2
            else
              a=0d0
            endif

! APW-LO storage initialisation
      if (allocated(haloij)) deallocate(haloij)
      if (allocated(haloijSize)) deallocate(haloijSize)
      allocate(haloijSize(nspecies))
      maxnlo=0
      Do is = 1, nspecies
        ias=idxas (1, is)
        ilo=nlorb (is)
        l1 = lorbl (ilo, is)
        lm1 = idxlm (l1, l1)
        l3 = lorbl (1, is)
        lm3 = idxlm (l3, -l3)
        haloijSize(is)=idxlo (lm1, ilo, ias)- idxlo (lm3, 1, ias)+1
        if (maxnlo.lt.haloijSize(is)) maxnlo=haloijSize(is)
      Enddo
      allocate(haloij(maxnlo,haaijSize,natmtot))
      haloij=dcmplx(0d0,0d0)
      Allocate (halointegrals(lmmaxvr, apwordmax, 0:input%groundstate%lmaxmat, nlomax))

      Do is = 1, nspecies
         nr = nrmt (is)
         Do ir = 1, nr
            r2 (ir) = spr (ir, is) ** 2
         End Do
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
           if (input%groundstate%SymmetricKineticEnergy) then
! If kinetic energy is calculated as nabla*nabla
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
! Radial integrals first
            Do l1 = 0, input%groundstate%lmaxmat
               Do io1 = 1, apword (l1, is)
                  Do l3 = 0, input%groundstate%lmaxmat
                     Do io2 = 1, apword (l3, is)
                        If (l1 .Eq. l3) Then
                           angular=dble(l1*(l1+1))
                           Do ir = 1, nr
                              rm=1d0/(1d0+a*(energyref-veffmt (1, ir, ias)*y00))
                              t1=apwfr(ir, 1, io1, l1, ias)*apwfr(ir, 1, io2, l3, ias)
                              t2=apwfr(ir, 2, io1, l1, ias)*apwfr(ir, 2, io2, l3, ias)
                              fr (ir) = (0.5d0*t2*rm + 0.5d0*angular*t1*rm/spr(ir,is)**2 + t1*veffmt(1, ir, ias)* y00)*r2 (ir)
                           End Do
                           Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                           hintegrals (1, io2, l3, io1, l1)= gr (nr) / y00
! calculate more integrals if linearized Koelling-Harmon is demanded
                           if (input%groundstate%ValenceRelativity.eq.'lkh') then
                             Do ir = 1, nr
                               rm=1d0/(1d0+a*(energyref-veffmt (1, ir, ias)*y00))
                               t1=apwfr(ir, 1, io1, l1, ias)*apwfr(ir, 1, io2, l3, ias)
                               t2=apwfr(ir, 2, io1, l1, ias)*apwfr(ir, 2, io2, l3, ias)
                               fr (ir) = a*(0.5d0*t2*rm**2 + 0.5d0*angular*t1*rm**2/spr(ir,is)**2)*r2 (ir)
                             End Do
                             Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                             h1aa (io1, io2, l1, ias) = gr (nr) / y00
                           endif
                        Else
                           hintegrals (1, io2, l3, io1, l1) = 0.d0
                        End If 
                           Do l2 = 1, input%groundstate%lmaxvr
                              Do m2 = - l2, l2
                                 lm2 = idxlm (l2, m2)
                                 Do ir = 1, nr
                                    t1=apwfr(ir,1,io1,l1,ias)*apwfr(ir,1,io2,l3,ias)*r2(ir)
                                    fr (ir) = t1 * veffmt (lm2, ir, ias)
                                 End Do
                                 Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                                 hintegrals (lm2, io2, l3, io1, l1)=gr (nr)
                              End Do
                           End Do
                     End Do
                  End Do
               End Do
            End Do
! Angular integrals
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
                            zsum=zsum+gntnonz(inonz)*hintegrals (gntnonzlm2(inonz), io2, l3, io1, l1)
                            inonz=inonz+1
                          enddo
                          haaij(if1,if3,ias)=zsum
                          if (io2.ne.apword (l3, is)) inonz=ireset3
                        End Do
                  End Do
               End Do
              if (io1.ne.apword (l1, is)) inonz=ireset1
            End Do
         End Do
      End Do
!--------------------------------------!
!     local-orbital-APW integtrals     !
!--------------------------------------!
! Radial integrals
            Do ilo = 1, nlorb (is)
               l1 = lorbl (ilo, is)
               Do l3 = 0, input%groundstate%lmaxmat
                  Do io = 1, apword (l3, is)
                     If (l1 .Eq. l3) Then
                        angular=dble(l1*(l1+1))
                        Do ir = 1, nr
                           rm=1d0/(1d0+a*(energyref-veffmt (1, ir, ias)*y00))
                           t1=apwfr(ir, 1, io, l1, ias)*lofr(ir, 1, ilo, ias)
                           t2=apwfr(ir, 2, io, l1, ias)*lofr(ir, 2, ilo, ias)
                           fr (ir) = (0.5d0*t2*rm + 0.5d0*angular*t1*rm/spr(ir,is)**2 + t1*veffmt(1, ir, ias)* y00)*r2 (ir)
                        End Do
                        Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                        halointegrals (1, io, l3, ilo) = gr (nr) / y00
! calculate more integrals if linearized Koelling-Harmon is demanded
                        if (input%groundstate%ValenceRelativity.eq.'lkh') then
                          Do ir = 1, nr
                            rm=1d0/(1d0+a*(energyref-veffmt (1, ir, ias)*y00))
                            t1=apwfr(ir, 1, io, l1, ias)*lofr(ir, 1, ilo, ias)
                            t2=apwfr(ir, 2, io, l1, ias)*lofr(ir, 2, ilo, ias)
                            fr (ir) = a*(0.5d0*t2*rm**2 + 0.5d0*angular*t1*rm**2/spr(ir,is)**2)*r2 (ir)
                          End Do
                          Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                          h1loa (ilo, io, l1, ias) = gr (nr) / y00
                        endif
                     Else
                        halointegrals (1, io, l3, ilo) = 0.d0
                     End If
                     Do l2 = 1, input%groundstate%lmaxvr
                        Do m2 = - l2, l2
                           lm2 = idxlm (l2, m2)
                           Do ir = 1, nr
                              t1 = lofr (ir, 1, ilo, ias) * apwfr (ir, &
                             & 1, io, l3, ias) * r2 (ir)
                              fr (ir) = t1 * veffmt (lm2, ir, ias)
                           End Do
                           Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                           halointegrals (lm2, io, l3, ilo) = gr (nr)
                        End Do
                     End Do
                  End Do
               End Do
            End Do
! Angular integrals
     nalo1=haloijSize(is)
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
                     haloij(if1,if3,ias)=zsum
                     if (io.ne.apword(l3,is)) inonz=ireset3
                  End Do
               End Do 
            End Do 
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
                        rm=1d0/(1d0+a*(energyref-veffmt (1, ir, ias)*y00))
                        t1=lofr(ir, 1, ilo1, ias)*lofr(ir, 1, ilo2, ias)
                        t2=lofr(ir, 2, ilo1, ias)*lofr(ir, 2, ilo2, ias)
                        fr (ir) = (0.5d0*t2*rm + 0.5d0*angular*t1*rm/spr(ir,is)**2 + t1*veffmt(1, ir, ias)* y00)*r2 (ir)
                     End Do
                     Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                     hlolo (ilo1, ilo2, 1, ias) = gr (nr) / y00
                     if (input%groundstate%ValenceRelativity.eq.'lkh') then
                       Do ir = 1, nr
                         rm=1d0/(1d0+a*(energyref-veffmt (1, ir, ias)*y00))
                         t1=lofr(ir, 1, ilo1, ias)*lofr(ir, 1, ilo2, ias)
                         t2=lofr(ir, 2, ilo1, ias)*lofr(ir, 2, ilo2, ias)
                         fr (ir) = a*(0.5d0*t2*rm**2 + 0.5d0*angular*t1*rm**2/spr(ir,is)**2)*r2 (ir)
                       End Do
                       Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                       h1lolo (ilo1, ilo2, ias) = gr (nr) / y00
                     endif
                  Else
                     hlolo (ilo1, ilo2, 1, ias) = 0.d0
                  End If
                  Do l2 = 1, input%groundstate%lmaxvr
                     Do m2 = - l2, l2
                        lm2 = idxlm (l2, m2)
                        Do ir = 1, nr
                           t1 = lofr (ir, 1, ilo1, ias) * lofr (ir, 1, &
                          & ilo2, ias) * r2 (ir)
                           fr (ir) = t1 * veffmt (lm2, ir, ias)
                        End Do
                        Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                        hlolo (ilo1, ilo2, lm2, ias) = gr (nr)
                     End Do
                  End Do
               End Do
            End Do

           else
! If kinetic energy is calculated as Laplacian
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
! Radial part
            Do l1 = 0, input%groundstate%lmaxmat
               Do io1 = 1, apword (l1, is)
                  Do l3 = 0, input%groundstate%lmaxapw
                     Do io2 = 1, apword (l3, is)
                        If (l1 .Eq. l3) Then
                           Do ir = 1, nr
                              fr (ir) = apwfr (ir, 1, io1, l1, ias) * apwfr (ir, 2, io2, l3, ias) * r2 (ir) 
                           End Do
                           Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                           hintegrals (1, io2, l3, io1, l1)= gr (nr) / y00
                        Else
                           hintegrals (1, io2, l3, io1, l1) = 0.d0
                        End If
                           Do l2 = 1, input%groundstate%lmaxvr
                              Do m2 = - l2, l2
                                 lm2 = idxlm (l2, m2)
                                 Do ir = 1, nr
                                    t1=apwfr(ir,1,io1,l1,ias)*apwfr(ir,1,io2,l3,ias)*r2(ir)
                                    fr (ir) = t1 * veffmt (lm2, ir, ias)
                                 End Do
                                 Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                                 hintegrals (lm2, io2, l3, io1, l1)=gr (nr)
                              End Do
                           End Do
                     End Do
                  End Do
               End Do
            End Do
! Angular integrals including the surface contribution to the kinetic energy
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
                            zsum=zsum+gntnonz(inonz)*hintegrals (gntnonzlm2(inonz), io2, l3, io1, l1)!  haa(gntnonzlm2(inonz), io2, l3, io1, l1, ias)
                            inonz=inonz+1
                          enddo
                          haaij(if1,if3,ias)=zsum
                          if (lm1.eq.lm3) then
                            haaij(if1,if3,ias)=haaij(if1,if3,ias)+t1*apwfr(nrmt(is),1,io1,l1,ias)*apwdfr(io2,l1,ias)*1d0/(1d0+(energyref-veffmt(1,nrmt(is),ias)*y00)*a)
                          endif
                          if (io2.ne.apword (l3, is)) inonz=ireset3
                        End Do
                  End Do
               End Do
              if (io1.ne.apword (l1, is)) inonz=ireset1
            End Do
         End Do
      End Do

!--------------------------------------!
!     local-orbital-APW integtrals     !
!--------------------------------------!
! Radial integrals
            Do ilo = 1, nlorb (is)
               l1 = lorbl (ilo, is)
               Do l3 = 0, input%groundstate%lmaxmat
                  Do io = 1, apword (l3, is)
                     If (l1 .Eq. l3) Then
                        Do ir = 1, nr
                           fr (ir) = lofr (ir, 1, ilo, ias) * apwfr(ir, 2, io, l3, ias) * r2 (ir)
                        End Do
                        Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                        halointegrals (1, io, l3, ilo) = gr (nr) / y00
                     Else
                        halointegrals (1, io, l3, ilo) = 0.d0
                     End If
                     Do l2 = 1, input%groundstate%lmaxvr
                        Do m2 = - l2, l2
                           lm2 = idxlm (l2, m2)
                           Do ir = 1, nr
                              t1 = lofr (ir, 1, ilo, ias) * apwfr (ir, 1, io, l3, ias) * r2 (ir)
                              fr (ir) = t1 * veffmt (lm2, ir, ias)
                           End Do
                           Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                           halointegrals (lm2, io, l3, ilo) = gr (nr)
                        End Do
                     End Do 
                  End Do 
               End Do 
            End Do
! Angular integrals 
     nalo1=haloijSize(is)
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
                     haloij(if1,if3,ias)=zsum
                     if (io.ne.apword(l3,is)) inonz=ireset3
                  End Do
               End Do
            End Do
         End Do
      End Do

!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
! Radial integrals
            Do ilo1 = 1, nlorb (is)
               l1 = lorbl (ilo1, is)
               Do ilo2 = 1, nlorb (is)
                  l3 = lorbl (ilo2, is)
                  If (l1 .Eq. l3) Then
                     Do ir = 1, nr
                        fr (ir) = lofr (ir, 1, ilo1, ias) * lofr (ir, &
                       & 2, ilo2, ias) * r2 (ir)
                     End Do
                     Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                     hlolo (ilo1, ilo2, 1, ias) = gr (nr) / y00
                  Else
                     hlolo (ilo1, ilo2, 1, ias) = 0.d0
                  End If
                  Do l2 = 1, input%groundstate%lmaxvr
                     Do m2 = - l2, l2
                        lm2 = idxlm (l2, m2)
                        Do ir = 1, nr
                           t1 = lofr (ir, 1, ilo1, ias) * lofr (ir, 1, &
                          & ilo2, ias) * r2 (ir)
                           fr (ir) = t1 * veffmt (lm2, ir, ias)
                        End Do
                        Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                        hlolo (ilo1, ilo2, lm2, ias) = gr (nr)
                     End Do
                  End Do
               End Do
            End Do
           endif
! end loops over atoms and species
         End Do
      End Do
      deallocate(hintegrals)
      Return
End Subroutine
!EOC
