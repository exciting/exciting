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
      Integer :: is, ia, ias, nr, ir
      Integer :: l1, l2, l3, m2, lm2
      Integer :: ilo, ilo1, ilo2, io, io1, io2
      Real (8) :: t1,t2,angular
! automatic arrays
      Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax),a,rm,alpha
      parameter (alpha=1d0 / 137.03599911d0)

      If (allocated(haa)) deallocate (haa)
      Allocate (haa(apwordmax, 0:input%groundstate%lmaxmat, apwordmax, 0:input%groundstate%lmaxapw, lmmaxvr, natmtot))
      If (allocated(hloa)) deallocate (hloa) 
      Allocate (hloa(nlomax, apwordmax, 0:input%groundstate%lmaxmat, lmmaxvr, natmtot))
      If (allocated(hlolo)) deallocate (hlolo)
      Allocate (hlolo(nlomax, nlomax, lmmaxvr, natmtot))


! begin loops over atoms and species
      Do is = 1, nspecies
         nr = nrmt (is)
         Do ir = 1, nr
            r2 (ir) = spr (ir, is) ** 2
         End Do
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
           if (input%groundstate%SymmetricKineticEnergy) then
            if (input%groundstate%ValenceRelativity.ne.'none') then
              a=0.5d0*alpha**2
            else
              a=0d0
            endif
            Do l1 = 0, input%groundstate%lmaxmat
               Do io1 = 1, apword (l1, is)
                  Do l3 = 0, input%groundstate%lmaxapw
                     Do io2 = 1, apword (l3, is)
                        If (l1 .Eq. l3) Then
                           angular=dble(l1*(l1+1))
                           Do ir = 1, nr
                              rm=1d0/(1d0-a*veffmt (1, ir, ias)*y00)
                              t1=apwfr(ir, 1, io1, l1, ias)*apwfr(ir, 1, io2, l3, ias)
                              t2=apwfr(ir, 2, io1, l1, ias)*apwfr(ir, 2, io2, l3, ias)
                              fr (ir) = (0.5d0*t2*rm + 0.5d0*angular*t1*rm/spr(ir,is)**2 + t1*veffmt(1, ir, ias)* y00)*r2 (ir)
                           End Do
                           Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                           haa (io1, l1, io2, l3, 1, ias) = gr (nr) / y00
! calculate more integrals if linearized Koelling-Harmon is demanded
                           if (input%groundstate%ValenceRelativity.eq.'iora*') then
                             Do ir = 1, nr
                               rm=1d0/(1d0-a*veffmt (1, ir, ias)*y00)
                               t1=apwfr(ir, 1, io1, l1, ias)*apwfr(ir, 1, io2, l3, ias)
                               t2=apwfr(ir, 2, io1, l1, ias)*apwfr(ir, 2, io2, l3, ias)
                               fr (ir) = a*(0.5d0*t2*rm**2 + 0.5d0*angular*t1*rm**2/spr(ir,is)**2)*r2 (ir)
                             End Do
                             Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                             h1aa (io1, io2, l1, ias) = gr (nr) / y00
                           endif
                        Else
                           haa (io1, l1, io2, l3, 1, ias) = 0.d0
                        End If 
                        If (l1 .Ge. l3) Then
                           Do l2 = 1, input%groundstate%lmaxvr
                              Do m2 = - l2, l2
                                 lm2 = idxlm (l2, m2)
                                 Do ir = 1, nr
                                    t1=apwfr(ir,1,io1,l1,ias)*apwfr(ir,1,io2,l3,ias)*r2(ir)
                                    fr (ir) = t1 * veffmt (lm2, ir, ias)
                                 End Do
                                 Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                                 haa (io1, l1, io2, l3, lm2, ias) = gr (nr)
                              End Do
                           End Do
                        End If
                     End Do
                  End Do
               End Do
            End Do
!--------------------------------------!
!     local-orbital-APW integtrals     !
!--------------------------------------!
            Do ilo = 1, nlorb (is)
               l1 = lorbl (ilo, is)
               Do l3 = 0, input%groundstate%lmaxmat
                  Do io = 1, apword (l3, is)
                     If (l1 .Eq. l3) Then
                        angular=dble(l1*(l1+1))
                        Do ir = 1, nr
                           rm=1d0/(1d0-a*veffmt (1, ir, ias)*y00)
                           t1=apwfr(ir, 1, io, l1, ias)*lofr(ir, 1, ilo, ias)
                           t2=apwfr(ir, 2, io, l1, ias)*lofr(ir, 2, ilo, ias)
                           fr (ir) = (0.5d0*t2*rm + 0.5d0*angular*t1*rm/spr(ir,is)**2 + t1*veffmt(1, ir, ias)* y00)*r2 (ir)
                        End Do
                        Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                        hloa (ilo, io, l3, 1, ias) = gr (nr) / y00
! calculate more integrals if linearized Koelling-Harmon is demanded
                        if (input%groundstate%ValenceRelativity.eq.'iora*') then
                          Do ir = 1, nr
                            rm=1d0/(1d0-a*veffmt (1, ir, ias)*y00)
                            t1=apwfr(ir, 1, io, l1, ias)*lofr(ir, 1, ilo, ias)
                            t2=apwfr(ir, 2, io, l1, ias)*lofr(ir, 2, ilo, ias)
                            fr (ir) = a*(0.5d0*t2*rm**2 + 0.5d0*angular*t1*rm**2/spr(ir,is)**2)*r2 (ir)
                          End Do
                          Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                          h1loa (io, ilo, ias) = gr (nr) / y00
                        endif
                     Else
                        hloa (ilo, io, l3, 1, ias) = 0.d0
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
                           hloa (ilo, io, l3, lm2, ias) = gr (nr)
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
                        rm=1d0/(1d0-a*veffmt (1, ir, ias)*y00)
                        t1=lofr(ir, 1, ilo1, ias)*lofr(ir, 1, ilo2, ias)
                        t2=lofr(ir, 2, ilo1, ias)*lofr(ir, 2, ilo2, ias)
                        fr (ir) = (0.5d0*t2*rm + 0.5d0*angular*t1*rm/spr(ir,is)**2 + t1*veffmt(1, ir, ias)* y00)*r2 (ir)
                     End Do
                     Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                     hlolo (ilo1, ilo2, 1, ias) = gr (nr) / y00
                     if (input%groundstate%ValenceRelativity.eq.'iora*') then
                       Do ir = 1, nr
                         rm=1d0/(1d0-a*veffmt (1, ir, ias)*y00)
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

!           write(*,*) haa (1, 0, 1, 0, 1, 1)
           else
            Do l1 = 0, input%groundstate%lmaxmat
               Do io1 = 1, apword (l1, is)
                  Do l3 = 0, input%groundstate%lmaxapw
                     Do io2 = 1, apword (l3, is)
                        If (l1 .Eq. l3) Then
                           Do ir = 1, nr
                              fr (ir) = apwfr (ir, 1, io1, l1, ias) * apwfr (ir, 2, io2, l3, ias) * r2 (ir)
                           End Do
                           Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                           haa (io1, l1, io2, l3, 1, ias) = gr (nr) / y00
                        Else
                           haa (io1, l1, io2, l3, 1, ias) = 0.d0
                        End If
                        If (l1 .Ge. l3) Then
                           Do l2 = 1, input%groundstate%lmaxvr
                              Do m2 = - l2, l2
                                 lm2 = idxlm (l2, m2)
                                 Do ir = 1, nr
                                    t1=apwfr(ir,1,io1,l1,ias)*apwfr(ir,1,io2,l3,ias)*r2(ir)
                                    fr (ir) = t1 * veffmt (lm2, ir, ias)
                                 End Do
                                 Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                                 haa (io1, l1, io2, l3, lm2, ias) = gr (nr)
                              End Do
                           End Do
                        End If
                     End Do
                  End Do
               End Do
            End Do
!            write(*,*) haa (1, 0, 1, 0, 1, 1), 0.25d0 * rmt (is) ** 2 * apwfr (nrmt(1), 1, 1, 0, 1) * apwdfr (1, 0, 1)
!--------------------------------------!
!     local-orbital-APW integtrals     !
!--------------------------------------!
            Do ilo = 1, nlorb (is)
               l1 = lorbl (ilo, is)
               Do l3 = 0, input%groundstate%lmaxmat
                  Do io = 1, apword (l3, is)
                     If (l1 .Eq. l3) Then
                        Do ir = 1, nr
                           fr (ir) = lofr (ir, 1, ilo, ias) * apwfr &
                          & (ir, 2, io, l3, ias) * r2 (ir)
                        End Do
                        Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                        hloa (ilo, io, l3, 1, ias) = gr (nr) / y00
                     Else
                        hloa (ilo, io, l3, 1, ias) = 0.d0
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
                           hloa (ilo, io, l3, lm2, ias) = gr (nr)
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
      Return
End Subroutine
!EOC
