!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: olprad
! !INTERFACE:
!
!
Subroutine olprad
! !USES:
      Use modmain
! !DESCRIPTION:
!   Calculates the radial overlap integrals of the APW and local-orbital basis
!   functions. In other words, for spin $\sigma$ and atom $j$ of species $i$, it
!   computes integrals of the form
!   $$ o^{\sigma;ij}_{qp}=\int_0^{R_i}u^{\sigma;ij}_{q;l_p}(r)v^{\sigma;ij}_p(r)
!    r^2dr $$
!   and
!   $$ o^{\sigma;ij}_{pp'}=\int_0^{R_i}v^{\sigma;ij}_p(r)v^{\sigma;ij}_{p'}(r)
!    r^2dr,\quad l_p=l_{p'} $$
!   where $u^{\sigma;ij}_{q;l}$ is the $q$th APW radial function for angular
!   momentum $l$; and $v^{\sigma;ij}_p$ is the $p$th local-orbital radial
!   function and has angular momentum $l_p$.
!
! !REVISION HISTORY:
!   Created November 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, ir, nr
      Integer :: l, ilo, ilo1, ilo2, io,io1,io2
! automatic arrays
      Real (8) :: r2 (nrmtmax), fr (nrmtmax), gr (nrmtmax), cf (3, &
     & nrmtmax)
      Real (8) :: angular,t1,t2,rm,a,alpha
      parameter (alpha=1d0 / 137.03599911d0)

      h1aa=0d0
      h1loa=0d0
      h1lolo=0d0

      if (input%groundstate%ValenceRelativity.ne.'none') then
        a=0.5d0*alpha**2
      else
        a=0d0
      endif


      Do is = 1, nspecies
         nr = nrmt (is)
         Do ir = 1, nr
            r2 (ir) = spr (ir, is) ** 2
         End Do
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)

            if (input%groundstate%ValenceRelativity.eq.'iora*') then
              Do l = 0, input%groundstate%lmaxmat
                angular=dble(l*(l+1))
                Do io1 = 1, apword (l, is)
                  Do io2 = 1, apword (l, is)
! calculate more integrals if linearized Koelling-Harmon is demanded
                    Do ir = 1, nr
                      rm=1d0/(1d0-a*veffmt (1, ir, ias)*y00)
                      t1=apwfr(ir, 1, io1, l, ias)*apwfr(ir, 1, io2, l, ias)
                      t2=apwfr(ir, 2, io1, l, ias)*apwfr(ir, 2, io2, l, ias)
                      fr (ir) = a*(0.5d0*t2*rm**2 + 0.5d0*angular*t1*rm**2/spr(ir,is)**2)*r2 (ir)
                    End Do
                    Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                    h1aa (io1, io2, l, ias) = gr(nr)
!                    if (io1.eq.io2) then
!                      h1aa (io1, io2, l, ias) = 1d0+gr(nr) 
!                    else
!                      h1aa (io1, io2, l, ias) = gr(nr)
!                    endif
                  End Do
!                 h1aa(io1,io1,l,ias)=1d0+h1aa(io1,io1,l,ias)
                End Do
              End Do
            endif

!--------------------------------------!
!     APW-local-orbital integtrals     !
!--------------------------------------!
            Do ilo = 1, nlorb (is)
               l = lorbl (ilo, is)
               Do io = 1, apword (l, is)
!                 if (input%groundstate%ValenceRelativity.ne.'lkh') then
                   Do ir = 1, nr
                     fr (ir) = apwfr (ir, 1, io, l, ias) * lofr (ir, 1, ilo, ias) * r2 (ir)
                   End Do
                  Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                  oalo (io, ilo, ias) = gr (nr) 
!                 else
                  if (input%groundstate%ValenceRelativity.eq.'iora*') then
                   angular=dble(l*(l+1))
                   Do ir = 1, nr
                     rm=1d0/(1d0-a*veffmt (1, ir, ias)*y00)
                     t1=apwfr(ir, 1, io, l, ias)*lofr(ir, 1, ilo, ias)
                     t2=apwfr(ir, 2, io, l, ias)*lofr(ir, 2, ilo, ias)
                     fr (ir) = (a*(0.5d0*t2*rm**2 + 0.5d0*angular*t1*rm**2/spr(ir,is)**2))*r2 (ir)
                   End Do
                   Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                   h1loa (io, ilo, ias) = gr (nr)
                 endif
               End Do
            End Do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
            Do ilo1 = 1, nlorb (is)
               l = lorbl (ilo1, is)
               Do ilo2 = 1, nlorb (is)
                  If (lorbl(ilo2, is) .Eq. l) Then
!                    if (input%groundstate%ValenceRelativity.ne.'lkh') then
                      Do ir = 1, nr
                        fr (ir) = lofr (ir, 1, ilo1, ias) * lofr (ir, 1, ilo2, ias) * r2 (ir)
                      End Do
                      Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                      ololo (ilo1, ilo2, ias) = gr (nr)
                    if (input%groundstate%ValenceRelativity.eq.'iora*') then
                      angular=dble(l*(l+1))
                      Do ir = 1, nr
                        rm=1d0/(1d0-a*veffmt (1, ir, ias)*y00)
                        t1=lofr(ir, 1, ilo1, ias)*lofr(ir, 1, ilo2, ias)
                        t2=lofr(ir, 2, ilo1, ias)*lofr(ir, 2, ilo2, ias)
                        fr (ir) = (a*(0.5d0*t2*rm**2 + 0.5d0*angular*t1*rm**2/spr(ir,is)**2))*r2 (ir)
                      End Do
                      Call fderiv (-1, nr, spr(:, is), fr, gr, cf)
                      h1lolo (ilo1, ilo2, ias) = gr (nr)
                    endif
                  End If
               End Do
            End Do
         End Do
      End Do
      Return
End Subroutine
!EOC
