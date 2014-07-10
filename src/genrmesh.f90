!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genrmesh
! !INTERFACE:
!
!
Subroutine genrmesh
! !USES:
      Use modinput
      Use mod_atoms
      Use mod_muffin_tin
! !DESCRIPTION:
!   Generates the coarse and fine radial meshes for each atomic species in the
!   crystal. Also determines which points are in the inner part of the
!   muffin-tin using the value of {\tt fracinr}.
!   For species $i$ the radial mesh starts from the value $R_0=r(1)$, hits the
!   muffin-tin surface at $R_{\rm MT}=r(N)$, and ends at the effective infinity
!   value $R_{\infty}=r(N_{\infty})$. The number of points up to the effective
!   infinity are determined by the number of points $N$ up to the muffin-tin
!   radius as well as by the smallest and largest mesh point and the muffin-tin
!   radius, and is given by
!   $$ N_{\infty}={\rm round}
!   \left[ \frac{(N-1)\ln(R_{\infty}/R_0)}{\ln(R_{\rm MT}/R_0)} \right] +1.
!   $$
!   The radial mesh points are finally defined by
!   $$ r(j)= R_0\left(\frac{R_{\rm MT}}{R_0}\right)^\frac{j-1}{N-1},
!   $$
!   for $j=1,\ldots,N_\infty$. \\
!   Note: The number of mesh points initially defined in species file is adapted
!   to be commensurate with the coarse mesh of step size {\tt lradstep}
!   $$ N = N^*-\mod(N^*,{\tt lradstep}),
!   $$
!   if $N^*$ was the number of points defined in the species file.
!   The number of mesh points $\tilde{N} $ of the coarse mesh $\tilde{r}(j)$
!   reads
!   $$ \tilde{N} = \left\lfloor \frac{N-1}{\tt lradstep}\right\rfloor + 1.
!   $$
!   It is given by
!   $$ \tilde{r}(j) = r([j-1]*{\tt lradstep}+1),
!   $$
!   for $j=1,\ldots,\tilde{N}$ and has the properties $\tilde{r}(1)=r(1)=R_0$
!   and $\tilde{r}(\tilde{N})=r(N)=R_{\rm MT}$.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!   Revised and updated documentation, April 2011 (S. Sagmeister)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ir, irc,cutoff
      Real (8) :: t1, t2
      Real (8) :: erf
      external erf
! estimate the number of radial mesh points to infinity
      if ((input%groundstate%radialgridtype.ne."cubic").and. &
          (input%groundstate%radialgridtype.ne."expocubic").and. &
          (input%groundstate%radialgridtype.ne."exponential")) then 
         write(*,*) 'Wrong radialGridType.'
         write(*,*) 'Choose between cubic, expocubic and exponential!'
         write(*,*) 'Terminating...'
         stop
      endif
      spnrmax = 1
      Do is = 1, nspecies
         if (input%groundstate%radialgridtype.eq."exponential") then
           t1 = Log (sprmax(is)/sprmin(is)) / Log (rmt(is)/sprmin(is))
           t2 = dble (nrmt(is)-1) * t1
           spnr (is) = Nint (t2) + 1
         else
           spnr (is) =2+Nint(dble(nrmt(is)-1)*((sprmax(is)-sprmin(is))/(rmt(is)-sprmin(is)))**0.333333333333333d0)
         endif
           spnrmax = Max (spnrmax, spnr(is))
      End Do
! allocate the global radial mesh arrays
      If (allocated(spr)) deallocate (spr)
      Allocate (spr(spnrmax, nspecies))
      If (allocated(rcmt)) deallocate (rcmt)
      Allocate (rcmt(nrcmtmax, nspecies))
! generate the radial meshes
      Do is = 1, nspecies
         t1 = 1.d0 / dble (nrmt(is)-1)
! logarithmic mesh
         t2 = Log (rmt(is)/sprmin(is))
!         cutoff=min(Nint(nrmt(is)*0.5d0),150)
         cutoff=Nint(nrmt(is)*0.5d0)

         Do ir = 1, spnr (is)
           if (input%groundstate%radialgridtype.eq."cubic") then
            spr (ir, is) = sprmin (is)+(dble(ir-1)/dble(nrmt(is)-1))**3*(rmt(is)-sprmin (is))
!            spr (ir, is) = sprmin (is)*dble(ir)+(dble(ir-1)/dble(nrmt(is)-1))**3*(rmt(is)-sprmin (is)*dble(nrmt(is)))
!             spr (ir, is) = sprmin (is)*dble(2*ir-1)+(dble(ir-1)/dble(nrmt(is)-1))**3*(rmt(is)-sprmin (is)*dble(2*nrmt(is)-1))
           elseif (input%groundstate%radialgridtype.eq."exponential") then
             spr (ir, is) = sprmin (is) * Exp (dble(ir-1)*t1*t2)
           else
             spr (ir, is) = 0.5d0*(erf(5d0*dble(ir-cutoff)/nrmt(is))+1d0)* &
            &  (sprmin (is)+(dble(ir-1)/dble(nrmt(is)-1))**3*(rmt(is)-sprmin (is)))+ &
            &  (1d0-0.5d0*(erf(5d0*dble(ir-cutoff)/nrmt(is))+1d0))*sprmin (is) * Exp(dble(ir-1)*t1*t2)
           endif
!            write(*,*) spr (ir, is)
         End Do
!         stop
      End Do
! find the inner part of the muffin-tin (where rho is calculated with lmaxinr)
      Do is = 1, nspecies
         t1 = input%groundstate%fracinr * rmt (is)
         nrmtinr (is) = 0
         Do ir = 1, nrmt (is)
            If (spr(ir, is) .Gt. t1) Then
               nrmtinr (is) = ir
               Go To 10
            End If
         End Do
10       Continue
      End Do
! set up the coarse radial meshes
      Do is = 1, nspecies
         irc = 0
         Do ir = 1, nrmt (is), input%groundstate%lradstep
            irc = irc + 1
            rcmt (irc, is) = spr (ir, is)
         End Do
      End Do
      Return
End Subroutine
!EOC
