!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: poteff
! !INTERFACE:
!
!
Subroutine poteff
! !USES:
      Use modmain
! !DESCRIPTION:
!   Computes the effective potential by adding together the Coulomb and
!   exchange-correlation potentials. See routines {\tt potcoul} and {\tt potxc}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, ir, lm, lmmax
      Real (8) :: ts0, ts1
      Call timesec (ts0)
! compute the Coulomb potential
      Call potcoul
! compute the exchange-correlation potential
      Call potxc
! add Coulomb and exchange-correlation potentials together
! muffin-tin part
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
!            write(*,*) vclmt(:,nrmt(is),ias)
!            write(*,*)
            lmmax = lmmaxinr
            Do ir = 1, nrmt (is)
               If (ir .Gt. nrmtinr(is)) lmmax = lmmaxvr
               Do lm = 1, lmmax
                  veffmt (lm, ir, ias) = vclmt (lm, ir, ias) + vxcmt &
                 & (lm, ir, ias)
               End Do
               Do lm = lmmax + 1, lmmaxvr
                  veffmt (lm, ir, ias) = 0.d0
               End Do
            End Do
!            write(*,*) veffmt (1, nrmt(is), ias)
         End Do
      End Do
      
!      veffmt (1, :, :)=veffmt (1, :, :)+input%groundstate%energyref*y00
!      write(*,*) 
! interstitial part
      veffir (:) = vclir (:) + vxcir (:)!+input%groundstate%energyref
!      do ir=1,nrmt(1)
!        write(*,*) spr(ir,1),veffmt (1, ir, 1)
!      enddo
      
!      do ir=1,60
!        write(*,*)  4.185743066d0*sqrt(2d0)*dble(ir-1),veffir (ir)
!      enddo
!      stop
      Call timesec (ts1)
!      timepot = timepot + ts1 - ts0
      Return
End Subroutine
!EOC
