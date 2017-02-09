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
      real (8) :: shift
      Real (8) :: ts0, ts1, ta, tb
      
      Call timesec (ts0)
      
!---------------------------------------------
! compute the exchange-correlation potential
!---------------------------------------------
      Call potxc

!---------------------------------
! compute the Coulomb potential
!---------------------------------
      Call potcoul
      shift=input%groundstate%energyref

!----------------------------------------------------------
! add Coulomb and exchange-correlation potentials together
!----------------------------------------------------------
      
      ! muffin-tin part
      vclmt(1,:,:) = vclmt(1,:,:)+shift/y00
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            lmmax = lmmaxinr
            Do ir = 1, nrmt (is)
               If (ir .Gt. nrmtinr(is)) lmmax = lmmaxvr
               Do lm = 1, lmmax
                  if (associated(input%groundstate%dfthalf)) then
                    veffmt(lm,ir,ias) = vclmt(lm,ir,ias) + vxcmt(lm,ir,ias) + vhalfmt (lm, ir, ias)
                  else
                    veffmt(lm,ir,ias) = vclmt(lm,ir,ias) + vxcmt(lm,ir,ias)
                  endif
               End Do
               Do lm = lmmax + 1, lmmaxvr
                  veffmt(lm,ir,ias) = 0.d0
               End Do
            End Do
         End Do
      End Do
      
      ! interstitial part
      vclir(:) = vclir(:) + shift
      
      if (associated(input%groundstate%dfthalf)) then
        veffir(:) = vclir(:) + vxcir(:) + vhalfir(:)
      else
        veffir(:) = vclir(:) + vxcir(:)
      endif
      
      Call timesec (ts1)
      timepot = timepot + ts1 - ts0
      
      Return
End Subroutine
!EOC
