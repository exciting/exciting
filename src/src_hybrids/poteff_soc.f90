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
Subroutine poteff_soc(veffmt_pbe)
! !USES:
      Use modmain
      Use mod_muffin_tin
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
      real(8), intent(out):: veffmt_pbe(lmmaxvr,nrmtmax,natmtot)      
      integer::xctype_
      logical :: flag

!      if (allocated(veffmt_pbe)) deallocate(veffmt_pbe)
!      write(*,*) lmmaxvr, nrmtmax, natmtot
!      allocate(veffmt_pbe(lmmaxvr,nrmtmax,natmtot))     
!---------------------------------------------
! compute the exchange-correlation potential
!---------------------------------------------
      ex_coef = 0.d0
      ec_coef = 1.d0
      
      Call potxc
      call potcoul
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
                  veffmt_pbe(lm,ir,ias) = vclmt(lm,ir,ias) + vxcmt(lm,ir,ias)
               End Do
               Do lm = lmmax + 1, lmmaxvr
                  veffmt_pbe(lm,ir,ias) = 0.d0
               End Do
            End Do
         End Do
      End Do
      ex_coef = input%groundstate%Hybrid%excoeff  
      ec_coef = input%groundstate%Hybrid%eccoeff  
      Return
End Subroutine

