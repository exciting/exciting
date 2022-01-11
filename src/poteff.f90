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
      Use vx_enums, only: HYB_PBE0, HYB_HSE
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
      call timesec(ta)
!---------------------------------
! compute the Coulomb potential
!---------------------------------
      Call potcoul
      call timesec(tb)
      shift=input%groundstate%energyref

!        write(*,*) 'potxc', ta-ts0
!        write(*,*) 'potcoul', tb-ta
!----------------------------------------------------------
! add Coulomb and exchange-correlation potentials together
!----------------------------------------------------------
      
      ! muffin-tin part
      vclmt(1,:,:) = vclmt(1,:,:)+shift/y00
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            lmmax = lmmaxinr
! The loops should be reordered here. Cheers, Andris.
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

! Should have been more general, something like "if (hybrid) then"
      If  (xctype(1)==HYB_PBE0 .or. xctype(1)==HYB_HSE) Then
! Is it fine to exclude vhalfir(:) from the effective potential used in the relativistic kinetic energy operator?
! How likely it is that we try to run DFT-1/2 combined with hybrids?
        vrelmt = vrelmt + vclmt 
      endif

      ! interstitial part
      vclir(:) = vclir(:) + shift
      
      if (associated(input%groundstate%dfthalf)) then
        veffir(:) = vclir(:) + vxcir(:) + vhalfir(:)
      else
        veffir(:) = vclir(:) + vxcir(:)
      endif

! Should have been more general, something like "if (hybrid) then"
      If  (xctype(1)==HYB_PBE0 .or. xctype(1)==HYB_HSE) Then
! Is it fine to exclude vhalfir(:) from the effective potential used in the relativistic kinetic energy operator?
        vrelir = vrelir + vclir
      endif


      Call timesec (ts1)
      timepot = timepot + ts1 - ts0
      
      Return
End Subroutine
!EOC
