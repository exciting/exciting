!
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: rhoplot
! !INTERFACE:
!
!
Subroutine rhoplot
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Outputs the charge density and the charge density gradients (modulus)
!   read in from {\tt STATE.OUT}, for 1D, 2D or 3D
!   plotting.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!   Density gradients are added, March 2011 (DIN)
!EOP
!BOC
      Implicit None
! Local variables for gradient calculation
      Real*8, Allocatable :: grhomt(:,:,:,:), grhoir(:,:)
      Real*8, Allocatable :: modgrhomt(:,:,:), modgrhoir(:)
      Integer :: i, j, k
! initialise universal variables
      Call init0
! read density from file
      Call readstate
! when gradients are needed
      If (input%properties%chargedensityplot%plotgradient) Then
! allocate density gradient arrays
          If (allocated(grhomt)) deallocate (grhomt)
          If (allocated(modgrhomt)) deallocate (modgrhomt)
          Allocate (grhomt(lmmaxvr, nrmtmax, natmtot, 3), &
                  & modgrhomt(lmmaxvr, nrmtmax, natmtot))
          If (allocated(grhoir)) deallocate (grhoir)
          If (allocated(modgrhoir)) deallocate (modgrhoir)
          Allocate (grhoir(ngrtot,3), modgrhoir(ngrtot))
! calculate gradients
          Call gradrf (rhomt, rhoir, grhomt, grhoir)
! get gradient modules
! muffin-tin area
          Do i = 1, lmmaxvr
             Do j = 1, nrmtmax
                Do k = 1, natmtot
                   modgrhomt (i, j, k) = dsqrt ( grhomt(i, j, k, 1)**2 + &
                                               & grhomt(i, j, k, 2)**2 + &
                                               & grhomt(i, j, k, 3)**2 )
                End Do 
             End Do
          End Do  
! interstitial
          Do i = 1, ngrtot
             modgrhoir (i) = dsqrt ( grhoir(i, 1)**2 + &
           & grhoir(i, 2)**2 + grhoir(i, 3)**2 )
          End Do
      End If
! write the density plot to file
      If (associated(input%properties%chargedensityplot%plot1d)) Then
!
         Call plot1d ("RHO", 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rhomt, rhoir, input%properties%chargedensityplot%plot1d)
!
         Write (*,*)
         Write (*, '("Info(rhoplot):")')
         Write (*, '(" 1D density plot written to RHO1D.OUT")')
         Write (*, '(" vertex location lines written to RHOLINES.OUT")')
!
! when gradients are needed      
         If (input%properties%chargedensityplot%plotgradient) Then
             Call plot1d ("GRHO", 1, input%groundstate%lmaxvr, lmmaxvr, &
             & modgrhomt, modgrhoir, input%properties%chargedensityplot%plot1d)
!
             Write (*,*)
             Write (*, '("Info(rhoplot): 1D module of density gradient plot written to GRHO1D.OUT")')
         End If
!      
      End If
      If (associated(input%properties%chargedensityplot%plot2d)) Then
!
         Call plot2d ("RHO", 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rhomt, rhoir, input%properties%chargedensityplot%plot2d)
!
         Write (*,*)
         Write (*, '("Info(rhoplot): 2D density plot written to RHO2D.OUT")')
! when gradients are needed      
         If (input%properties%chargedensityplot%plotgradient) Then
             Call plot2d ("GRHO", 1, input%groundstate%lmaxvr, lmmaxvr, &
             & modgrhomt, modgrhoir, input%properties%chargedensityplot%plot2d)
!
             Write (*,*)
             Write (*, '("Info(rhoplot): 2D module of density gradient plot written to GRHO2D.OUT")')
         End If
!      
      End If
      If (associated(input%properties%chargedensityplot%plot3d)) Then
         Call plot3d ("RHO", 1, input%groundstate%lmaxvr, lmmaxvr, &
        & rhomt, rhoir, input%properties%chargedensityplot%plot3d)
         Write (*,*)
         Write (*, '("Info(rhoplot): 3D density plot written to RHO3D.OUT")')
! when gradients are needed      
         If (input%properties%chargedensityplot%plotgradient) Then
             Call plot3d ("GRHO", 1, input%groundstate%lmaxvr, lmmaxvr, &
             & modgrhomt, modgrhoir, input%properties%chargedensityplot%plot3d)
!
             Write (*,*)
             Write (*, '("Info(rhoplot): 3D module of density gradient plot written to GRHO3D.OUT")')
         End If
!      
      End If
      Write (*,*)
      Deallocate (grhomt, grhoir)
      Return
End Subroutine
!EOC
