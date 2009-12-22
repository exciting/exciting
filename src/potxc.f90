!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: potxc
! !INTERFACE:
!
!
Subroutine potxc
! !USES:
      Use modinput
      Use modmain
      Use modxcifc
! !DESCRIPTION:
!   Computes the exchange-correlation potential and energy density. In the
!   muffin-tin, the density is transformed from spherical harmonic coefficients
!   $\rho_{lm}$ to spherical coordinates $(\theta,\phi)$ with a backward
!   spherical harmonic transformation (SHT). Once calculated, the
!   exchange-correlation potential and energy density are transformed with a
!   forward SHT.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: m, is, ia, ias, ir, itp, i
      Real (8) :: t1, t2, t3, t4
! allocatable arrays
      Real (8), Allocatable :: rftp1 (:, :)
      Real (8), Allocatable :: rftp2 (:, :)
      Real (8), Allocatable :: rfir (:, :)
      Real (8), Allocatable :: ex (:)
      Real (8), Allocatable :: ec (:)
      Real (8), Allocatable :: vx (:, :)
      Real (8), Allocatable :: vc (:, :)
      Real (8), Allocatable :: grhomt (:, :)
      Real (8), Allocatable :: gupmt (:, :)
      Real (8), Allocatable :: gdnmt (:, :)
      Real (8), Allocatable :: g2upmt (:, :)
      Real (8), Allocatable :: g2dnmt (:, :)
      Real (8), Allocatable :: g3rhomt (:, :)
      Real (8), Allocatable :: g3upmt (:, :)
      Real (8), Allocatable :: g3dnmt (:, :)
      Real (8), Allocatable :: grhoir (:)
      Real (8), Allocatable :: gupir (:)
      Real (8), Allocatable :: gdnir (:)
      Real (8), Allocatable :: g2upir (:)
      Real (8), Allocatable :: g2dnir (:)
      Real (8), Allocatable :: g3rhoir (:)
      Real (8), Allocatable :: g3upir (:)
      Real (8), Allocatable :: g3dnir (:)
      Allocate (ex(lmmaxvr))
      Allocate (ec(lmmaxvr))
      Allocate (rftp1(lmmaxvr, 4))
      Allocate (rftp2(lmmaxvr, 2))
      m = Max (lmmaxvr, ngrtot)
      If (associated(input%groundstate%spin)) Then
         Allocate (rfir(ngrtot, 2))
         Allocate (vx(m, 2), vc(m, 2))
      Else
         Allocate (vx(m, 1), vc(m, 1))
      End If
      If (xcgrad .Eq. 1) Then
         Allocate (grhomt(lmmaxvr, nrmtmax))
         Allocate (gupmt(lmmaxvr, nrmtmax))
         Allocate (gdnmt(lmmaxvr, nrmtmax))
         Allocate (g2upmt(lmmaxvr, nrmtmax))
         Allocate (g2dnmt(lmmaxvr, nrmtmax))
         Allocate (g3rhomt(lmmaxvr, nrmtmax))
         Allocate (g3upmt(lmmaxvr, nrmtmax))
         Allocate (g3dnmt(lmmaxvr, nrmtmax))
         Allocate (grhoir(ngrtot))
         Allocate (gupir(ngrtot))
         Allocate (gdnir(ngrtot))
         Allocate (g2upir(ngrtot))
         Allocate (g2dnir(ngrtot))
         Allocate (g3rhoir(ngrtot))
         Allocate (g3upir(ngrtot))
         Allocate (g3dnir(ngrtot))
      End If
!---------------------------------------!
!     muffin-tin potential and field    !
!---------------------------------------!
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! gradients for GGA if required
            If (xcgrad .Eq. 1) Call ggamt (is, ia, grhomt, gupmt, &
           & gdnmt, g2upmt, g2dnmt, g3rhomt, g3upmt, g3dnmt)
            Do ir = 1, nrmt (is)
               If (associated(input%groundstate%spin)) Then
!------------------------!
!     spin-polarised     !
!------------------------!
                  If (ncmag) Then
! non-collinear
                     Do i = 1, 3
                        Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, &
                       & rbshtvr, lmmaxvr, magmt(:, ir, ias, i), 1, &
                       & 0.d0, rftp1(:, i), 1)
                     End Do
                     Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                    & lmmaxvr, rhomt(:, ir, ias), 1, 0.d0, rftp1(:, 4), &
                    & 1)
! compute (rho+|m|)/2 and (rho-|m|)/2
                     Do itp = 1, lmmaxvr
                        t1 = rftp1 (itp, 4)
                        t2 = Sqrt (rftp1(itp, 1)**2+rftp1(itp, &
                       & 2)**2+rftp1(itp, 3)**2)
                        rftp2 (itp, 1) = 0.5d0 * (t1+t2)
                        rftp2 (itp, 2) = 0.5d0 * (t1-t2)
                     End Do
                  Else
! collinear
                     Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                    & lmmaxvr, rhomt(:, ir, ias), 1, 0.d0, rftp1(:, 1), &
                    & 1)
                     Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                    & lmmaxvr, magmt(:, ir, ias, 1), 1, 0.d0, rftp1(:, &
                    & 2), 1)
                     Do itp = 1, lmmaxvr
                        t1 = rftp1 (itp, 1)
                        t2 = rftp1 (itp, 2)
                        rftp2 (itp, 1) = 0.5d0 * (t1+t2)
                        rftp2 (itp, 2) = 0.5d0 * (t1-t2)
                     End Do
                  End If
                  If (xcgrad .Le. 0) Then
                     Call xcifc (input%groundstate%xctypenumber, &
                    & n=lmmaxvr, rhoup=rftp2(:, 1), rhodn=rftp2(:, 2), &
                    & ex=ex, ec=ec, vxup=vx(:, 1), vxdn=vx(:, 2), &
                    & vcup=vc(:, 1), vcdn=vc(:, 2))
                  Else
                     Call xcifc (input%groundstate%xctypenumber, &
                    & n=lmmaxvr, rhoup=rftp2(:, 1), rhodn=rftp2(:, 2), &
                    & ex=ex, ec=ec, vxup=vx(:, 1), vxdn=vx(:, 2), &
                    & vcup=vc(:, 1), vcdn=vc(:, 2), grho=grhomt(:, ir), &
                    & gup=gupmt(:, ir), gdn=gdnmt(:, ir), &
                    & g2up=g2upmt(:, ir), g2dn=g2dnmt(:, ir), &
                    & g3rho=g3rhomt(:, ir), g3up=g3upmt(:, ir), &
                    & g3dn=g3dnmt(:, ir))
                  End If
                  If (ncmag) Then
! non-collinear: spin rotate the local exchange-correlation potential
                     Do itp = 1, lmmaxvr
                        t1 = vx (itp, 1) + vc (itp, 1)
                        t2 = vx (itp, 2) + vc (itp, 2)
                        t3 = 0.5d0 * (t1-t2)
! determine |m| again
                        t4 = 2.d0 * rftp2 (itp, 1) - rftp1 (itp, 4)
                        If (t4 .Gt. 1.d-8) t4 = t3 / t4
                        rftp1 (itp, 1:3) = rftp1 (itp, 1:3) * t4
                        rftp1 (itp, 4) = 0.5d0 * (t1+t2)
                     End Do
                     Do i = 1, 3
                        Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, &
                       & rfshtvr, lmmaxvr, rftp1(:, i), 1, 0.d0, &
                       & bxcmt(:, ir, ias, i), 1)
                     End Do
                     Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, &
                    & lmmaxvr, rftp1(:, 4), 1, 0.d0, vxcmt(:, ir, ias), &
                    & 1)
                  Else
! collinear
                     Do itp = 1, lmmaxvr
                        t1 = vx (itp, 1) + vc (itp, 1)
                        t2 = vx (itp, 2) + vc (itp, 2)
                        rftp1 (itp, 1) = 0.5d0 * (t1+t2)
                        rftp1 (itp, 2) = 0.5d0 * (t1-t2)
                     End Do
                     Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, &
                    & lmmaxvr, rftp1(:, 1), 1, 0.d0, vxcmt(:, ir, ias), &
                    & 1)
                     Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, &
                    & lmmaxvr, rftp1(:, 2), 1, 0.d0, bxcmt(:, ir, ias, &
                    & 1), 1)
                  End If
               Else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
                  Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rbshtvr, &
                 & lmmaxvr, rhomt(:, ir, ias), 1, 0.d0, rftp1, 1)
                  If (xcgrad .Le. 0) Then
                     Call xcifc (input%groundstate%xctypenumber, &
                    & n=lmmaxvr, rho=rftp1, ex=ex, ec=ec, vx=vx, vc=vc)
                  Else
                     Call xcifc (input%groundstate%xctypenumber, &
                    & n=lmmaxvr, rho=rftp1, ex=ex, ec=ec, vx=vx, vc=vc, &
                    & grho=gupmt(:, ir), g2rho=g2upmt(:, ir), &
                    & g3rho=g3upmt(:, ir))
                  End If
                  Do itp = 1, lmmaxvr
                     rftp1 (itp, 1) = vx (itp, 1) + vc (itp, 1)
                  End Do
! exchange-correlation potential
                  Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, &
                 & lmmaxvr, rftp1, 1, 0.d0, vxcmt(:, ir, ias), 1)
               End If
! convert energy densities from spherical coordinates to spherical harmonics
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, &
              & lmmaxvr, ex, 1, 0.d0, exmt(:, ir, ias), 1)
               Call dgemv ('N', lmmaxvr, lmmaxvr, 1.d0, rfshtvr, &
              & lmmaxvr, ec, 1, 0.d0, ecmt(:, ir, ias), 1)
            End Do
         End Do
      End Do
!--------------------------------!
!     interstitial potential     !
!--------------------------------!
! gradients for GGA if required
      If (xcgrad .Eq. 1) Call ggair (grhoir, gupir, gdnir, g2upir, &
     & g2dnir, g3rhoir, g3upir, g3dnir)
      If (associated(input%groundstate%spin)) Then
!------------------------!
!     spin-polarised     !
!------------------------!
         Do ir = 1, ngrtot
            If (ncmag) Then
! non-collinear
               t1 = rhoir (ir)
! compute |m|
               t2 = Sqrt (magir(ir, 1)**2+magir(ir, 2)**2+magir(ir, &
              & 3)**2)
            Else
! collinear
               t1 = rhoir (ir)
               t2 = magir (ir, 1)
            End If
            rfir (ir, 1) = 0.5d0 * (t1+t2)
            rfir (ir, 2) = 0.5d0 * (t1-t2)
         End Do
         If (xcgrad .Le. 0) Then
            Call xcifc (input%groundstate%xctypenumber, n=ngrtot, &
           & rhoup=rfir(:, 1), rhodn=rfir(:, 2), ex=exir, ec=ecir, &
           & vxup=vx(:, 1), vxdn=vx(:, 2), vcup=vc(:, 1), vcdn=vc(:, &
           & 2))
         Else
            Call xcifc (input%groundstate%xctypenumber, n=ngrtot, &
           & rhoup=rfir(:, 1), rhodn=rfir(:, 2), ex=exir, ec=ecir, &
           & vxup=vx(:, 1), vxdn=vx(:, 2), vcup=vc(:, 1), vcdn=vc(:, &
           & 2), grho=grhoir, gup=gupir, gdn=gdnir, g2up=g2upir, &
           & g2dn=g2dnir, g3rho=g3rhoir, g3up=g3upir, g3dn=g3dnir)
         End If
         If (ncmag) Then
! non-collinear: spin rotate the local exchange potential
            Do ir = 1, ngrtot
               t1 = vx (ir, 1) + vc (ir, 1)
               t2 = vx (ir, 2) + vc (ir, 2)
               t3 = 0.5d0 * (t1-t2)
! determine |m| again
               t4 = 2.d0 * rfir (ir, 1) - rhoir (ir)
               If (t4 .Gt. 1.d-8) t4 = t3 / t4
               bxcir (ir, :) = magir (ir, :) * t4
               vxcir (ir) = 0.5d0 * (t1+t2)
            End Do
         Else
! collinear
            Do ir = 1, ngrtot
               t1 = vx (ir, 1) + vc (ir, 1)
               t2 = vx (ir, 2) + vc (ir, 2)
               vxcir (ir) = 0.5d0 * (t1+t2)
               bxcir (ir, 1) = 0.5d0 * (t1-t2)
            End Do
         End If
      Else
!--------------------------!
!     spin-unpolarised     !
!--------------------------!
         If (xcgrad .Le. 0) Then
            Call xcifc (input%groundstate%xctypenumber, n=ngrtot, &
           & rho=rhoir, ex=exir, ec=ecir, vx=vx, vc=vc)
         Else
            Call xcifc (input%groundstate%xctypenumber, n=ngrtot, &
           & rho=rhoir, ex=exir, ec=ecir, vx=vx, vc=vc, grho=gupir, &
           & g2rho=g2upir, g3rho=g3upir)
         End If
         vxcir (1:ngrtot) = vx (1:ngrtot, 1) + vc (1:ngrtot, 1)
      End If
! optimised effective potential
      If (input%groundstate%xctypenumber .Lt. 0) Call oepmain
! symmetrise the exchange-correlation potential
      Call symrf (1, vxcmt, vxcir)
      If (associated(input%groundstate%spin)) Then
! remove the source contribution if required
         If (input%groundstate%nosource) Call projsbf
! symmetrise the exchange-correlation effective field
         Call symrvf (1, bxcmt, bxcir)
      End If
      Deallocate (ex, ec, vx, vc, rftp1, rftp2)
      If (associated(input%groundstate%spin)) deallocate (rfir)
      If (xcgrad .Eq. 1) Then
         Deallocate (grhomt, gupmt, gdnmt, g2upmt, g2dnmt, g3rhomt, &
        & g3upmt, g3dnmt)
         Deallocate (grhoir, gupir, gdnir, g2upir, g2dnir, g3rhoir, &
        & g3upir, g3dnir)
      End If
      Return
End Subroutine
!EOC
