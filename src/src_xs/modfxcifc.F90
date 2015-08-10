! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module modfxcifc
      Implicit None
Contains
!
!BOP
! !ROUTINE: fxcifc
! !INTERFACE:
!
!
      Subroutine fxcifc (fxctype, w, iw, iq, ng, oct, alrc, alrcd, &
     & blrcd, fxcg)
         Use m_fxc_lrc, Only: fxc_lrc
         Use m_fxc_lrcd, Only: fxc_lrcd
         Use m_fxc_alda, Only: fxc_alda
         Use m_fxc_bse_ma03, Only: fxc_bse_ma03
! !INPUT/OUTPUT PARAMETERS:
!   fxctype : type of exchange-correlation functional (in,integer)
! !DESCRIPTION:
!   Interface to the exchange-correlation kernel routines. This makes it
!   relatively
!   simple to add new functionals which do not necessarily depend only on
!   all input parameters. Based upon the routine {\tt modxcifc}.
!
! !REVISION HISTORY:
!   Created October 2007 (Sagmeister)
!EOP
!BOC
         Implicit None
    ! mandatory arguments
         Integer, Intent (In) :: fxctype
    ! optional arguments
         Complex (8), Optional, Intent (In) :: w
         Integer, Optional, Intent (In) :: iw
         Integer, Optional, Intent (In) :: iq
         Integer, Optional, Intent (In) :: ng
         Integer, Optional, Intent (In) :: oct
         Real (8), Optional, Intent (In) :: alrc
         Real (8), Optional, Intent (In) :: alrcd
         Real (8), Optional, Intent (In) :: blrcd
         Complex (8), Optional, Intent (Out) :: fxcg (:, :)
    ! local variables
    ! automatic arrays
         Select Case (Abs(fxctype))
         Case (0)
       ! RPA case fxc is zero
            If (present(fxcg) .And. present(ng)) Then
               fxcg = (0.d0, 0.d0)
            Else
               Go To 10
            End If
         Case (1)
       ! static long-range kernel without local field effects
       ! L. Reining, Phys. Rev. Lett. 88, 06404 (2002)
            If (present(fxcg) .And. present(ng) .And. present(alrc)) &
           & Then
               Call fxc_lrc (ng, .False., alrc, fxcg)
            Else
               Go To 10
            End If
         Case (2)
       ! static long-range kernel including local field effects
       ! L. Reining, Phys. Rev. Lett. 88, 06404 (2002)
            If (present(fxcg) .And. present(ng) .And. present(alrc)) &
           & Then
               Call fxc_lrc (ng, .True., alrc, fxcg)
            Else
               Go To 10
            End If
         Case (3)
       ! dynamical long-range kernel without local field effects
       ! L. Reining, Phys. Rev. B 72, 125203 (2005)
            If (present(fxcg) .And. present(ng) .And. present(alrcd) &
           & .And. present(w) .And. present(blrcd)) Then
               Call fxc_lrcd (ng, .False., alrcd, blrcd, w, fxcg)
            Else
               Go To 10
            End If
         Case (4)
       ! static long-range kernel including local field effects
       ! L. Reining, Phys. Rev. B 72, 125203 (2005)
            If (present(fxcg) .And. present(ng) .And. present(alrcd) &
           & .And. present(w) .And. present(blrcd)) Then
               Call fxc_lrcd (ng, .True., alrcd, blrcd, w, fxcg)
            Else
               Go To 10
            End If
         Case (5)
       ! ALDA kernel, [Reference]
            If (present(fxcg) .And. present(ng) .And. present(iq)) Then
               Call fxc_alda (iq, ng, fxcg)
            Else
               Go To 10
            End If
         Case (7)
       ! xc-kernel derived from the Bethe-Salpeter equation
       ! no local field effects
       ! A. Marini, Phys. Rev. Lett. 91, 256402 (2003)
            If (present(fxcg) .And. present(oct) .And. present(ng) &
           & .And. present(iw)) Then
               Call fxc_bse_ma03 (ng, oct, .False., iw, fxcg)
            Else
               Go To 10
            End If
         Case (8)
       ! xc-kernel derived from the Bethe-Salpeter equation
       ! inclusion of local field effects
       ! A. Marini, Phys. Rev. Lett. 91, 256402 (2003)
            If (present(fxcg) .And. present(ng) .And. present(oct) &
           & .And. present(iw)) Then
               Call fxc_bse_ma03 (ng, oct, .True., iw, fxcg)
            Else
               Go To 10
            End If
         Case Default
            Write (*,*)
            Write (*, '("Error(fxcifc): fxctype not defined : ", I8)') &
           & fxctype
            Write (*,*)
            Stop
         End Select
         Return
    ! error treatment
10       Continue
         Write (*,*)
         Write (*, '("Error(fxcifc): missing arguments for exchange-cor&
        &relation kernel type ", I5)') fxctype
         Write (*,*)
         Stop
      End Subroutine fxcifc
!EOC
!
!BOP
! !ROUTINE: getfxcdata
! !INTERFACE:
!
!
      Subroutine getfxcdata (fxctype, fxcdescr, fxcspin)
! !INPUT/OUTPUT PARAMETERS:
         Use modinput
!   fxctype  : type of exchange-correlation functional (in,integer)
!   fxcdescr : description of functional (out,character(256))
!   fxcspin  : spin treatment (out,integer)
!   fxcgrad  : gradient treatment (out,integer)
! !DESCRIPTION:
!   Returns data on the exchange-correlation functional labelled by
!   {\tt fxctype}. The character array {\tt fxctype} contains a short
!   description
!   of the functional including journal references. The variable
!   {\tt fxcspin} is
!   set to 1 or 0 for spin-polarised or -unpolarised functionals,
!   respectively.
!
! !REVISION HISTORY:
!   Created October 2007 (Sagmeister)
!EOP
!BOC
         Implicit None
         Integer, Intent (In) :: fxctype
         Character (256), Intent (Out) :: fxcdescr
         Integer, Intent (Out) :: fxcspin
         Select Case (Abs(fxctype))
         Case (0)
            fxcdescr = 'xc-kernel set to zero (RPA case)'
       ! spin-polarisation not required
            fxcspin = - 1
            Return
         Case (1)
            fxcdescr = 'long-range xc-kernel, no local field effects'
       ! spin-polarisation not required
            fxcspin = - 1
            Return
         Case (2)
            fxcdescr = 'long-range xc-kernel, including local field eff&
           &ects'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (3)
            fxcdescr = 'dynamical long-range xc-kernel, no local field &
           &effects'
       ! spin-polarisation not required
            fxcspin = - 1
            Return
         Case (4)
            fxcdescr = 'dynamical long-range xc-kernel, including local&
           & field effects'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (5)
            fxcdescr = 'ALDA kernel, including local field effects'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (7)
            fxcdescr = 'BSE kernel, A. Marini, Phys. Rev. Lett. 91, 256&
           &402 (2003), no local field effects'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (8)
            fxcdescr = 'BSE kernel, A. Marini, Phys. Rev. Lett. 91, 256&
           &402 (2003), including local field effects'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case Default
            Write (*,*)
            Write (*, '("Error(getfxcdata): fxctype not defined : ", I8&
           &)') fxctype
            Write (*,*)
            Stop
         End Select
      End Subroutine getfxcdata
!EOC
!
End Module modfxcifc
