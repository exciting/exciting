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
      Subroutine fxcifc (fxctype, nw, w, iw, iq, ng, ngtot, oct, alrc, &
      & alrcd, blrcd, bootstrapit, kernsca, kerndiag, chim, chim_w, &
      & fxcg, fxcg_w)
         Use m_fxc_lrc, Only: fxc_lrc
         Use m_fxc_lrcd, Only: fxc_lrcd
         Use m_fxc_alda, Only: fxc_alda
         Use m_fxc_bse_ma03, Only: fxc_bse_ma03
         Use m_fxc_bootstrap, Only: fxc_bootstrap
         Use m_fxc_bootstrap_sc, Only: fxc_bootstrap_sc
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
         Logical, Optional, Intent (In) :: bootstrapit
         Logical, Optional, Intent (In) :: kernsca
         Logical, Optional, Intent (In) :: kerndiag
         Integer, Optional, Intent (In) :: nw
         Integer, Optional, Intent (In) :: iw
         Integer, Optional, Intent (In) :: iq
         Integer, Optional, Intent (In) :: ng
         Integer, Optional, Intent (In) :: ngtot
         Integer, Optional, Intent (In) :: oct
         Real (8), Optional, Intent (In) :: alrc
         Real (8), Optional, Intent (In) :: alrcd
         Real (8), Optional, Intent (In) :: blrcd
         Complex (8), Optional, Intent (In) :: chim (:, :)
         Complex (8), Optional, Intent (In) :: chim_w (:, :, :)
         Complex (8), Optional, Intent (Out) :: fxcg (:, :)
         Complex (8), Optional, Intent (Out) :: fxcg_w (:, :, :)
    ! local variables
        Logical :: ks, kd, bi, dyn
        Integer :: i
    ! automatic arrays
        If(present(kernsca)) then
            ks = kernsca
        Else
            ks = .false.
        End If
        If(present(kerndiag)) then
            kd = kerndiag
        Else
            kd = .false.
        End If
        If(present(bootstrapit)) then
            bi = bootstrapit
        Else
            bi = .false.
        End If
        If(present(fxcg_w)) then
            dyn = .true.
        Else
            dyn = .false.
        End If

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
          Case (9:18,20)

       ! bootstrap kernel, Phys. Rev. Lett. 107, 186401 (2011)
            If (present(ng)) Then
                If (.not.bi) then
                   If(dyn) then

                        Call fxc_bootstrap (fxctype, ng, &
                        & ngtot=ngtot, nw=nw, chim_w=chim_w, &
                        & fxc_w=fxcg_w)
                   else
                        Call fxc_bootstrap (fxctype, ng, &
                        & chim=chim, fxc=fxcg)
                   End If
                Else
                   Call fxc_bootstrap_sc(fxctype,.true., ng, chim, fxcg)
                End If
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

        If ((ks .or. kd) .and. (.not.dyn))  Call cropfxc(fxcg, ng, ks ,kd)
        If ((ks .or. kd) .and. dyn) Then
            do i = 1,nw
                Call cropfxc(fxcg_w(:,:,i), ng, ks ,kd)
            end do
        end if

        if(kd .and. ((fxctype.eq.13).or.(fxctype.eq.14))) Then
            if(ng .gt. 1) then
                do i = 2, ng
                    fxcg(i,i) = fxcg(1,1)
                end do
            end if
        end if

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
         Case (9)
            fxcdescr = 'Static bootstrap xc-kernel'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (10)
            fxcdescr = 'Static first-order bootstrap xc-kernel'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (11)
            fxcdescr = 'Static zero-order bootstrap xc-kernel'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (12)
            fxcdescr = 'Static first-order bootstrap xc-kernel with RPA&
            & df'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (13)
            fxcdescr = 'Static scalar first-order bootstrap xc-kernel &
            &with modified chi'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (14)
            fxcdescr = 'Static scalar first-order bootstrap xc-kernel &
            &with modified chi and RPA df'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (15)
            fxcdescr = 'Static zero-order bootstrap xc-kernel. &
            &Head component from NLF formula, rest of the fxc matrix&
            &built with full LF formula. This tries to correctly match&
            & both the position and hight of the excitonic peak.'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (16)
            fxcdescr = 'Dynamic zero-order bootstrap xc-kernel. &
            &Head component from NLF formula plus LF correction&
            &, rest of the fxc matrix&
            &built with full LF formula. This tries to correctly match&
            & both the position and hight of the excitonic peak.'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (17)
            fxcdescr = 'Matrix version of the fmr kernel: &
            &The reciprocal of the macroscopic dielectric function&
            &is replaced by the full rpa inverse dielectric matrix.'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (18)
            fxcdescr = 'Same as the original bootstrap but with the&
            &dielectric function replaced by its rpa version.'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (19)
            fxcdescr = 'Mimics xc-kernel. &
            &Head component from NLF formula plus LF correction&
            &, rest of the fxc matrix&
            &built with full LF formula. This tries to correctly match&
            & both the position and hight of the excitonic peak.'
       ! spin-polarisation not required
            fxcspin = 0
            Return
         Case (20)
            fxcdescr = 'Scalar bootstrap kernel. &
            &fxc = 1/(eps_M * \chi_0,00).'
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

!EOC
!
!BOP
! !ROUTINE: getfxcdata
! !INTERFACE:
!
!
      Subroutine cropfxc (fxcg, ng, kernsca, kerndiag)
! !INPUT/OUTPUT PARAMETERS:
         Use modinput
!   fxcg : exchange correlation kernel
! !DESCRIPTION:
!   Crops the xc kernel matrix according to value of attributes
!   kerndiag (set to zero all non-diagonal matrix elements of fxc)
!   or kernsca (set to zero all non-head matrix elements of fxc).
! !REVISION HISTORY:
!   Created August 2013 (SR)
!EOP
!BOC
         Implicit None
         Complex (8), Intent (InOut) :: fxcg (:, :)
         Integer, Intent (In) :: ng
         Logical, Intent (In) :: kernsca
         Logical, Intent (In) :: kerndiag
         Integer :: i, j

         if(kernsca) then
            write(*,*) 'making fxc scalar'
            do i = 1, ng
                do j = 1, ng
                    if ((i .ne. 1) .or. (j .ne. 1)) then
                        fxcg(i,j) = 0.d0
                    end if
                end do
            end do
         else if(kerndiag) then
            write(*,*) 'making fxc diagonal'
            do i = 1, ng
                do j = 1, ng
                    if (i .ne. j) then
                        fxcg(i,j) = 0.d0
                    end if
                end do
            end do
         end if

      End Subroutine cropfxc
!EOC
!
End Module modfxcifc
