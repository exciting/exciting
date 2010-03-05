
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine hamiltonandoverlapsetup (system, ngp, apwalm, igpig, vgpc)
      Use modfvsystem
      Use modinput
      Use mod_eigensystem
      Use mod_atoms
      Use mod_timing
      Use mod_muffin_tin
      Use mod_APW_LO
      Use mod_gkvector
!
      Implicit None
      Type (evsystem) :: system
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Integer :: n
      Character (256) :: prefix
!local variables
      Integer, Save :: ikc
      Real (8), Save :: cputot
      Real (8) :: cpuaa, cpualo, cpulolo, cpui, cpu00, cpu01
      Integer :: i, is, ia
      Complex (8) v (1)
      Real (8) :: cpu0, cpu1
      Real (8) :: threshold
!----------------------------------------!
!     Hamiltonian and overlap set up     !
!----------------------------------------!
!
!
      Call timesec (cpu0)
! set the matrices to zero
!
! muffin-tin contributions
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            Call hmlaan (system%hamilton, is, ia, ngp, apwalm)
            Call hmlalon (system%hamilton, is, ia, ngp, apwalm)
            Call hmllolon (system%hamilton, is, ia, ngp)
            Call olpaan (system%overlap, is, ia, ngp, apwalm)
            Call olpalon (system%overlap, is, ia, ngp, apwalm)
            Call olplolon (system%overlap, is, ia, ngp)
         End Do
      End Do
!
! interstitial contributions
      Call hmlistln (system%hamilton, ngp, igpig, vgpc)
      Call olpistln (system%overlap, ngp, igpig)
      threshold = 1e-16
!call HermiteanMatrixTruncate(system%hamilton,threshold)
!call HermiteanMatrixTruncate(system%overlap,threshold)
!
!
!
      If ( .Not. ispacked(system%hamilton)) Then
         Call hamiltonoverlapocopy_UL (system)
      End If
#ifdef DEBUGHO
      Write (*,*) "apwalm", apwalm
      prefix = "H"
      Call HermiteanMatrixToFiles (system%hamilton, prefix)
      prefix = "O"
      Call HermiteanMatrixToFiles (system%overlap, prefix)
      Write (*,*) "wrote"
      Stop
#endif
!
      Call timesec (cpu1)
      timemat = timemat + cpu1 - cpu0
!
!
!write(60,*)
!write(60,'("Muffin-tin Hamiltonian setup; Timings (CPU seconds) :")')
!write(60,'(" k-point",T40,": ",I8)') ikc
!write(60,'(" APW-APW",T40,": ",F12.2)') cpuaa
!write(60,'(" APW- lo",T40,": ",F12.2)') cpualo
!write(60,'(" lo - lo",T40,": ",F12.2)') cpulolo
!write(60,'(" interstitial",T40,": ",F12.2)') cpui
!write(60,'(" total",T40,": ",F12.2)') cpuaa+cpualo+cpulolo+cpui
!cputot=cputot+cpuaa+cpualo+cpulolo+cpui
!write(60,'(" cumulative total",T40,": ",F12.2)') cputot
!write(60,*)
!
!
!
End Subroutine
