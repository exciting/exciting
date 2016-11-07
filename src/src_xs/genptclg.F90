!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genptclg (cuttype, ngpmax, ngp, vgpc, gpc, sptcl)
      Use modmain
      Use modxs
      Implicit None
  ! arguments
      Character (*), Intent (In) :: cuttype
      Integer, Intent (In) :: ngpmax, ngp
      Real (8), Intent (In) :: vgpc (3, ngpmax), gpc (ngpmax)
      Real (8), Intent (Out) :: sptcl (ngpmax)
  ! local variables
      Integer :: igp
  ! external functions
      Real (8), External :: ptclg
      Do igp = 1, ngp
         sptcl (igp) = ptclg (cuttype, vgpc(:, igp), gpc(igp))
      End Do
End Subroutine genptclg
!
!
Real (8) Function ptclg (cuttype, vgpc, gpc)
      Use modmain, Only: fourpi
      use modmpi
      Implicit None
  ! arguments
      Character (*), Intent (In) :: cuttype
      Real (8), Intent (In) :: vgpc (3), gpc
  ! local variables
      real(8), parameter :: eps=1.d-8
      Real (8) :: t1
      Select Case (cuttype)
      Case ('nocutoff')
     ! set up the square root of the Coulomb potential from analytical
     ! expression (no cutoff)
         ptclg = 0.d0
         if (gpc .gt. eps) then
           ptclg = Sqrt (fourpi) / gpc
         else
           ! set sqrt of Coulomb potential to zero for |G+q| = 0
           ptclg = 0.d0
         end if
      Case ('0d')
     ! 0D spherical cutoff
         t1 = vgpc (1)
         Write (*,*)
         Write (*, '("Error(genptclg): 0D cutoff to be implemented")')
         Write (*,*)
         Call terminate
      Case ('1d')
     ! 1D infinite cylinder
         Write (*,*)
         Write (*, '("Error(genptclg): 1D cutoff to be implemented")')
         Write (*,*)
         Call terminate
      Case ('2d')
     ! 2D infinite slab
         Write (*,*)
         Write (*, '("Error(genptclg): 2D cutoff to be implemented")')
         Write (*,*)
         Call terminate
      Case Default
         Write (*,*)
         Write (*, '("Error(genptclg): unknown cutoff type for Coulomb &
        &potential: ", a)') cuttype
         Write (*,*)
         Call terminate
      End Select
!
End Function ptclg
