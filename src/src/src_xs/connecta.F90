!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: connecta
! !INTERFACE:
!
!
Subroutine connecta (cvec, nv, np, vvl, vpl, dv, dp)
! !INPUT/OUTPUT PARAMETERS:
!   cvec : matrix of (reciprocal) lattice vectors stored column-wise
!         (in,real(3,3))
!   nv   : number of vertices (in,integer)
!   np   : number of connecting points (in,integer)
!   vvl  : vertex vectors in lattice coordinates (in,real(3,nv))
!   vpl  : connecting point vectors in lattice coordinates (out,real(3,np))
!   dv   : cummulative distance to each vertex (out,real(nv))
!   dp   : cummulative distance to each connecting point (out,real(np))
! !DESCRIPTION:
!   Generates a set of points which interpolate between a given set of vertices.
!   Vertex points are supplied in lattice coordinates in the array {\tt vvl} and
!   converted to Cartesian coordinates with the matrix {\tt cvec}. Interpolating
!   points are stored in the array {\tt vpl}. The cummulative distances to the
!   vertices and points along the path are stored in arrays {\tt dv} and
!   {\tt dp}, respectively. The given vertex points are contained in the set of
!   output vectors. Based upon the routine {\tt connect}.
!
! !REVISION HISTORY:
!   Created June 2007 (Sagmeister)
!EOP
!BOC
      Implicit None
! arguments
      Real (8), Intent (In) :: cvec (3, 3)
      Integer, Intent (In) :: nv
      Integer, Intent (In) :: np
      Real (8), Intent (In) :: vvl (3, nv)
      Real (8), Intent (Out) :: vpl (3, np)
      Real (8), Intent (Out) :: dv (nv)
      Real (8), Intent (Out) :: dp (np)
! local variables
      Integer :: iv, ip, j, c, spts, npi
      Real (8) :: st, vl (3), vc (3), v1 (3), v2 (3)
! alloctable arrays
      Real (8), Allocatable :: seg (:)
      Integer, Allocatable :: idx (:), segpts (:)
      If (nv .Lt. 1) Then
         Write (*,*)
         Write (*, '("Error(connecta): nv < 1 : ",I8)') nv
         Write (*,*)
         Stop
      End If
      If (np .Lt. nv) Then
         Write (*,*)
         Write (*, '("Error(connecta): np < nv : ",2I8)') np, nv
         Write (*,*)
         Stop
      End If
      If (np .Eq. 1) Then
         vpl (:, 1) = vvl (:, 1)
         dv (1) = 0.d0
         dp (1) = 0.d0
         Return
      End If
      Allocate (seg(nv-1), idx(nv-1), segpts(nv-1))
! find the total distance and the length of each segment
      st = 0.d0
      Do iv = 1, nv - 1
         dv (iv) = st
         vl (:) = vvl (:, iv+1) - vvl (:, iv)
         vc (:) = vl (1) * cvec (:, 1) + vl (2) * cvec (:, 2) + vl (3) &
        & * cvec (:, 3)
         seg (iv) = Sqrt (vc(1)**2+vc(2)**2+vc(3)**2)
         st = st + seg (iv)
      End Do
      dv (nv) = st
! sort segments according to their length in descending order
      Call sortidx (nv-1, seg, idx)
      npi = np - nv
      segpts (:) = 0
      Do ip = 1, npi
         iv = Mod (ip, nv-1)
         If (iv .Eq. 0) iv = nv - 1
         segpts (idx(nv-iv)) = segpts (idx(nv-iv)) + 1
      End Do
! loop over vertices
      c = 1
      Do iv = 1, nv - 1
         v1 (:) = vvl (:, iv)
         v2 (:) = vvl (:, iv+1)
         vpl (:, c) = v1 (:)
         dp (c) = dv (iv)
         c = c + 1
         spts = segpts (iv)
   ! linear interplation along segment
         Call linterplin (spts, v1, v2, vpl(1, c))
         dp (c:c+spts-1) = (/ (dp(c-1)+seg(iv)*dble(j)/dble(spts+1), &
        & j=1, spts) /)
         c = c + spts
         If (iv .Eq. nv-1) Then
            vpl (:, c) = v2 (:)
            dp (c) = dv (iv+1)
            c = c + 1
         End If
      End Do
      Deallocate (seg, idx, segpts)
End Subroutine
!EOC
!
!
Subroutine linterplin (np, v1, v2, vintp)
      Implicit None
  ! arguments
      Integer, Intent (In) :: np
      Real (8), Intent (In) :: v1 (3), v2 (3)
      Real (8), Intent (Out) :: vintp (3, np)
  ! local variables
      Real (8) :: np2, lam
      Integer :: j
      np2 = np + 2
      Do j = 1, np
         lam = dble (j) / (np2-1)
         vintp (:, j) = v1 (:) * (1-lam) + v2 (:) * lam
      End Do
End Subroutine linterplin
