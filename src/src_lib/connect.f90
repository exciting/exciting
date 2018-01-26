!
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: connect
! !INTERFACE:
!
!
Subroutine connect (cvec, plotdef, nv, np, vpl, dv, dp)
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
!   {\tt dp}, respectively.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!   Improved September 2007 (JKD)
!EOP
!BOC
      Use modinput
      Implicit None
! arguments
      Real (8), Intent (In) :: cvec (3, 3)
      Type (plot1d_type), Intent (In) :: plotdef
      Integer, Intent (In) :: nv
      Integer, Intent (In) :: np
      Real (8), Intent (Out) :: vpl (3, np)
      Real (8), Intent (Out) :: dv (nv)
      Real (8), Intent (Out) :: dp (np)
!
!
! local variables
!
      Real (8) :: vvl (3, size(plotdef%path%pointarray))
!
      Integer :: iv, ip, ip0, ip1, n
      Real (8) :: vl (3), vc (3)
      Real (8) :: dt, f, t1
! alloctable arrays
      Real (8), Allocatable :: seg (:)
!
      Do iv = 1, nv
         vvl (:, iv) = plotdef%path%pointarray(iv)%point%coord
      End Do
!
      If (nv .Lt. 1) Then
         Write (*,*)
         Write (*, '("Error(connect): nv < 1 : ", I8)') nv
         Write (*,*)
         Stop
      End If
      If (np .Lt. nv) Then
         Write (*,*)
         Write (*, '("Error(connect): np < nv : ", 2I8)') np, nv
         Write (*,*)
         Stop
      End If
      If (np .Eq. 1) Then
         vpl (:, 1) = vvl (:, 1)
         dv (1) = 0.d0
         dp (1) = 0.d0
         Return
      End If
      Allocate (seg(nv))
! find the total distance and the length of each segment
      dt = 0.d0
      Do iv = 1, nv - 1
         dv (iv) = dt
         vl (:) = vvl (:, iv+1) - vvl (:, iv)
         Call r3mv (cvec, vl, vc)
         seg (iv) = Sqrt (vc(1)**2+vc(2)**2+vc(3)**2)
         dt = dt + seg (iv)
      End Do
      dv (nv) = dt
      If (dt .Lt. 1.d-8) Then
         Do ip = 1, np
            vpl (:, ip) = vvl (:, 1)
            dp (ip) = 0.d0
         End Do
      Else
         dp( 1) = 0.d0
         vpl( :, 1) = vvl( :, 1)
         n = 1
         Do iv = 1, nv - 1
            t1 = (np-nv)*seg( iv)/dt 
            !t1 = dble (np) * dv (iv) / dt
            !ip0 = Nint (t1) + 1
            !If (ip0 .Lt. 1) ip0 = 1
            !t1 = dble (np) * dv (iv+1) / dt
            !ip1 = Nint (t1)
            !If (ip1 .Gt. np) ip1 = np
            !n = ip1 - ip0
            !If (n .Le. 0) n = 1
            ip1 = 1+nint( t1)
            if( iv .eq. nv-1) ip1 = np - n

            Do ip = 1, ip1
               f = dble( ip)/ip1
               n = n + 1
               dp( n) = f*seg( iv) + dv( iv)
               vpl( :, n) = vvl( :, iv)*(1.d0 - f) + vvl( :, iv+1)*f
            End Do
         End Do
      End If
      Deallocate (seg)
      Return
End Subroutine
!EOC
