!
!
!
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: genppts
! !INTERFACE:
!
!
Subroutine genppts (reducep, tfbz, ngridp, boxl, nppt, ipmap, ivp, vpl, &
& vpc, wppt)
! !USES:
      Use modinput
      Use modmain
#ifdef XS
      Use modxs
#endif
! !INPUT/OUTPUT PARAMETERS:
!   reducep : .true. if p-point set is to be reduced (in,logical)
!   tfbz    : .true. if vpl and vpc should be mapped to the first Brillouin
!             zone (in,logical)
!   ngridp  : p-point grid size (in,integer(3))
!   boxl    : corners of box containing p-points in lattice coordinates, the
!             first vector is the origin (in,real(3,4))
!   nppt    : total number of p-points (out,integer)
!   ipmap   : map from integer grid to p-point index
!             (out,integer(0:ngridp(1)-1,0:ngridp(2)-1,0:ngridp(3)-1))
!   ivp     : integer coordinates of the p-points
!             (out,integer(3,ngridp(1)*ngridp(2)*ngridp(3)))
!   vpl     : lattice coordinates of each p-point
!             (out,real(3,ngridp(1)*ngridp(2)*ngridp(3)))
!   vpc     : Cartesian coordinates of each p-point
!             (out,real(3,ngridp(1)*ngridp(2)*ngridp(3)))
!   wppt    : weights of each p-point (out,real(ngridp(1)*ngridp(2)*ngridp(3)))
! !DESCRIPTION:
!   This routine is used for generating $k$-point or $q$-point sets. Since these
!   are stored in global arrays, the points passed to this and other routines
!   are referred to as $p$-points. If {\tt reducep} is {\tt .true.} the set is
!   reduced with the spatial part of the crystal symmetries. In lattice
!   coordinates, the ${\bf p}$ vectors are given by
!   $$ {\bf p}=\left(\begin{matrix} & & \\
!     {\bf B}_2-{\bf B}_1 & {\bf B}_3-{\bf B}_1 & {\bf B}_4-{\bf B}_1 \\
!       & & \end{matrix}\right)
!     \left(\begin{matrix}i_1/n_1 \\ i_2/n_2 \\ i_3/n_3 \end{matrix}\right)
!     +{\bf B}_1 $$
!   where $i_j$ runs from 0 to $n_j-1$, and the ${\bf B}$ vectors define the
!   corners of a box with ${\bf B}_1$ as the origin. If {\tt tfbz} is
!   {\tt .true.} then the vectors {\tt vpl} (and {\tt vpc}) are mapped to the
!   first Brillouin zone. If {\tt tfbz} is {\tt .false.} and {\tt reducep} is
!   {\tt .true.} then the coordinates of {\tt vpl} are mapped to the $[0,1)$
!   interval. The $p$-point weights are stored in {\tt wppt} and the array
!   {\tt ipmap} contains the map from the integer coordinates to the reduced
!   index.
!
! !REVISION HISTORY:
!   Created August 2002 (JKD)
!   Updated April 2007 (JKD)
!   Modifications for excited states, November 2007 (Sagmeister)
!   Added mapping to the first Brillouin zone, September 2008 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: reducep
      Logical, Intent (In) :: tfbz
      Integer, Intent (In) :: ngridp (3)
      Real (8), Intent (In) :: boxl (3, 4)
      Integer, Intent (Out) :: nppt
      Integer, Intent (Out) :: ipmap (0:ngridp(1)-1, 0:ngridp(2)-1, &
     & 0:ngridp(3)-1)
      Integer, Intent (Out) :: ivp (3, ngridp(1)*ngridp(2)*ngridp(3))
      Real (8), Intent (Out) :: vpl (3, ngridp(1)*ngridp(2)*ngridp(3))
      Real (8), Intent (Out) :: vpc (3, ngridp(1)*ngridp(2)*ngridp(3))
      Real (8), Intent (Out) :: wppt (ngridp(1)*ngridp(2)*ngridp(3))
! local variables
      Integer :: i1, i2, i3, ip, jp
      Integer :: isym, lspl, iv (3)
      Real (8) :: v1 (3), v2 (3), v3 (3)
      Real (8) :: b (3, 3), s (3, 3), t1, t2
#ifdef XS
      Integer :: jsym, nsymcrys_, lsplsymc_ (maxsymcrys), lsplsymct &
     & (maxsymcrys)
  ! use symmetries of little group of q
      If ((iqcu .Ne. 0) .And. reducep) Then
         If (nsymcrys .Ne. nsymcrysq(iqcu)) Then
            Write (*, '(a)') 'Info(genppts): using associated(input%str&
           &ucture%symmetries) of the (little/small) group of q only'
         End If
     ! save global variables
         nsymcrys_ = nsymcrys
         lsplsymc_ (:) = lsplsymc (:)
     ! map to point group elements
         lsplsymct (:) = 0
         jsym = 0
         Do isym = 1, nsymcrysq (iqcu)
            jsym = jsym + 1
            lsplsymct (jsym) = lsplsymc (scqmap(isym, iqcu))
         End Do
     ! update global variables
         nsymcrys = nsymcrysq (iqcu)
         lsplsymc (:) = lsplsymct (:)
      End If
  ! now we are working with the point group symmetries of the small group of q
#endif
      If ((ngridp(1) .Le. 0) .Or. (ngridp(2) .Le. 0) .Or. (ngridp(3) &
     & .Le. 0)) Then
         Write (*,*)
         Write (*, '("Error(genppts): invalid ngridp : ", 3I8)') ngridp
         Write (*,*)
         Stop
      End If
! box vector matrix
      b (:, 1) = boxl (:, 2) - boxl (:, 1)
      b (:, 2) = boxl (:, 3) - boxl (:, 1)
      b (:, 3) = boxl (:, 4) - boxl (:, 1)
      t1 = 1.d0 / dble (ngridp(1)*ngridp(2)*ngridp(3))
      ip = 0
      Do i3 = 0, ngridp (3) - 1
         v1 (3) = dble (i3) / dble (ngridp(3))
         Do i2 = 0, ngridp (2) - 1
            v1 (2) = dble (i2) / dble (ngridp(2))
            Do i1 = 0, ngridp (1) - 1
               v1 (1) = dble (i1) / dble (ngridp(1))
               Call r3mv (b, v1, v2)
               v2 (:) = v2 (:) + boxl (:, 1)
               If (reducep) Then
                  Call r3frac (input%structure%epslat, v2, iv)
! determine if this point is equivalent to one already in the set
                  Do isym = 1, nsymcrys
                     lspl = lsplsymc (isym)
                     s (:, :) = dble (symlat(:, :, lspl))
                     Call r3mtv (s, v2, v3)
                     Call r3frac (input%structure%epslat, v3, iv)
                     Do jp = 1, ip
                        t2 = Abs (vpl(1, jp)-v3(1)) + Abs (vpl(2, &
                       & jp)-v3(2)) + Abs (vpl(3, jp)-v3(3))
                        If (t2 .Lt. input%structure%epslat) Then
! equivalent k-point found so add to current weight
                           ipmap (i1, i2, i3) = jp
                           wppt (jp) = wppt (jp) + t1
                           Go To 10
                        End If
                     End Do
                  End Do
               End If
! add new point to set
               ip = ip + 1
               ipmap (i1, i2, i3) = ip
               ivp (1, ip) = i1
               ivp (2, ip) = i2
               ivp (3, ip) = i3
               vpl (:, ip) = v2 (:)
               wppt (ip) = t1
10             Continue
            End Do
         End Do
      End Do
      nppt = ip
      Do ip = 1, nppt
! map vpl to the first Brillouin zone if required
         If (tfbz) Call vecfbz (input%structure%epslat, bvec, vpl(:, &
        & ip), iv)
! determine the Cartesian coordinates of the p-points
         Call r3mv (bvec, vpl(:, ip), vpc(:, ip))
      End Do
#ifdef XS
      If ((iqcu .Ne. 0) .And. reducep) Then
     ! restore global variables
         nsymcrys = nsymcrys_
         lsplsymc (:) = lsplsymc_ (:)
      End If
#endif
      Return
End Subroutine
!EOC
