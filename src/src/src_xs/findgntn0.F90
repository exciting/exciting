!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_findgntn0
      Implicit None
      Save
!
  ! shapes
      Integer :: l1shape
      Integer, Allocatable :: m1shape (:), l2shape (:, :), m2shape (:, &
     & :, :)
      Integer, Allocatable :: l3shape (:, :, :, :), m3shape (:, :, :, &
     & :, :)
!
  ! maps
      Integer, Allocatable :: l1map (:), m1map (:, :), l2map (:, :, :), &
     & m2map (:, :, :, :)
      Integer, Allocatable :: l3map (:, :, :, :, :), m3map (:, :, :, :, &
     & :, :)
!
  ! number of non-zero Gaunt coefficients
      Integer :: ngauntnz
  ! number of Gaunt coefficients
      Integer :: ngaunt
!
Contains
!
!
      Subroutine findgntn0 (lmax1, lmax2, lmax3, gntc)
!
         Implicit None
!
         Character (*), Parameter :: thisnam = 'findgntn0'
!
    ! maximum l values for checking Gaunts
         Integer, Intent (In) :: lmax1, lmax2, lmax3
!
    ! array containing Gaunt coefficients like gntc(lm1,lm2,lm3)
    ! where lm(l,m) runs from l=0,lmax; m=-l,l; m fastest!
         Real (8), Intent (In) :: gntc (:, :, :)
!
    ! local variables
         Integer :: l1, m1, lm1, l2, m2, lm2, l3, m3, lm3, l, m, lmax, &
        & lmmax, co, coo
         Integer :: cl1, cm1, cl2, cm2, cl3, cm3
         Integer :: l1len, m1len, l2len, m2len, l3len, m3len
         Real (8), Parameter :: eps = 1.d-20
         Integer, Allocatable :: idxlm (:, :)
         Logical, Allocatable :: l1mask (:), m1mask (:, :)
         Logical, Allocatable :: l2mask (:, :, :), m2mask (:, :, :, :)
         Logical, Allocatable :: l3mask (:, :, :, :, :), m3mask (:, :, &
        & :, :, :, :)
!
         lmax = Max (lmax1, lmax2, lmax3)
         lmmax = (lmax+1) ** 2
!
         Allocate (idxlm(0:lmax,-lmax:lmax))
!
         co = 0
         Do l = 0, lmax
            Do m = - l, l
               co = co + 1
               idxlm (l, m) = co
            End Do
         End Do
!
         Allocate (l1mask(0:lmax1))
         Allocate (m1mask(0:lmax1,-lmax1:lmax1))
         Allocate (l2mask(0:lmax1,-lmax1:lmax1, 0:lmax2))
         Allocate (m2mask(0:lmax1,-lmax1:lmax1, 0:lmax2,-lmax2:lmax2))
         Allocate (l3mask(0:lmax1,-lmax1:lmax1, 0:lmax2,-lmax2:lmax2, &
        & 0:lmax3))
         Allocate (m3mask(0:lmax1,-lmax1:lmax1, 0:lmax2,-lmax2:lmax2, &
        & 0:lmax3,-lmax3:lmax3))
!
         l1mask (:) = .False.
         m1mask (:, :) = .False.
         l2mask (:, :, :) = .False.
         m2mask (:, :, :, :) = .False.
         l3mask (:, :, :, :, :) = .False.
         m3mask (:, :, :, :, :, :) = .False.
!
    !
    ! second part of loops
    !
!
         co = 0
         coo = 0
         Do l1 = 0, lmax1
            Do m1 = - l1, l1
               lm1 = idxlm (l1, m1)
               Do l2 = 0, lmax2
                  Do m2 = - l2, l2
                     lm2 = idxlm (l2, m2)
                     Do l3 = 0, lmax3
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
!
                           coo = coo + 1
                           If (Abs(gntc(lm1, lm3, lm2)) > eps) Then ! *** lm2 <-> lm3
                              co = co + 1
                              l1mask (l1) = .True.
                              m1mask (l1, m1) = .True.
                              l2mask (l1, m1, l2) = .True.
                              m2mask (l1, m1, l2, m2) = .True.
                              l3mask (l1, m1, l2, m2, l3) = .True.
                              m3mask (l1, m1, l2, m2, l3, m3) = .True.
                           End If
!
                        End Do ! m3
                     End Do ! l3
                  End Do ! m2
               End Do ! l2
            End Do ! m1
         End Do ! l1
!
         ngaunt = coo
         ngauntnz = co
!
         If (allocated(m1shape)) deallocate (m1shape)
         If (allocated(l2shape)) deallocate (l2shape)
         If (allocated(m2shape)) deallocate (m2shape)
         If (allocated(l3shape)) deallocate (l3shape)
         If (allocated(m3shape)) deallocate (m3shape)
         Allocate (m1shape(0:lmax1))
         Allocate (l2shape(0:lmax1,-lmax1:lmax1))
         Allocate (m2shape(0:lmax1,-lmax1:lmax1, 0:lmax2))
         Allocate (l3shape(0:lmax1,-lmax1:lmax1, 0:lmax2,-lmax2:lmax2))
         Allocate (m3shape(0:lmax1,-lmax1:lmax1, 0:lmax2,-lmax2:lmax2, &
        & 0:lmax3))
!
    !
    ! shape arrays of counting variables
    !
!
         l1shape = count (l1mask, 1)
         m1shape = count (m1mask, 2)
         l2shape = count (l2mask, 3)
         m2shape = count (m2mask, 4)
         l3shape = count (l3mask, 5)
         m3shape = count (m3mask, 6)
!
    !
    ! dimensions of maps
    !
!
         l1len = l1shape
         m1len = maxval (m1shape)
         l2len = maxval (l2shape)
         m2len = maxval (m2shape)
         l3len = maxval (l3shape)
         m3len = maxval (m3shape)
!
         If (allocated(l1map)) deallocate (l1map)
         If (allocated(m1map)) deallocate (m1map)
         If (allocated(l2map)) deallocate (l2map)
         If (allocated(m2map)) deallocate (m2map)
         If (allocated(l3map)) deallocate (l3map)
         If (allocated(m3map)) deallocate (m3map)
         Allocate (l1map(l1len))
         Allocate (m1map(0:lmax1, m1len))
         Allocate (l2map(0:lmax1,-lmax1:lmax1, l2len))
         Allocate (m2map(0:lmax1,-lmax1:lmax1, 0:lmax2, m2len))
         Allocate (l3map(0:lmax1,-lmax1:lmax1, 0:lmax2,-lmax2:lmax2, &
        & l3len))
         Allocate (m3map(0:lmax1,-lmax1:lmax1, 0:lmax2,-lmax2:lmax2, &
        & 0:lmax3, m3len))
!
         cl1 = 0
         Do l1 = 0, lmax1
!
            If (l1mask(l1)) Then
               cl1 = cl1 + 1
               l1map (cl1) = l1
            End If
            cm1 = 0
!
            Do m1 = - l1, l1
!
               If (m1mask(l1, m1)) Then
                  cm1 = cm1 + 1
                  m1map (l1, cm1) = m1
               End If
               cl2 = 0
!
               Do l2 = 0, lmax2
!
                  If (l2mask(l1, m1, l2)) Then
                     cl2 = cl2 + 1
                     l2map (l1, m1, cl2) = l2
                  End If
                  cm2 = 0
!
                  Do m2 = - l2, l2
!
                     If (m2mask(l1, m1, l2, m2)) Then
                        cm2 = cm2 + 1
                        m2map (l1, m1, l2, cm2) = m2
                     End If
                     cl3 = 0
!
                     Do l3 = 0, lmax3
!
                        If (l3mask(l1, m1, l2, m2, l3)) Then
                           cl3 = cl3 + 1
                           l3map (l1, m1, l2, m2, cl3) = l3
                        End If
                        cm3 = 0
!
                        Do m3 = - l3, l3
!
                           If (m3mask(l1, m1, l2, m2, l3, m3)) Then
                              cm3 = cm3 + 1
                              m3map (l1, m1, l2, m2, l3, cm3) = m3
                           End If
!
                        End Do ! cm3
                     End Do ! cl3
                  End Do ! cm2
               End Do ! cl2
            End Do ! cm1
         End Do ! cl1
!
         Deallocate (l1mask, m1mask, l2mask, m2mask, l3mask, m3mask)
         Deallocate (idxlm)
!
      End Subroutine findgntn0
!
!
      Subroutine findgntn0_clear
         Implicit None
         l1shape = 0
         If (allocated(m1shape)) deallocate (m1shape)
         If (allocated(l2shape)) deallocate (l2shape)
         If (allocated(m2shape)) deallocate (m2shape)
         If (allocated(l3shape)) deallocate (l3shape)
         If (allocated(m3shape)) deallocate (m3shape)
         If (allocated(l1map)) deallocate (l1map)
         If (allocated(m1map)) deallocate (m1map)
         If (allocated(l2map)) deallocate (l2map)
         If (allocated(m2map)) deallocate (m2map)
         If (allocated(l3map)) deallocate (l3map)
         If (allocated(m3map)) deallocate (m3map)
      End Subroutine findgntn0_clear
!
!
End Module m_findgntn0
