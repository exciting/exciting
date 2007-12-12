
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_findgntn0
  implicit none
  save

  ! shapes
  integer :: l1shape
  integer, allocatable :: m1shape(:), l2shape(:,:), m2shape(:,:,:)
  integer, allocatable :: l3shape(:,:,:,:), m3shape(:,:,:,:,:)

  ! maps
  integer, allocatable :: l1map(:), m1map(:,:), l2map(:,:,:), m2map(:,:,:,:)
  integer, allocatable :: l3map(:,:,:,:,:), m3map(:,:,:,:,:,:)

  ! number of non-zero Gaunt coefficients
  integer :: ngauntnz
  ! number of Gaunt coefficients
  integer :: ngaunt

contains

  subroutine findgntn0 ( lmax1, lmax2, lmax3, gntc )

    implicit none

    character(*), parameter :: thisnam = 'findgntn0'

    ! maximum l values for checking Gaunts
    integer, intent(in) :: lmax1, lmax2, lmax3

    ! array containing Gaunt coefficients like gntc(lm1,lm2,lm3)
    ! where lm(l,m) runs from l=0,lmax; m=-l,l; m fastest!
    real(8), intent(in) :: gntc(:,:,:)

    ! local variables
    integer :: l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,l,m,lmax,lmmax,co,coo
    integer :: cl1,cm1,cl2,cm2,cl3,cm3
    integer :: l1len,m1len,l2len,m2len,l3len,m3len
    real(8), parameter :: eps = 1.d-20
    integer, allocatable :: idxlm(:,:)
    logical, allocatable :: l1mask(:), m1mask(:,:)
    logical, allocatable :: l2mask(:,:,:), m2mask(:,:,:,:)
    logical, allocatable :: l3mask(:,:,:,:,:), m3mask(:,:,:,:,:,:)

    lmax = max(lmax1, lmax2, lmax3)
    lmmax = (lmax+1)**2

    allocate(idxlm(0:lmax,-lmax:lmax))

    co=0
    do l=0,lmax
       do m=-l,l
          co = co + 1
          idxlm(l,m) = co
       end do
    end do

    allocate(l1mask(0:lmax1))
    allocate(m1mask(0:lmax1,-lmax1:lmax1))
    allocate(l2mask(0:lmax1,-lmax1:lmax1,0:lmax2))
    allocate(m2mask(0:lmax1,-lmax1:lmax1,0:lmax2,-lmax2:lmax2))
    allocate(l3mask(0:lmax1,-lmax1:lmax1,0:lmax2,-lmax2:lmax2,0:lmax3))
    allocate(m3mask(0:lmax1,-lmax1:lmax1,0:lmax2,-lmax2:lmax2,0:lmax3,&
         -lmax3:lmax3))

    l1mask(:) = .false.
    m1mask(:,:) = .false.
    l2mask(:,:,:) = .false.
    m2mask(:,:,:,:) = .false.
    l3mask(:,:,:,:,:) = .false.
    m3mask(:,:,:,:,:,:) = .false.

    !
    ! second part of loops
    !

    co = 0
    coo = 0
    do l1=0,lmax1
       do m1=-l1,l1
          lm1=idxlm(l1,m1)
          do l2=0,lmax2
             do m2=-l2,l2
                lm2=idxlm(l2,m2)
                do l3=0,lmax3
                   do m3=-l3,l3
                      lm3=idxlm(l3,m3)

                      coo = coo + 1
                      if (abs(gntc(lm1,lm3,lm2)) > eps) then ! *** lm2 <-> lm3
                         co = co + 1
                         l1mask(l1) = .true.
                         m1mask(l1,m1) = .true.
                         l2mask(l1,m1,l2) = .true.
                         m2mask(l1,m1,l2,m2) = .true.
                         l3mask(l1,m1,l2,m2,l3) = .true.
                         m3mask(l1,m1,l2,m2,l3,m3) = .true.
                      end if

                   end do ! m3
                end do ! l3
             end do ! m2
          end do ! l2
       end do ! m1
    end do ! l1

    ngaunt = coo
    ngauntnz = co

    if (allocated(m1shape)) deallocate(m1shape)
    if (allocated(l2shape)) deallocate(l2shape)
    if (allocated(m2shape)) deallocate(m2shape)
    if (allocated(l3shape)) deallocate(l3shape)
    if (allocated(m3shape)) deallocate(m3shape)
    allocate(m1shape(0:lmax1))
    allocate(l2shape(0:lmax1,-lmax1:lmax1))
    allocate(m2shape(0:lmax1,-lmax1:lmax1,0:lmax2))
    allocate(l3shape(0:lmax1,-lmax1:lmax1,0:lmax2,-lmax2:lmax2))
    allocate(m3shape(0:lmax1,-lmax1:lmax1,0:lmax2,-lmax2:lmax2,0:lmax3))

    !
    ! shape arrays of counting variables
    !

    l1shape = count(l1mask,1)
    m1shape = count(m1mask,2)
    l2shape = count(l2mask,3)
    m2shape = count(m2mask,4)
    l3shape = count(l3mask,5)
    m3shape = count(m3mask,6)

    !
    ! dimensions of maps
    !

    l1len = l1shape
    m1len = maxval(m1shape)
    l2len = maxval(l2shape)
    m2len = maxval(m2shape)
    l3len = maxval(l3shape)
    m3len = maxval(m3shape)

    if (allocated(l1map)) deallocate(l1map)
    if (allocated(m1map)) deallocate(m1map)
    if (allocated(l2map)) deallocate(l2map)
    if (allocated(m2map)) deallocate(m2map)
    if (allocated(l3map)) deallocate(l3map)
    if (allocated(m3map)) deallocate(m3map)
    allocate(l1map(l1len))
    allocate(m1map(0:lmax1,m1len))
    allocate(l2map(0:lmax1,-lmax1:lmax1,l2len))
    allocate(m2map(0:lmax1,-lmax1:lmax1,0:lmax2,m2len))
    allocate(l3map(0:lmax1,-lmax1:lmax1,0:lmax2,-lmax2:lmax2,l3len))
    allocate(m3map(0:lmax1,-lmax1:lmax1,0:lmax2,-lmax2:lmax2,0:lmax3,m3len))

    cl1 = 0
    do l1 =0, lmax1

       if (l1mask(l1)) then
          cl1 = cl1 + 1
          l1map(cl1) = l1
       end if
       cm1 = 0

       do m1 = -l1, l1

          if (m1mask(l1,m1)) then
             cm1 = cm1 + 1
             m1map(l1,cm1) = m1
          end if
          cl2 = 0

          do l2 = 0, lmax2

             if (l2mask(l1,m1,l2)) then
                cl2 = cl2 + 1
                l2map(l1,m1,cl2) = l2
             end if
             cm2 = 0

             do m2 = -l2, l2

                if (m2mask(l1,m1,l2,m2)) then
                   cm2 = cm2 + 1
                   m2map(l1,m1,l2,cm2) = m2
                end if
                cl3 = 0

                do l3 = 0, lmax3

                   if ( l3mask(l1,m1,l2,m2,l3)) then
                      cl3 = cl3 + 1
                      l3map(l1,m1,l2,m2,cl3) = l3
                   end if
                   cm3 = 0

                   do m3 = -l3, l3

                      if ( m3mask(l1,m1,l2,m2,l3,m3) ) then
                         cm3 = cm3 + 1
                         m3map(l1,m1,l2,m2,l3,cm3) = m3
                      end if

                   end do ! cm3
                end do ! cl3
             end do ! cm2
          end do ! cl2
       end do ! cm1
    end do ! cl1

    deallocate(l1mask,m1mask,l2mask,m2mask,l3mask,m3mask)
    deallocate(idxlm)

  end subroutine findgntn0


  subroutine findgntn0_clear
    implicit none
    l1shape = 0
    if (allocated(m1shape)) deallocate(m1shape)
    if (allocated(l2shape)) deallocate(l2shape)
    if (allocated(m2shape)) deallocate(m2shape)
    if (allocated(l3shape)) deallocate(l3shape)
    if (allocated(m3shape)) deallocate(m3shape)
    if (allocated(l1map)) deallocate(l1map)
    if (allocated(m1map)) deallocate(m1map)
    if (allocated(l2map)) deallocate(l2map)
    if (allocated(m2map)) deallocate(m2map)
    if (allocated(l3map)) deallocate(l3map)
    if (allocated(m3map)) deallocate(m3map)
  end subroutine findgntn0_clear


end module m_findgntn0
