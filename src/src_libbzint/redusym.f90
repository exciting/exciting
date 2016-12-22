!BOP
!
! !ROUTINE: redusym
!
! !INTERFACE:
      subroutine redusym(nsym,symmat,dvof,nsymr)

! !DESCRIPTION:
!
! This subroutine reduces the symmetry operations in the BZ according to
! the offset. 
!
! !USES:

      use kgen_internals

! !INPUT PARAMETERS:

      implicit none

      integer, intent(in) :: nsym
      integer, intent(in) :: symmat(3,3,nsym)
      integer, intent(in) :: dvof

! !OUTPUT PARAMETERS:

      integer, intent(out) :: nsymr

! !LOCAL VARIABLES:

      integer :: i
      integer :: isym
      integer :: isymr

      integer :: t1(3)
      integer :: t2(3)
      integer :: tsdif(3)
      integer :: toff(3)
      integer :: tmpsym(3,3,nsym)
 
      logical :: sgroupsh, tsdink

! !REVISION HISTORY:
!
! Created May 2006 by RGA
!
!EOP
!BOC
      tmpsym=0
      isymr=0

      do isym = 1, nsym

        ! Apply symmetry operation number isym
        ! to offset vector in a scaled k-grid.
        ! Note:
        ! vkloff = shift/dvof
        ! div = ngridk
        ! So vkloff lies on a grid point 
        ! of the scaled ngridk with scaling factor dvof
        do i = 1, 3

          ! Apply symmetry operation on the offset in the scaled k grid
          ! S*shift = t1
          t1(i) = symmat(i,1,isym)*shift(1)&
               &+ symmat(i,2,isym)*shift(2)&
               &+ symmat(i,3,isym)*shift(3)

          ! Map resulting vector back into scaled k grid (why not modulo??)
          t2(i) = mod(t1(i),dvof*div(i))
          ! Add lattive vector in k-ponint coordinates to the elements of t2
          ! which are negative 
          toff(i) = t2(i) + (1-isign(1,t2(i)))/2*dvof*div(i)

        end do

        ! Get part of the difference between original offset in
        ! scaled coordinates and symmetry equivalent offset that
        ! lies on k points not present in the unscaled k-grid
        do i = 1, 3
          tsdif(i) = mod(abs(toff(i)-shift(i)), dvof)
        end do

        ! Check if symmetry operation leaves offset unchanged.
        sgroupsh = all(toff == shift)

        ! Check if symmetry operation maps the offset to a point
        ! present in the original coarse grid.
        tsdink = all(tsdif == 0)

        ! Save the preserved symmetry operation
        if(tsdink .or. sgroupsh) then
          isymr=isymr+1
          tmpsym(1:3,1:3,isymr)=symmat(1:3,1:3,isym)
        end if

      end do

      nsymr = isymr

      if(nsymr .lt. nsym) then
        write(6,*)'WARNING in libbzint%redusym :', &
        & '  The k-point offset selected reduces the symmetry group of',&
        & '  the mesh, resulting in a larger number of irred. k-points'
      end if

      ! Save valid symmetry operations to module variable.
      allocate(iio(3,3,nsymr))
      iio(1:3,1:3,1:nsymr)=tmpsym(1:3,1:3,1:nsymr)

      return

      end subroutine redusym
!EOC
