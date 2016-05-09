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
      do isym=1,nsym
        do i=1,3
          t1(i)=symmat(i,1,isym)*shift(1)+symmat(i,2,isym)*shift(2)+   &
     &          symmat(i,3,isym)*shift(3)
          t2(i)=mod(t1(i),dvof*div(i))
          toff(i)=t2(i)+(1-isign(1,t2(i)))*dvof*div(i)/2
        enddo
        do i=1,3
          tsdif(i)=mod(abs(toff(i)-shift(i)),dvof)
        enddo
        sgroupsh=((toff(1).eq.shift(1)).and.(toff(2).eq.shift(2)).and.  &
     &       (toff(3).eq.shift(3)))
        tsdink=((tsdif(1).eq.0).and.(tsdif(2).eq.0).and.(tsdif(3).eq.0))
        if(tsdink.or.sgroupsh)then
          isymr=isymr+1
          tmpsym(1:3,1:3,isymr)=symmat(1:3,1:3,isym)
        endif
      enddo
      nsymr=isymr
      if(nsymr.lt.nsym) write(6,*)'WARNING in libbzint%redusym :', &
     & '  The k-point offset selected reduces the symmetry group of',&
     & '  the mesh, resulting in a larger number of irred. k-points'

      allocate(iio(3,3,nsymr))
      iio(1:3,1:3,1:nsymr)=tmpsym(1:3,1:3,1:nsymr)

      return

      end subroutine redusym
!EOC
