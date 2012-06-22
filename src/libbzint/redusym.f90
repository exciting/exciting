
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

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
  !<sag>
  use control, only: tetradbglv, tetraifc
  !</sag>

  ! !INPUT PARAMETERS:

  implicit none

  integer(4), intent(in) :: nsym
  integer(4), intent(in) :: symmat(3,3,nsym)
  integer(4), intent(in) :: dvof

! !OUTPUT PARAMETERS:

  integer(4), intent(out) :: nsymr

! !LOCAL VARIABLES:

  integer(4) :: i
  integer(4) :: isym
  integer(4) :: isymr

  integer(4) :: t1(3)
  integer(4) :: t2(3)
  integer(4) :: tsdif(3)
  integer(4) :: toff(3)
  integer(4) :: tmpsym(3,3,nsym)

  logical :: sgroupsh, tsdink

! !REVISION HISTORY:
!
! Created May 2006 by RGA
!
!EOP
!BOC
  !<sag>
  if (tetraifc.eq."wien2k") then

!!! *** OLD CODE ***
!!$     tmpsym=0
!!$     isymr=0
!!$     do isym=1,nsym
!!$        do i=1,3
!!$           t1(i)=symmat(i,1,isym)*shift(1)+symmat(i,2,isym)*shift(2)+   &
!!$                &          symmat(i,3,isym)*shift(3)
!!$           t2(i)=mod(t1(i),dvof*div(i))
!!$           toff(i)=t2(i)+(1-isign(1,t2(i)))*dvof*div(i)/2
!!$        enddo
!!$        do i=1,3
!!$           tsdif(i)=mod(abs(toff(i)-shift(i)),dvof)
!!$           !          write(24,'(3i4,5i6)')symmat(1:3,i,isym),shift(i),    &
!!$           !     &           t1(i),t2(i),toff(i),tsdif(i)
!!$        enddo
!!$        sgroupsh=((toff(1).eq.shift(1)).and.(toff(2).eq.shift(2)).and.  &
!!$             &       (toff(3).eq.shift(3)))
!!$        tsdink=((tsdif(1).eq.0).and.(tsdif(2).eq.0).and.(tsdif(3).eq.0))
!!$        !        write(24,'(i4,2l4)')isym,sgroupsh,tsdink
!!$        if(tsdink.or.sgroupsh)then
!!$           isymr=isymr+1
!!$           tmpsym(1:3,1:3,isymr)=symmat(1:3,1:3,isym)
!!$        endif
!!$     enddo
!!$     nsymr=isymr
!!$     !<sag>
!!$     if (tetradbglv > 0) then
!!$        !</sag>
!!$        if(nsymr.lt.nsym)write(*,*)'WARNING: The k-point offset ',       &
!!$             & 'selected reduces the symmetry group of the mesh, resulting in '&
!!$             & ,'a larger number of irreducible k-points'
!!$        !<sag>
!!$     end if
!!$     !</sag>
!!$     allocate(iio(3,3,nsymr))
!!$     iio(1:3,1:3,1:nsymr)=tmpsym(1:3,1:3,1:nsymr)

     ! return all symmetry operations and let kgen do the considerations
     ! for the offset of the k-mesh
     nsymr=nsym
     allocate(iio(3,3,nsymr))
     iio(1:3,1:3,1:nsymr)=symmat(1:3,1:3,1:nsymr)

  else if (tetraifc.eq."exciting") then

     ! return all symmetry operations and let kgen do the considerations
     ! for the offset of the k-mesh
     nsymr=nsym
     allocate(iio(3,3,nsymr))
     iio(1:3,1:3,1:nsymr)=symmat(1:3,1:3,1:nsymr)

  end if
  !</sag>
  return

end subroutine redusym
!EOC
