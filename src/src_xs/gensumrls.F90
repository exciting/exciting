
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_gensumrls
  implicit none
contains

!BOP
! !ROUTINE: gensumrls
! !INTERFACE:
  subroutine gensumrls(w,eps,sumrls)
! !USES:
    use modmain
    use modxs
! !INPUT/OUTPUT PARAMETERS:
!   w      : frequency grid (in,real(:))
!   eps    : dielectric function tensor component (in,complex(:))
!   sumrls : values of three different sumrules (out,real(3))
! !DESCRIPTION:
!   Expressions for the sumrules taken from Ambrosch-Draxl,
!   CPC 175 (2006) 1-14, p5, Eq. 26.
!
! !REVISION HISTORY:
!   Created March 2006 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    real(8), intent(in) :: w(:)
    complex(8), intent(in) :: eps(:)
    real(8), intent(out) :: sumrls(3)
    ! local variables
    character(*), parameter :: thisnam = 'gensumrls'
    real(8), allocatable :: f(:), cf(:,:),g(:)
    integer :: n1(1),n

    if (any(shape(w).ne.shape(eps))) then
       write(unitout,'(a)') 'Error('//thisnam//'): input arrays have &
            &diffenrent shape'
       call terminate
    end if

    n1=shape(w)
    n=n1(1)
    allocate(f(n),g(n),cf(3,n))

    ! zeroth frequency moment sumrule
    f(:)=aimag(eps(:))
    call fderiv(-1,n,w,f,g,cf)
    sumrls(1)=g(n)

    ! first frequency moment sumrule
    f(:)=aimag(-1/eps(:))*w(:)
    call fderiv(-1,n,w,f,g,cf)
    sumrls(2)=g(n)

    ! one over frequency sumrule (pi half sumrule)
    f(1)=0.d0
    if (n.gt.1) f(2:)=aimag(-1/eps(2:))/w(2:)
    call fderiv(-1,n,w,f,g,cf)
    sumrls(3)=g(n)

    deallocate(f,g,cf)

  end subroutine gensumrls
!EOC

end module m_gensumrls
