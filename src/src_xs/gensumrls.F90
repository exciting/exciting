
module m_gensumrls
  implicit none
contains

  subroutine gensumrls(w,eps,sumrls)
    !
    ! Expressions for the sumrules taken from  
    ! [cad, CPC 175 (2006) 1-14, p5, eq. 26]
    !
    use modmain
    use modxs
    implicit none
    ! arguments
    real(8), intent(in) :: w(:)
    complex(8), intent(in) :: eps(:)
    real(8), intent(out) :: sumrls(3)
    ! local variables
    character(*), parameter :: thisnam = 'gensumrls'
    real(8), allocatable :: f(:), cf(:,:),g(:), om(:)
    real(8) :: delt
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
!!!/// add analytic stuff ///
    f(:)=aimag(-1/eps(:))/(w(:)+1.d-7)
    call fderiv(-1,n,w,f,g,cf)
    sumrls(3)=g(n)

    deallocate(f,g,cf)

  end subroutine gensumrls

end module m_gensumrls
