subroutine wannier_bessel_int( lmax, barg, rpt, nrpt, rmax, f1, f2, bint)
    use modmain
    implicit none
! arguments
    real(8), intent( in) :: f1( nrpt), f2( nrpt)
    integer, intent( in) :: nrpt, lmax, rmax
    real(8), intent( in) :: rpt( nrpt), barg
    real(8), intent( out) :: bint( 0:lmax)
! local variables
    integer :: ir, l
    real(8) :: x
! allocatable arrays
    real(8), allocatable :: jlkr(:,:)
    real(8), allocatable :: uf(:), gr(:), cf(:,:)
    
! compute spherical bessel functions
    allocate( jlkr( lmax, nrpt))
    do ir = 1, nrpt
      x = barg*rpt( ir)
      Call sbessel( lmax, x, jlkr( :, ir))
    end do
! do integration
    allocate( uf( nrpt), gr( nrpt), cf( 3, nrpt))
    do l = 0, lmax
      do ir = 1, nrpt
        uf( ir) = f1( ir)*jlkr( l, ir)*f2( ir)*rpt( ir)**2
      end do
      call fderiv( -1, nrpt, rpt, uf, gr, cf)
      bint( l) = gr( rmax)
    end do
    return
end subroutine
