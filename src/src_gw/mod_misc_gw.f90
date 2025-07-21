!> miscellaneous variables used in GW part just for convenience  
module mod_misc_gw
    use modinput
    use modmain
    use constants, only: maxatoms, maxspecies, pi
    use precision, only: dp
    
    implicit none

! shortcut for atomic position array
    real(8) :: atposl(3, maxatoms, maxspecies)
! lengths of the basis vectors
    real(8) :: alat(3)
! 2*pi/a, 2*pi/b, 2*pi/c
    real(8) :: pia(3)
! 1/omega - reciprocal unitcell volume
    real(8) :: vi
! spin-polarized flag
    logical :: spinpol
! Ha --> eV
    real(8), parameter :: hev=27.21138505d0
! Check if the point is \Gamma
    logical :: Gamma

    ! Characteristic function in real space
    complex(8), allocatable :: zfunir(:)

contains

!-------------------------------------------------------------------------------

    subroutine init_misc_gw
        use modmain, only: avec
        use m_zfftifc, only: zfftifc
        implicit none
        integer :: i, is, ia, ias
        integer :: ig, ifg

! reciprocal cell volume
        vi = 1.0d0/omega

! reciprocal lattice basis lengths
        do i = 1, 3
            alat(i) = dsqrt(avec(1,i)*avec(1,i)+ &
           &                avec(2,i)*avec(2,i)+ &
           &                avec(3,i)*avec(3,i))
            pia(i) = 2.0d0*pi/alat(i)
        end do

! additional arrays
        do is = 1, nspecies
            do ia = 1, natoms(is)
                ias = idxas(ia,is)
! shortcut for atomic positions
                atposl(:,ia,is) = &
               &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            end do
        end do

        ! FFT of the characteristic function
        if (allocated(zfunir)) deallocate(zfunir)
        allocate(zfunir(ngrtot))
        zfunir(:) = zzero
        do ig = 1, ngrtot
          ifg = igfft(ig)
          zfunir(ifg) = conjg(cfunig(ig))
        end do
        call zfftifc(3,ngrid,1,zfunir)

        return
     end subroutine

!-------------------------------------------------------------------------------

     !> Check if all coordinates of `vec` are close to zero.
     !> If `vec` is a k-point, that means it can be interpreted as the \(\Gamma\)-point
     pure logical function gammapoint(vec, tol)
        real(dp), intent(in) :: vec(3)
        real(dp), intent(in), optional :: tol

        real(dp), parameter :: default_tol = 1.e-6_dp
        real(dp) :: tolerance 

        tolerance = default_tol
        if( present(tol) ) tolerance = tol
        gammapoint = ( norm2(vec) <= tolerance )
     end function gammapoint


end module
