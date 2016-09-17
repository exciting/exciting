!-------------------------------------------------------------------!
!     miscellaneous variables used in GW part just for convenience  !
!-------------------------------------------------------------------!

#include "maxdefinitions.inc"

module mod_misc_gw

    use modinput
    use modmain
    implicit none

! shortcut for atomic position array
    real(8) :: atposl(3,_MAXATOMS_,_MAXSPECIES_)
! shortcut for basis vectors
    real(8) :: avec(3,3)
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
        implicit none
        integer :: i, is, ia, ias
        integer :: ig, ifg
    
! reciprocal cell volume
        vi = 1.0d0/omega

! shortcut for basis vectors 
        avec(:,1) = input%structure%crystal%basevect(:,1)
        avec(:,2) = input%structure%crystal%basevect(:,2)
        avec(:,3) = input%structure%crystal%basevect(:,3)

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

     logical function gammapoint(vec)
        implicit none
        real(8), intent(in) :: vec(3)
        real(8) :: len
        len = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
        if (len > 1.0d-6) then
          gammapoint = .false.
        else
          gammapoint = .true.
        endif
        return
     end function gammapoint


end module
