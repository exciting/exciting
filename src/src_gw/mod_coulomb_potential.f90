!--------------------------------------------!
!     Bare Coulomb potential related data    !
!--------------------------------------------!

module mod_coulomb_potential
    
    ! The lattice summations matrix      
    complex(8), allocatable :: sgm(:,:,:)
    
    ! The matrix representation of the bare coulomb potential in the mixed basis            
    complex(8), allocatable :: barc(:,:)

    ! full set of the eigenvalues of barcoul matrix
    real(8), allocatable :: barcev(:)
      
    ! full set of eigenvectors of barcoul matrix        
    complex(8), allocatable :: vmat(:,:)
    
    ! Matrix elements between mixed functions and constant function
    complex(8), allocatable :: wi0(:)
    
    ! use geometry dependent cutoff for the Coulomb potential 
    logical :: vccut
    
    ! spherical integral over the Coulomb singularity
    real(8) :: rccut, i_sz
    
contains
 
    subroutine delete_coulomb_potential()
      implicit none
      if (allocated(vmat)) deallocate(vmat)
      if (allocated(barcev)) deallocate(barcev)
    end subroutine
    
end module
