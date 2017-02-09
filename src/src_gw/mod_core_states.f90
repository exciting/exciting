module mod_core_states

! maximum number of core states per atom
    integer(4) :: ncmax

! Max. num of core states including lm
    integer(4) :: nclm
    
! Max. L of core states
    integer(4) :: lcoremax
    
! Total num of core states over all atoms including lm
    integer(4) :: ncg
    
! number of core states of each species
    integer(4), allocatable :: ncore(:)
    
! indexes of the core states
    integer(4), allocatable :: corind(:,:)   
!   corind(:,1) = species
!   corind(:,2) = atom of species
!   corind(:,3) = core state
!   corind(:,4) = l
!   corind(:,5) = m

! Modified core radial wave-functions
    real(8), allocatable :: ucore(:,:,:,:)
       
contains

!-------------------------------------------------------------------------------

    subroutine init_core_states()
        use modmain
        implicit none
        integer :: is, ia, ias, ir, l, m, ic, ist
        real(8) :: norm
      
! determine the number of core states for each species (auxiliary arrays)
        if (allocated(ncore)) deallocate(ncore)
        allocate(ncore(nspecies))
        ncmax = 0
        nclm = 0
        ncg = 0
        lcoremax = 0
        do is = 1, nspecies
            ncore(is) = 0
            ic = 0
            do ist = 1, spnst(is)
                if (spcore(ist,is)) then
                    ncore(is) = ncore(is)+1
                    l = spl(ist,is)
                    lcoremax = max(lcoremax,l)
                    do m = -spk(ist,is), spk(ist,is)-1
                        ic = ic+1
                    end do
                end if
            end do
            ncmax = max(ncmax,ncore(is))
            nclm = max(nclm,ic)
            ncg = ncg+ic*natoms(is)
        end do

! setting a unique index for all the core states of all atoms
        if (allocated(corind)) deallocate(corind)
        allocate(corind(ncg,5))
        ic = 0
        do is = 1, nspecies
            do ia = 1, natoms(is)
                ias = idxas(ia,is)
                do ist = 1, ncore(is)
                    l = spl(ist,is)
                    do m = -l, l
                        ic = ic+1
                        corind(ic,1) = is
                        corind(ic,2) = ia
                        corind(ic,3) = ist
                        corind(ic,4) = l
                        corind(ic,5) = m
                    enddo
                enddo
            enddo
        enddo  
        ncg = ic
        
        if (allocated(ucore)) deallocate(ucore)
        allocate(ucore(nrmtmax,2,spnstmax,natmtot))
        ucore(:,:,:,:) = 0.d0

!       According to the definition of core wafefunction in FHIgap code [Eq.(1.1.3)],
!       one has to include the following prefactor into radial part.
     
        do is = 1, nspecies
            do ia = 1, natoms(is)
                ias = idxas(ia,is)
                do ist = 1, ncore(is)
                    l = spl(ist,is)
                    norm = sqrt(0.5d0*spocc(ist,is)/dble(2*l+1))
                    do ir = 1, nrmt(is)
                        ucore(ir,1,ist,ias) = norm*rwfcr(ir,1,ist,ias)/spr(ir,is)
                    end do
                enddo ! ist
            end do ! ia 
        end do ! is
        
    end subroutine
    
!-------------------------------------------------------------------------------

    subroutine delete_core_states()
        if (allocated(ncore)) deallocate(ncore)
        if (allocated(corind)) deallocate(corind)
        if (allocated(ucore)) deallocate(ucore)
    end subroutine
      
end module
