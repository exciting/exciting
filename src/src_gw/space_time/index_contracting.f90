module index_contracting
    use mod_atoms, only: nspecies, natmtot, natoms, idxas
    use mod_product_basis, only: locmatsiz, nmix, bigl
    implicit none
    private 
    public :: nlm_index_contracting, atomic_id_to_species_and_atom

    contains

    subroutine nlm_index_contracting(nlm, nlm_indices)
        !> Total number for combined N, L, M 
        integer, intent(out) :: nlm
        !> Combined N, L, M index
        integer, allocatable, intent(out) :: nlm_indices(:, :, :)
        
        !! Local variables
        integer :: is, ia, ias
        integer :: n, l, m
        integer, allocatable :: nlm_indices_tmp(:, :, :)

        !! With only ias information
        ! allocate(nlm_indices_tmp(locmatsiz, natmtot, 3))
        ! do ias = 1, natmtot
        !   nlm = 0
        !   do n = 1, nmix(ias)
        !     l = bigl(n, ias)
        !     do m = -l, l
        !       nlm = nlm + 1
        !       nlm_indices_tmp(nlm, ias, 1) = n
        !       nlm_indices_tmp(nlm, ias, 2) = l
        !       nlm_indices_tmp(nlm, ias, 3) = m
        !     enddo
        !   enddo
        ! enddo
        
        ! allocate(nlm_indices(nlm, natmtot, 3), source = nlm_indices_tmp(1:nlm, :, :))
        ! deallocate(nlm_indices_tmp)


        !! With is and ia and ias information as well
        allocate(nlm_indices_tmp(locmatsiz, natmtot, 5))
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia,is)
            nlm = 0
            do n = 1, nmix(ias)
              l = bigl(n, ias)
              do m = -l, l
                nlm = nlm + 1
                nlm_indices_tmp(nlm, ias, 1) = is
                nlm_indices_tmp(nlm, ias, 2) = ia
                nlm_indices_tmp(nlm, ias, 3) = n
                nlm_indices_tmp(nlm, ias, 4) = l
                nlm_indices_tmp(nlm, ias, 5) = m
              enddo 
            enddo
          enddo
        enddo
        
        allocate(nlm_indices(nlm, natmtot, 5), source = nlm_indices_tmp(1:nlm, :, :))
        deallocate(nlm_indices_tmp)

    end subroutine nlm_index_contracting

    subroutine atomic_id_to_species_and_atom(nspecies, natoms, natmtot, isia)
    !> Total number of species 
    integer, intent(in) :: nspecies
    !> Number of atoms per each species
    integer, intent(in) :: natoms(:)
    !> Total number of atoms including all species
    integer, intent(in) :: natmtot
    !> Map atomic id (ias) to individual species and atom
    integer, intent(out) :: isia(natmtot, 2)

    !! Local variables 
    integer :: is, ia, ias 

    ias = 0
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = ias + 1
        isia(ias, 1) = is
        isia(ias, 2) = ia
      enddo 
    enddo 

    end subroutine atomic_id_to_species_and_atom

end module index_contracting