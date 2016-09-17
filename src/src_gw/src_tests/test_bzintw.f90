
subroutine test_bzintw

    use modinput
    use modmain
    use modgw, only   : freq, nomax, numin 
    use mod_kpointset
    use mod_bzintw
    implicit none

    type(k_set)     :: kset
    type(kq_set)    :: kqset
    type(k_bzintw)  :: kbzw
    type(kq_bzintw) :: kqbzw
    integer :: iq

    !----------------
    ! Initialization
    !----------------
    call init_gw

    !--------------------------------------------!
    ! Generate k-point set
    !--------------------------------------------!
    call generate_k_vectors(kset, &
    &                       bvec, &
    &                       input%gw%ngridq, &
    &                       input%gw%vqloff, &
    &                       input%gw%reduceq)
    call print_k_vectors(kset,6)
    
    !--------------------------------------------!
    ! Generate k-weights
    !--------------------------------------------!
    call generate_k_bzintw(kbzw,kset,.true.)
    call print_k_bzintw(kbzw,6)
    call delete_k_bzintw(kbzw)
    
    !--------------------------------------------!
    ! Generate k/q-vectors
    !--------------------------------------------!
    call generate_kq_vectors(kqset, &
    &                        bvec, &
    &                        input%gw%ngridq, &
    &                        input%gw%vqloff, &
    &                        input%gw%reduceq)
    call print_kq_vectors(kqset,6)
    
    !--------------------------------------------!
    ! Generate q-dependent k-weights
    !--------------------------------------------!
    iq = 1
    call generate_kq_bzintw(kqbzw,iq,kset,kqset,freq,nomax,numin,.true.)
    call print_kq_bzintw(kqbzw,iq,6)
    call delete_kq_bzintw(kqbzw)
    
    call delete_k_vectors(kset)
    call delete_kq_vectors(kqset)
    
    return
end subroutine
