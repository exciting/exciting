
subroutine test_kqpts
    
    use modinput
    use modmain
    use mod_kpointset
    
    implicit none

    type(k_set)  :: kset, ksetnr
    type(kq_set) :: kqset
    type(G_set)  :: Gset
    type(Gk_set) :: Gkset
    
    call init0
    
!--------------------------------------------!
!   Generate reduced k-point set
!--------------------------------------------!
    call generate_k_vectors(kset, &
    &                       bvec, &
    &                       input%gw%ngridq, &
    &                       input%gw%vqloff, &
    &                       .true.)
    call print_k_vectors(kset,6)
    call delete_k_vectors(kset)
    
!--------------------------------------------!
!   Generate non-reduced k-point set
!--------------------------------------------!
    call generate_k_vectors(ksetnr, &
    &                       bvec, &
    &                       input%gw%ngridq, &
    &                       input%gw%vqloff, &
    &                       .false.)
    call print_k_vectors(ksetnr,6)
    call delete_k_vectors(ksetnr)
    
!--------------------------------------------!
!   Generate G-space
!--------------------------------------------!
    write(*,*) 'gmaxvr=', input%groundstate%gmaxvr
    call generate_G_vectors(Gset, &
    &                       bvec, &
    &                       intgv, &
    &                       input%groundstate%gmaxvr)
    call print_G_vectors(Gset,6)
    call delete_G_vectors(Gset)
    
!--------------------------------------------!
!   Generate G+k-vectors
!--------------------------------------------!
    call generate_k_vectors(kset, &
    &                       bvec, &
    &                       input%gw%ngridq, &
    &                       input%gw%vqloff, &
    &                       .true.)
    call print_k_vectors(kset,6)
    
    write(*,*) 'gmaxvr=', input%groundstate%gmaxvr
    call generate_G_vectors(Gset, &
    &                       bvec, &
    &                       intgv, &
    &                       input%groundstate%gmaxvr)
    call print_G_vectors(Gset,6)
    
    gkmax = input%groundstate%rgkmax / rmt(input%groundstate%isgkmax)
    write(*,*) 'gkmax=', gkmax
    
    call generate_Gk_vectors(Gkset,kset,Gset,gkmax)
    call print_Gk_vectors(Gkset,1,6)
    
    call delete_k_vectors(kset)
    call delete_G_vectors(Gset)
    call delete_Gk_vectors(Gkset)
    
!--------------------------------------------!
!   Generate k/q-vectors
!--------------------------------------------!
    call generate_kq_vectors(kqset, &
    &                        bvec, &
    &                        input%gw%ngridq, &
    &                        input%gw%vqloff, &
    &                        .true.)
    call print_kq_vectors(kqset,6)
    call delete_kq_vectors(kqset)
    
    return
end subroutine

