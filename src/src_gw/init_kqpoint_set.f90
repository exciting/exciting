
subroutine init_kqpoint_set()

    use modinput
    use modmain ,   only : gkmax, bvec, intgv, nspnfv
    use modgw,      only : kset, Gset, Gkset, kqset, Gqset, Gqbarc, fgw
    use mod_kpointset
    use mod_frequency
    use modmpi,     only : rank
    use m_getunit
    implicit none
    
    ! local variables
    integer :: fid, ig, ik, igk, ispn
    type(k_set) :: ksetnr
    real(8) :: gqmax, gqmaxbarc, v(3), t
    
    !call delete_init1_variables
    
    !====================
    ! k-points (reduced)
    !====================
    call generate_k_vectors(kset, &
    &                       bvec, &
    &                       input%gw%ngridq, &
    &                       input%gw%vqloff, &
    &                       input%groundstate%reducek)
    if (rank==0) call print_k_vectors(kset,fgw)
    
    if ((input%gw%debug).and.(rank==0)) then
      call getunit(fid)
      open(fid,file='GW_KPOINTS.OUT',action='Write',status='Unknown')
      call print_k_vectors(kset,fid)
    end if
    
    !========================
    ! k-points (non-reduced)
    !========================
    call generate_k_vectors(ksetnr, &
    &                       bvec, &
    &                       input%gw%ngridq, &
    &                       input%gw%vqloff, &
    &                       .false.)
    
    !==========
    ! G-points
    !==========
    call generate_G_vectors(Gset, &
    &                       bvec, &
    &                       intgv, &
    &                       input%groundstate%gmaxvr)
    if ((input%gw%debug).and.(rank==0)) call print_G_vectors(Gset,fid)
    
    !==========================
    ! G+k-points (non-reduced)
    !==========================
    call generate_Gk_vectors(Gkset, &
    &                        ksetnr, &
    &                        Gset, &
    &                        gkmax)
    
    !============
    ! k/q-points
    !============
    call generate_kq_vectors(kqset, &
    &                        bvec, &
    &                        input%gw%ngridq, &
    &                        input%gw%vqloff, &
    &                        input%gw%reduceq)
    if ((input%gw%debug).and.(rank==0)) call print_kq_vectors(kqset,fid)
    
    !==========================
    ! G+q-points (non-reduced)
    !==========================
    gqmax = gkmax*input%gw%mixbasis%gmb
    if ((input%gw%debug).and.(rank==0)) then
      write(fid,*)
      write(fid,*) 'Plane-wave cutoff for G+q <gqmax>: ', gqmax
      write(fid,*)
    end if
    
    call generate_Gk_vectors(Gqset, &
    &                        ksetnr, &
    &                        Gset, &
    &                        gqmax)
    !if ((input%gw%debug).and.(rank==0)) call print_Gk_vectors(Gqset,fid)
    
    !=============================================================
    ! PW basis set used for calculating the bare Coulomb potential
    !=============================================================
    gqmaxbarc = min(input%gw%BareCoul%pwm*gqmax, &
    &              input%groundstate%gmaxvr)
    if ((input%gw%debug).and.(rank==0)) then
      write(fid,*)
      write(fid,*) 'Plane-wave cutoff for the bare Coulomb potential &
      &<gqmaxbarc>: ', gqmaxbarc
    end if
    
    ! determine the number of G+q combinations which satisfy |G+q|<gqmaxbarc
    if (allocated(Gqbarc%ngk)) deallocate(Gqbarc%ngk)
    allocate(Gqbarc%ngk(nspnfv,kqset%nkpt))
    if (allocated(Gqbarc%igkig)) deallocate(Gqbarc%igkig)
    allocate(Gqbarc%igkig(Gset%ngrtot,nspnfv,kqset%nkpt))
    Gqbarc%igkig(:,:,:) = 0
    if (allocated(Gqbarc%igigk)) deallocate(Gqbarc%igigk)
    allocate(Gqbarc%igigk(Gset%ngrtot,nspnfv,kqset%nkpt))
    Gqbarc%igigk(:,:,:) = 0
    do ispn = 1, nspnfv
      do ik = 1, kqset%nkpt
        igk = 0
        do ig = 1, Gset%ngrtot
          v(:) = Gset%vgc(:,ig)+kqset%vkc(:,ik)
          t = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
          if (t < gqmaxbarc) then
            igk = igk+1
            Gqbarc%igigk(ig,ispn,ik) = igk
            Gqbarc%igkig(igk,ispn,ik) = ig
          end if
        end do ! ig
        Gqbarc%ngk(ispn,ik) = igk
      end do ! ik
    end do ! ispn

    ! maximum number of G+k vectors
    Gqbarc%ngkmax = maxval(Gqbarc%ngk)
    if (Gqbarc%ngkmax > Gset%ngrtot) then
      write(*,*) 'ERROR(init_kqpoint_set) ngkmax > ngrtot'
      stop
    end if
      
    !============================    
    ! delete the dummy k-set
    call delete_k_vectors(ksetnr)
    ! close the k-point info file
    if ((input%gw%debug).and.(rank==0)) close(fid)

contains

!===============================================================================
    subroutine delete_init1_variables()
        use modmain
        
        if (allocated(ik2ikp))  deallocate(ik2ikp)
        if (allocated(ikp2ik))  deallocate(ikp2ik)
        if (allocated(iwkp))    deallocate(iwkp)
        if (allocated(wtet))    deallocate(wtet)
        if (allocated(tnodes))  deallocate(tnodes)
        
        if (allocated(ivknr))   deallocate (ivknr)
        if (allocated(vklnr))   deallocate (vklnr)
        if (allocated(vkcnr))   deallocate (vkcnr)
        if (allocated(wkptnr))  deallocate (wkptnr)
        if (allocated(ikmapnr)) deallocate (ikmapnr)
        
        if (allocated(ngk))     deallocate (ngk)
        if (allocated(igkig))   deallocate (igkig)
        if (allocated(vgkl))    deallocate (vgkl)
        if (allocated(vgkc))    deallocate (vgkc)
        if (allocated(gkc))     deallocate (gkc)
        if (allocated(tpgkc))   deallocate (tpgkc)
        if (allocated(sfacgk))  deallocate (sfacgk)
    
    end subroutine
    
end subroutine
