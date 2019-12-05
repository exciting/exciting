
subroutine init_kqpoint_set()

    use modinput
    use modmain ,   only : gkmax, bvec, intgv, nspnfv
    use modgw,      only : kset, Gset, Gkset, kqset, Gkqset, Gqset, Gqbarc, fgw
    use mod_kpointset
    use mod_frequency
    use modmpi,     only : rank
    use m_getunit
    implicit none

    ! local variables
    integer :: fid, ig, ik, igk, ispn
    type(k_set) :: ksetnr
    real(8) :: gqmax, gqmaxbarc, v(3), t

    ! call delete_init1_variables()

    !====================
    ! k-points (reduced)
    !====================
    call generate_k_vectors(kset, &
    &                       bvec, &
    &                       input%gw%ngridq, &
    &                       input%gw%vqloff, &
    &                       input%groundstate%reducek)

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
    ! G+k-points (reduced)
    !==========================
    call generate_Gk_vectors(Gkset, &
    &                        kset, &
    &                        Gset, &
    &                        gkmax)

    !==========================
    ! G+kq-points (non-reduced)
    !==========================
    call generate_Gk_vectors(Gkqset, &
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
    ! if ((input%gw%debug).and.(rank==0)) call print_Gk_vectors(Gqset,1,fid)

    !=============================================================
    ! PW basis set used for calculating the bare Coulomb potential
    !=============================================================
    gqmaxbarc = min(input%gw%BareCoul%pwm*gqmax, &
                    input%groundstate%gmaxvr)
    if ((input%gw%debug).and.(rank==0)) then
      write(fid,*)
      write(fid,*) 'Plane-wave cutoff for the bare Coulomb potential &
      &<gqmaxbarc>: ', gqmaxbarc
    end if


    call generate_Gk_vectors(Gqbarc, &
    &                        ksetnr, &
    &                        Gset, &
    &                        gqmaxbarc)
    ! if ((input%gw%debug).and.(rank==0)) call print_Gk_vectors(Gqbarc,1,fid)

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
