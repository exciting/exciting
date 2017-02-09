      
subroutine task_chi0_q

    use modinput
    use modmain,               only : zzero, evalsv, efermi
    use modgw
    use mod_mpi_gw
    use m_getunit
    use mod_hdf5
    use mod_rpath
            
    implicit none
    integer(4) :: ikp, iq, fid, ik, im, iom, nempty
    real(8)    :: t0, t1
    integer(4) :: recl
    character(len=64) :: filename
    complex(8) :: trace

    !===========================================================================
    ! Initialization
    !===========================================================================
    
    iq = input%gw%iik
    
    ! prepare GW global data
    call init_gw
    
    ! clean not used anymore global exciting variables
    call clean_gndstate
    
    ! occupancy dependent BZ integration weights
    call kintw
    
    !===========================================================================
    ! Main loop: BZ integration
    !===========================================================================    

    iomstart = 1
    iomend = freq%nomeg

    write(*,*)
    write(*,*) '(task_chi0_q): q-point cycle, iq = ', iq
    
    Gamma = .false.
    
    !========================================
    ! Calculate interstitial basis functions
    !========================================
    matsiz = locmatsiz+Gqset%ngk(1,iq)
    call diagsgi(iq)
    call calcmpwipw(iq)
    
    !======================================
    ! Coulomb potential (not used)
    !======================================
    mbsiz = matsiz
    if (allocated(barc)) deallocate(barc)
    allocate(barc(matsiz,mbsiz))
    barc(:,:) = 0.d0
    do im = 1, matsiz
      barc(im,im) = 1.d0
    end do
          
    !===================================
    ! Calculate the chi0 function
    !===================================
    call init_dielectric_function(mbsiz,iomstart,iomend,Gamma)
    call calcepsilon(iq,iomstart,iomend)
    do im = 1, mbsiz
      epsilon(im,im,iomstart:iomend) = epsilon(im,im,iomstart:iomend)-zone
    end do
    
    if (rank==0) then
      ! Set the name of the output file
      5 format('chi0-q_',i4,'-nempty_',i8,'.out')
      nempty = nstsv
      if (nstsv>input%gw%nempty) nempty = nempty-int(chgval/2+1)
      write(filename,5) iq, nempty
      call str_strip(filename)
      open(unit=77,file=filename,status='unknown')
      write(77,'("# iq=",i8,4x,"nempty=",i8)') iq, nempty
      do iom = iomstart, iomend
        trace = 0.d0 
        do im = 1, mbsiz
          trace = trace+epsilon(im,im,iom)
        end do
        trace = trace/dble(mbsiz)
        write(77,'(i8,2f18.6)') iom, trace
      end do
      close(77)
    end if
    
    call delete_dielectric_function(Gamma)
    if (allocated(kcw)) deallocate(kcw)
    if (allocated(unw)) deallocate(unw)
    
    ! clean unused data
    if (allocated(mpwipw)) deallocate(mpwipw)
    if (allocated(barc)) deallocate(barc)
    
    if (allocated(kiw)) deallocate(kiw)
    if (allocated(ciw)) deallocate(ciw)
    
    if (allocated(evalsv)) deallocate(evalsv)
    call delete_freqgrid(freq)
    call delete_k_vectors(kset)
    call delete_G_vectors(Gset)
    call delete_Gk_vectors(Gkset)
    call delete_kq_vectors(kqset)
    call delete_Gk_vectors(Gqset)
    call delete_Gk_vectors(Gqbarc)
    
    return
end subroutine

