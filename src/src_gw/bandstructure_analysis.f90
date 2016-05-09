subroutine bandstructure_analysis(title,ib,nb,nkpt,eband,efermi)

    use modgw,     only : kset, fgw, hev
    use mod_bands, only : nomax, numin, ikvbm, ikcbm, ikvcm, metallic
    use modmpi,    only : rank

    implicit none
    character(len=*), intent(in) :: title
    integer, intent(in) :: ib, nb, nkpt
    real(8), intent(in) :: eband(ib:nb,nkpt)
    real(8), intent(in) :: efermi
    ! local variables
    real(8) :: ebmax, ebmin, egf, ego
    real(8) :: fermidos
    
    real(8), external :: dostet_exciting
    
    if (rank==0) then
      call boxmsg(fgw,'-',trim(title))
      write(fgw,'(a,f10.4)') " Fermi energy: ", efermi
    end if
    
    !-------------------------------------------------------------------
    ! check Fermi energy for correspondence to the specified band range 
    !-------------------------------------------------------------------
    ebmin = minval(eband)
    ebmax = maxval(eband)
    if ((ebmax < efermi) .or. (ebmin > efermi)) then 
        write(*,*) "ERROR(bandstructure_analysis): Fermi energy is outside the specified electronic bands energy range!"
        stop
    end if
    if (rank==0) then
      write(fgw,'(a,2f10.4)') " Energy range: ", ebmin, ebmax
    end if

    !----------------------------------------
    ! Search for the indices of VBM and CBM 
    !----------------------------------------
    call find_vbm_cbm(ib,nb,nkpt,eband,efermi,nomax,numin,ikvbm,ikcbm,ikvcm)
    if (rank==0) then
      write(fgw,'(a,i4)') " Band index of VBM:", nomax
      write(fgw,'(a,i4)') " Band index of CBM:", numin
      write(fgw,*)
    end if
    
    ! Calculate DOS at the fermi level
    fermidos = dostet_exciting(nb-ib+1,nkpt,eband, &
    &                          kset%ntet,kset%tnodes,kset%wtet,kset%tvol, &
    &                          efermi)
    
    ! check for VBM and CBM overlap (metal)
    if ((nomax >= numin).or.(dabs(fermidos)>1.d-4)) then
        metallic = .true.
    else
        metallic = .false.
    end if
    
    if (rank==0) then
      if (metallic) then
        write(fgw,'(a,f8.4)') " DOS at Fermi level: ", fermidos
        write(fgw,*) "WARNING(bandstructure_analysis): Valence and Conduction bands overlap (metal)!"
      else
        egf = eband(numin,ikcbm)-eband(nomax,ikvbm)
        !ego = eband(numin,ikvcm)-eband(nomax,ikvcm)
        ego = eband(numin,1)-eband(nomax,1)
        if (ikvbm == ikcbm) then 
          ! direct gap
          write(fgw,10) ' Direct BandGap (eV):', egf*hev
          write(fgw,11) kset%vkl(:,ikvbm), ikvbm
          write(fgw,10) ' Optical BandGap (eV):', ego*hev
        else
          ! indirect gap
          write(fgw,10) ' Fundamental BandGap (eV):', egf*hev
          write(fgw,12) kset%vkl(:,ikvbm), ikvbm, kset%vkl(:,ikcbm), ikcbm
          write(fgw,10) ' Optical BandGap (eV):', ego*hev
          !write(fgw,11) kset%vkl(:,ikvcm), ikvcm
        end if
      end if
      call linmsg(fgw,'-','')
      call flushifc(fgw)
    end if
    
    10 format(a,T40,f10.4)
    11 format(' at k      = ',3f8.3,' ik = ',i5)
    12 format(' at k(VBM) = ',3f8.3,' ik = ',i5,/,&
    &         '    k(CBM) = ',3f8.3,' ik = ',i5)
    
    return
end subroutine
