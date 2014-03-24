subroutine bandanalysis(title,ib,nb,eband,efermi)

    use modmain, only: nkpt, vkl
    use modgw,   only: nomax, numin, fgw, hev

    implicit none
    character(len=*), intent(in) :: title
    integer, intent(in) :: ib, nb
    real(8), intent(in) :: eband(ib:nb,nkpt), efermi

    integer :: i,ik,nc,nv,ikv,ikc
    integer :: ikvm,ikcm
    real(8) :: kvec(3),w,kvec0(3)
    real(8) :: evbm,ecbm,egap(3),ebmin,ebmax
    logical :: lmetal

    call boxmsg(fgw,'-',trim(title)//" Band Analysis")
    write(fgw,'(a,2i5)'  ) "  Range of bands considered: ", ib, nb
    write(fgw,'(a,f10.4,a)') "  Fermi Energy = ", efermi
!
!   Find Valence Band Maximum (VBM) and Cond. Band Minimum (CBM)
!
    ebmin = minval(eband)
    ebmax = maxval(eband)
    if ((ebmax.lt.efermi) .or. (ebmin.gt.efermi)) then 
        write(fgw,*) "WARNING(bandanaly): "
        write(fgw,*) " - Fermi energy outside the banstructure energy range!" 
        write(fgw,'(a,f10.4)') " - Fermi Energy  :     ", efermi
        write(fgw,'(a,f10.4)') " - Maximal energy:     ", ebmax
        write(fgw,'(a,f10.4)') " - Minimal energy:     ", ebmin
        return
    endif 

    nomax = 0
    numin = 1000
    ikvm = 1
    ikcm = 1
    lmetal = .false.
    do ik = 1, nkpt
        nv = 0
        nc = 1000
        do i = ib, nb
            if (eband(i,ik).lt.efermi) then
                if (i.gt.nv) nv = i
            else
                if (i.lt.nc) nc = i
            endif
        enddo
        if (nv.gt.nomax) nomax = nv 
        if (nc.lt.numin) numin = nc 
        nv = nomax
        nc = numin 
        ikv = ikvm
        ikc = ikcm
        if (eband(nv,ik).gt.eband(nv,ikv)) ikvm = ik
        if (eband(nc,ik).lt.eband(nc,ikc)) ikcm = ik
    enddo
      
    if (nomax.lt.ib) then
        write(fgw,*) "ERROR(bandanaly): nomax is out of the specified band ranges !!!"
        write(fgw,*) "nomax=", nomax, " < ib=", ib
        stop
    endif
    if (nomax.gt.nb) then
        write(fgw,*) "ERROR(bandanaly): nomax is out of the specified band ranges !!!"
        write(fgw,*) "nomax=", nomax, " > nb=", nb
        stop
    endif
    if (numin.lt.ib) then
        write(fgw,*) "ERROR(bandanaly): numin is out of the specified band ranges !!!"
        write(fgw,*) "numin=", numin, " < ib=", ib
        stop
    endif
    if (numin.gt.nb) then
        write(fgw,*) "ERROR(bandanaly): numin is out of the specified band ranges !!!"
        write(fgw,*) "numin=", numin, " > nb=", nb
        stop
    endif
    
    write(fgw,111) '  Band index for VBM and CBM=', nv, nc
    
    if (nomax >= numin) lmetal=.true.
    if (lmetal) then
      write(fgw,*) "WARNING(bandanaly): Valence and Conductance bands overlap (metal)!"
    else
      egap(1)= eband(numin,ikcm)-eband(nomax,ikvm)
      egap(2)= minval(eband(numin:nb,ikvm))-eband(nomax,ikvm)
      egap(3)= eband(numin,ikcm)-maxval(eband(ib:nomax,ikcm))
      if(ikvm.eq.ikcm) then ! direct gap 
        write(fgw,112) '[eV]', egap(1)*hev
        write(fgw,113) vkl(:,ikvm), ikvm
      else 
        write(fgw,114) '[eV]', egap(1:3)*hev
        write(fgw,115) vkl(:,ikvm), ikvm, vkl(:,ikcm),ikcm
      end if
      write(fgw,*) ' Range of each band [eV]:'
      write(fgw,'(a5,3a12)') 'n','Bottom','Top','Width'
      do i = ib, nb
        ebmin = minval(eband(i,1:nkpt))*hev
        ebmax = maxval(eband(i,1:nkpt))*hev
        write(fgw,'(i5,3F14.6)') i,ebmin,ebmax,ebmax-ebmin
      enddo
    end if

111 format(a,2i4)
112 format('  BandGap ', a4,' = ',f14.6)
113 format('  Direct gap at k = ',3f14.6,' ik=',i5)
114 format('  BandGap ', a4,' = ',3f14.6)
115 format('  Indirect gap, k(VBM) = ',3f8.3,' ik = ',i5,/,&
    &      '                k(CBM) = ',3f8.3,' ik = ',i5)

end subroutine 
