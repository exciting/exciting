subroutine bandanaly(ib,nb,nkp,kvecs,eband,efermi,title,file)

      implicit none
      integer, intent(in) :: ib,nb,nkp,file
      real(8), intent(in) :: eband(ib:nb,nkp),efermi,kvecs(3,nkp)
      character(len=*),intent(in):: title


      integer :: i,iat,ik,nc,nv,ikv,ikc     ! Indexes inequivalent atoms
      integer :: fid,info
      integer :: nat,nwf,nbk
      integer :: numin,nomax,ikvm,ikcm
      real(8) :: kvec(3),w,kvec0(3)
      real(8) :: evbm,ecbm,egap(3),ebmin,ebmax
      real(8) :: bande(ib:nb,nkp)
      logical :: lmetal

      real(8),parameter::hev=27.2113961

      call boxmsg(file,'-',trim(title)//" Band Analysis")
      write(file,'(a,2i5)'  ) "  Range of bands considered: ", ib, nb
      write(file,'(a,f10.4,a)') "  EFermi= ",efermi, '[Ha]'
      bande=eband
!
!     Find Valence Band Maximum (VBM) and Cond. Band Minimum (CBM)
!
      if(maxval(bande).lt.efermi .or. minval(bande).gt.efermi ) then 
        write(file,*) "WARNING from bandanaly: "
        write(file,*) " - Fermi energy outside the energy range of bande!" 
        write(file,'(a,f10.4)') "  Fermi Energy(eV):   ",efermi*hev
        write(file,'(a,f10.4)') " - Maximal energy:    ",maxval(bande)*heV
        write(file,'(a,f10.4)') " - Minimal energy:    ",minval(bande)*heV
        return 
      endif 

      nomax=0
      numin=1000
      ikvm = 1
      ikcm = 1
      lmetal = .false.
      do ik = 1, nkp
         nv=0
         nc=1000
         do i=ib,nb
            if(bande(i,ik).lt.efermi)then
               if(i.gt.nv) nv = i 
            else
               if(i.lt.nc) nc = i  
            endif
         enddo
         if(nv.gt.nomax) nomax = nv 
         if(nc.lt.numin) numin = nc 
          nv = nomax
          nc = numin 
          ikv=ikvm
          ikc=ikcm
          if(bande(nv,ik).gt.bande(nv,ikv)) ikvm = ik
          if(bande(nc,ik).lt.bande(nc,ikc)) ikcm = ik
      enddo
      
      if(nomax.lt.1) then
        write(file,*) "ERROR: nomax < 1  !!!"
        write(file,*) "  --- Check the Fermi energy"
        stop
      endif

      if(nomax.ge.nb) then
        write(file,*) "ERROR: nomax is larger than the number of available bands !!!"
        write(file,*) "  --- Check the Fermi energy"
        stop
      endif

      if(numin.lt.1) then
        write(file,*) "ERROR: numin <1  !!!"
        write(file,*) "  --- Check the Fermi energy"
        stop
      endif

      if(numin.ge.nb) then
        write(file,*) "ERROR: numin >= nb  !!!"
        write(file,*) "  --- Check the Fermi energy"
        stop
      endif

      if(nomax.ge.numin) then
        write(file,*) "Valence and Conductance bands overlap: metallic!"
        lmetal=.true.
      endif
!
!     set evbm
!
      if(lmetal) then
        evbm = efermi 
      else 
        evbm = bande(nomax,ikvm)
      endif 
      bande = (bande-evbm)*hev

      nv = nomax
      nc = numin
      ikv = ikvm 
      ikc = ikcm 
      write(file,101) "Band index for VBM and CBM=",nv,nc

      egap(1)= bande(nc,ikc)-bande(nv,ikv) 
      egap(2)= minval(bande(nc:nb,ikv))-bande(nv,ikv)
      egap(3)= bande(nc,ikc)-maxval(bande(ib:nv,ikc))

      if(ikvm.eq.ikcm) then ! direct gap 
        write(file,112) egap(1)
        write(file,113) kvecs(:,ikv),ikv
      else 
        write(file,114) egap(1:3)
        write(file,115) kvecs(:,ikv),ikv,kvecs(:,ikc),ikc
      endif

      write(file,*) "Range of each band with respect to VBM (eV):"
      write(file,'(a5,3a12)') 'n ','Bottom','Top','Width'
      do i=ib,nb
        ebmin=minval(bande(i,1:nkp))
        ebmax=maxval(bande(i,1:nkp))
        write(file,'(i5,3F12.3)') i,ebmin,ebmax,ebmax-ebmin
      enddo

  100 format(3e19.12,a10,2i6,f5.1)
  101 format(a,2i4)
 112  format(':BandGap =  ',f12.3,' eV')
 113  format("  Direct gap at k=  ",3f8.3,' ik=',i5)
 114  format(':BandGap =  ',3f12.3,' eV')
 115  format("  Indirect gap, k(VBM)=",3f8.3,' ik=',i5,/,&
     &       "                k(CBM)=",3f8.3,' ik=',i5)

end subroutine 
