
subroutine calcbarcmb_ipw_mt(iq)

    use modinput
    use modmain,               only : nspecies, nrmt, spr, natoms, &
    &                                 idxas, idxlm, nrmtmax
    use modgw,                 only : Gset, Gamma, kqset, Gqset, Gqbarc
    use mod_product_basis,     only : locmatsiz, maxbigl, nmix, bigl, umix, &
    &                                 mpwipw, mpwmix
    use mod_coulomb_potential, only : barc
    use mod_misc_gw,           only : vi, atposl
    use constants,             only : zi, zzero, zone, pi
    implicit none
    ! input variables
    integer, intent(in) :: iq
    ! local variables
    integer :: npw, ipw, ipw0  ! PW
    integer :: ngq, igq        ! IPW
    integer :: imix, is, ia, ias, irm, l1, m1, l1m1
    integer :: ir, nr
    real(8) :: const, x
    real(8) :: gvec(3), gqvec(3), gqlen, gpr
    real(8) :: janl
    real(8) :: vc
    real(8), allocatable :: fr(:), gr(:), cf(:,:)
    real(8), allocatable :: bessl(:,:)
    complex(8) :: expg, prefac
    complex(8) :: sph((maxbigl+1)*(maxbigl+1))
    complex(8), allocatable :: tmat1(:,:), tmat2(:,:)
    ! external routine 
    external zgemm
    integer :: polyord
    parameter (polyord=3)
    integer :: ipiv(polyord+1),io1,io2
    integer :: mid,info,lwork
    real(8), allocatable :: poly(:,:),weight(:,:),ints(:),abscissa(:),work(:)
    logical :: usesplines

    usesplines = .false.

    if (.not.usesplines) then
      allocate(poly(0:polyord,0:polyord))
      allocate(ints(0:polyord))
      allocate(abscissa(0:polyord))

      allocate(work(1))
      call dgetri(polyord+1,poly,polyord+1,ipiv,work,-1,info)
      lwork=int(work(1))
      deallocate(work)
      allocate(work(lwork))
     
      allocate(weight(nrmtmax,nspecies))
      weight=0d0

      do is=1,nspecies
        nr=nrmt(is)
        mid=polyord/2
        do ir=2,nr-polyord-1
          abscissa(0:polyord)=spr(ir:ir+polyord,is)-spr(ir,is)
          poly(:,0)=1d0
          do io2=1,polyord
            do io1=0,polyord
              poly(io1,io2)=poly(io1,io2-1)*abscissa(io1)
            enddo
          enddo
          do io1=0,polyord
            ints(io1)=(poly(mid+1,io1)*abscissa(mid+1)-poly(mid,io1)*abscissa(mid))/dble(io1+1)
          enddo
          call dgetrf(polyord+1,polyord+1,poly,polyord+1,ipiv,info)
          call dgetri(polyord+1,poly,polyord+1,ipiv,work,lwork,info)
          call dgemm('N','N',1,polyord+1,polyord+1,1d0,ints(0),1,poly(0,0),polyord+1,1d0,weight(ir,is),1)
        enddo

        abscissa(0:polyord)=spr(1:1+polyord,is)-spr(1,is)
        poly(:,0)=1d0
        do io2=1,polyord
          do io1=0,polyord
            poly(io1,io2)=poly(io1,io2-1)*abscissa(io1)
          enddo
        enddo
        do io1=0,polyord
          ints(io1)=(poly(mid+1,io1)*abscissa(mid+1)-poly(0,io1)*abscissa(0))/dble(io1+1)
        enddo
        call dgetrf(polyord+1,polyord+1,poly,polyord+1,ipiv,info)
        call dgetri(polyord+1,poly,polyord+1,ipiv,work,lwork,info)
        call dgemm('N','N',1,polyord+1,polyord+1,1d0,ints(0),1,poly(0,0),polyord+1,1d0,weight(1,is),1)

        abscissa(0:polyord)=spr(nr-polyord:nr,is)-spr(nr-polyord,is)
        poly(:,0)=1d0
        do io2=1,polyord
          do io1=0,polyord
            poly(io1,io2)=poly(io1,io2-1)*abscissa(io1)
          enddo
        enddo
        do io1=0,polyord
          ints(io1)=(poly(polyord,io1)*abscissa(polyord)-poly(mid,io1)*abscissa(mid))/dble(io1+1)
        enddo
        call dgetrf(polyord+1,polyord+1,poly,polyord+1,ipiv,info)
        call dgetri(polyord+1,poly,polyord+1,ipiv,work,lwork,info)
        call dgemm('N','N',1,polyord+1,polyord+1,1d0,ints(0),1,poly(0,0),polyord+1,1d0,weight(nr-polyord,is),1)
      enddo
      deallocate(poly,abscissa,ints,work)

    endif
    
    npw = Gqbarc%ngk(1,iq)
    ngq = Gqset%ngk(1,iq)
    
    const = 16.0*pi*pi*sqrt(vi)
    
    ! local array
    allocate(tmat1(npw,locmatsiz))
    tmat1(:,:) = zzero
    
    ipw0 = 1
    if (Gamma) ipw0 = 2
    
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipw,gvec,gqvec,gqlen,vc,sph,imix,is,nr,bessl,ir,x,fr,gr,cf,ia,ias,gpr,expg,irm,l1,janl,prefac,m1,l1m1)
#endif
#ifdef USEOMP
!$OMP DO
#endif
    do ipw = ipw0, npw
    
      ! G vector (lattice)
      gvec(1:3) = dble(Gset%ivg(1:3,Gqbarc%igkig(ipw,1,iq)))
      ! G+q vector
      gqvec(1:3) = Gset%vgc(1:3,Gqbarc%igkig(ipw,1,iq))+kqset%vqc(1:3,iq)
      ! length of G+q vector
      gqlen = dsqrt(gqvec(1)*gqvec(1)+gqvec(2)*gqvec(2)+gqvec(3)*gqvec(3))
      if (abs(gqlen) < 1.d-8) then
        write(*,*) 'WARNING(calcbarcmb_ipw_mt.f90): Zero length vector!' 
        cycle
      endif
     
      vc = 1.d0/(gqlen*gqlen)
      
      !----------------------
      ! Calculate Y_lm(q+G)
      !----------------------
      call ylm(gqvec,maxbigl,sph)
      
      ! loop over atoms
      imix = 0
      do is = 1, nspecies
        nr = nrmt(is)
        
        !------------------------------------------------------------
        ! Calculate the spherical Bessel function at each mesh point
        !------------------------------------------------------------
        allocate(bessl(nr,0:maxbigl))
        do ir = 1, nr
          x = spr(ir,is)*gqlen
          call sbessel(maxbigl,x,bessl(ir,0:maxbigl))
        end do ! ir
        allocate(fr(nr),gr(nr),cf(3,nr))
        
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          
          ! exp^{-G.r_a}
          gpr = gvec(1)*atposl(1,ia,is)+ &
          &     gvec(2)*atposl(2,ia,is)+ &
          &     gvec(3)*atposl(3,ia,is)
          expg = cmplx(cos(2.d0*pi*gpr),-sin(2.d0*pi*gpr),8)
          
          do irm = 1, nmix(ias)
            
            !---------------------------------------
            ! Calculate the radial integral J_{aNL}
            !---------------------------------------
            l1 = bigl(irm,ias)
            if (usesplines) then
              do ir = 1, nr
                fr(ir) = bessl(ir,l1)*umix(ir,irm,ias)*spr(ir,is)
              end do ! ir
              call fderiv(-1,nr,spr(1:nr,is),fr,gr,cf)
              janl = gr(nr)
            else
              janl=0d0
              Do ir = 1, nr
                janl=janl+bessl(ir,l1)*umix(ir,irm,ias)*spr(ir,is)*weight(ir,is)
              End Do              
            endif
            prefac = cmplx(const*janl*vc,0.d0,8)*expg*((-zi)**l1)
            do m1 = -l1, l1
              imix = imix + 1
              l1m1 = idxlm(l1,m1)
              tmat1(ipw,imix) = prefac*sph(l1m1)
            end do ! m1
          end do ! irm
          
        end do ! ia
        
        deallocate(bessl)
        deallocate(fr,gr,cf)
      end do ! is
      
    end do ! ig
#ifdef USEOMP
!$OMP END DO
#endif


#ifdef USEOMP
!$OMP END PARALLEL
#endif
    if (allocated(weight)) deallocate(weight)

    !-------------------------------------------------
    ! \Sum_G' \tilde{S}^{*}_{G'i}(q) J_{aNL}(|G'+q|) 
    !-------------------------------------------------
    allocate(tmat2(1:ngq,1:locmatsiz))  
    call zgemm( 'n','n',ngq,locmatsiz,npw, &
    &           zone,mpwipw,ngq, &
    &           tmat1,npw, &
    &           zzero,tmat2,ngq)
    deallocate(tmat1)
    
    do igq = 1, ngq
      do imix = 1, locmatsiz
        barc(locmatsiz+igq,imix) = tmat2(igq,imix)
        barc(imix,locmatsiz+igq) = conjg(tmat2(igq,imix))
      end do
    end do
    
    ! deallocate local arrays
    deallocate(tmat2)
    
    !write(*,*) 'IPW-MT'
    !do igq = 1, ngq, ngq/10
    !do imix = 1, locmatsiz, locmatsiz/10
    !  write(*,*) igq, imix, barc(locmatsiz+igq,imix)
    !end do
    !end do
    
end subroutine
