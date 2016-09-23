!
subroutine calcbarcmb_ipw_mt(iq)
!
    use modinput
    use modmain,               only : pi, nspecies, nrmt, spr, natoms, &
    &                                 idxas, idxlm, zi, zzero, zone
    use modgw,                 only : Gset, Gamma, kqset, Gqset, Gqbarc
    use mod_product_basis,     only : locmatsiz, maxbigl, nmix, bigl, umix, &
    &                                 mpwipw, mpwmix
    use mod_coulomb_potential, only : barc, rccut, i_sz, vccut
    use mod_misc_gw,           only : vi, atposl
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
    real(8) :: vc, kxy, kz
    real(8), allocatable :: fr(:), gr(:), cf(:,:)
    real(8), allocatable :: bessl(:,:)
    complex(8) :: expg, prefac
    complex(8) :: sph((maxbigl+1)*(maxbigl+1))
    complex(8), allocatable :: tmat1(:,:), tmat2(:,:)
    ! external routine 
    external zgemm
    
    npw = Gqbarc%ngk(1,iq)
    ngq = Gqset%ngk(1,iq)
    
    const = 16.0*pi*pi*sqrt(vi)
    
    ! local array
    allocate(tmat1(npw,locmatsiz))
    tmat1(:,:) = zzero
    
    ipw0 = 1
    if (Gamma) then
      ipw0 = 2
      !if (vccut) tmat1(1,1:locmatsiz) = i_sz*mpwmix(1:locmatsiz,1)
    end if
    
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipw,gvec,gqvec,gqlen,vc,kxy,kz,sph,imix,is,nr,bessl,ir,x,fr,gr,cf,ia,ias,gpr,expg,irm,l1,janl,prefac,m1)
!$OMP DO
#endif
    do ipw = ipw0, npw
    
      ! G vector (lattice)
      gvec(1:3) = dble(Gset%ivg(1:3,Gqbarc%igkig(ipw,1,iq)))
      ! G+q vector
      gqvec(1:3) = Gset%vgc(1:3,Gqbarc%igkig(ipw,1,iq))+kqset%vqc(1:3,iq)
      ! length of G+q vector
      gqlen = dsqrt(gqvec(1)*gqvec(1)+gqvec(2)*gqvec(2)+gqvec(3)*gqvec(3))
      
      if (vccut) then
     
        ! apply cutoff
        select case (trim(input%gw%barecoul%cutofftype))
      
          case('0d')
            vc = 1.d0/(gqlen*gqlen)
            vc = vc*(1.d0-dcos(gqlen*rccut))
      
          case ('2d')
            ! version by Ismail-Beigi (fixed rc = L_z/2)
            kxy = dsqrt(gqvec(1)*gqvec(1)+gqvec(2)*gqvec(2))
            kz = dabs(gqvec(3))
            vc = 1.d0/(gqlen*gqlen)
            vc = vc*(1.d0-dexp(-kxy*rccut)*dcos(kz*rccut))
          
          case default
            write(*,*) 'ERROR(calcbarcmb_ipw_mt): Specified cutoff type is not implemented!'
          
        end select
      
      else
      
        ! no cutoff
        vc = 1.d0/(gqlen*gqlen)
       
      end if
      
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
            do ir = 1, nr
              fr(ir) = bessl(ir,l1)*umix(ir,irm,ias)*spr(ir,is)
            end do ! ir
            call fderiv(-1,nr,spr(1:nr,is),fr,gr,cf)
            janl = gr(nr)
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
!$OMP END PARALLEL
#endif

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

