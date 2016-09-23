!BOP
!
!!ROUTINE: calcmpwmix
!
!!INTERFACE:
subroutine calcmpwmix(iq)
      
!!DESCRIPTION:
!
! This subroutine calculates the matrix elements between mixed basis
! functions and plane waves.
!
!!USES:
    use modmain,               only : pi, nspecies, nrmt, spr, natoms, &
    &                                 idxas, idxlm, zi, zzero, zone
    use modgw,                 only : Gset, Gamma, kqset, Gqset, Gqbarc
    use mod_product_basis,     only : matsiz, maxbigl, nmix, bigl, umix, &
    &                                 locmatsiz, mpwipw, mpwmix
    use mod_coulomb_potential, only : wi0
    use mod_misc_gw,           only : vi, atposl, alat
    
!!INPUT PARAMETERS:      
    implicit none
    integer(4), intent(in) :: iq

!!LOCAL VARIABLES:
    integer :: npw, ipw, ipw0  ! PW
    integer :: ngq, igq        ! IPW
    integer :: imix, is, ia, ias, irm, l1, m1, l1m1
    integer :: ir, nr
    real(8) :: const, x
    real(8) :: gvec(3), gqvec(3), gqlen, gpr
    real(8) :: janl
    real(8), allocatable :: fr(:), gr(:), cf(:,:)
    real(8), allocatable :: bessl(:,:)
    complex(8) :: expg, prefac
    complex(8) :: sph((maxbigl+1)*(maxbigl+1))
    complex(8), allocatable :: tmat1(:,:), tmat2(:,:)
 
!!REVISION HISTORY:
!   Created: Mar 2014 by DIN

!EOP
!BOC
    npw = Gqbarc%ngk(1,iq)

    if (allocated(mpwmix)) deallocate(mpwmix)
    allocate(mpwmix(matsiz,npw))
    mpwmix(:,:) = zzero
      
    const = 4.d0*pi*sqrt(vi)
      
    if (Gamma) then
      ipw0 = 2
      call calcwmix0
      mpwmix(1:matsiz,1) = wi0(1:matsiz)
      deallocate(wi0)
    else
      ipw0 = 1
    end if
    
  !===========================
  ! PW - MB-MT overlap matrix
  !===========================
    
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipw,gvec,gqvec,gqlen,sph,imix,is,nr,bessl,ir,x,fr,gr,cf,ia,ias,gpr,expg,irm,l1,janl,prefac,m1,l1m1)
!$OMP DO
#endif
    do ipw = ipw0, npw
    
      ! G vector (lattice)
      gvec(1:3) = dble(Gset%ivg(1:3,Gqbarc%igkig(ipw,1,iq)))
      ! G+q vector
      gqvec(1:3) = Gset%vgc(1:3,Gqbarc%igkig(ipw,1,iq))+kqset%vqc(1:3,iq)
      ! length of G+q vector
      gqlen = dsqrt(gqvec(1)*gqvec(1)+gqvec(2)*gqvec(2)+gqvec(3)*gqvec(3))
      
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
          expg = cmplx(cos(2.d0*pi*gpr),sin(2.d0*pi*gpr),8)
          
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
            
            prefac = cmplx(const*janl,0.d0,8)*expg*(zi**l1)
            do m1 = -l1, l1
              imix = imix + 1
              l1m1 = idxlm(l1,m1)
              mpwmix(imix,ipw) = prefac*conjg(sph(l1m1))
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
    
  !===========================
  ! PW - MB-PW overlap matrix
  !===========================

    do imix = locmatsiz+1, matsiz
      do ipw = 1, npw
        mpwmix(imix,ipw) = mpwipw(imix-locmatsiz,ipw)
      end do
    end do  
    
    return
end subroutine
!EOC
