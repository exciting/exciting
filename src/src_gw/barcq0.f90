!BOP
!
!!ROUTINE: barcq0
!
!!INTERFACE:
!
subroutine barcq0()
!
!!DESCRIPTION:
!
! This subroutine calculates the matrix of the bare coulomb potential for
! q=0 and atomic functions with L=0
!
!!USES:
    use modinput
    use modmain,               only : nspecies, natoms, idxas, zzero, pi, &
    &                                 gkmax, natmtot, nrmt, spr, bvec
    use mod_product_basis,     only : nmix, bigl, maxnmix, umix, mbindex, &
                                      locmatsiz
    use mod_coulomb_potential, only : barc
    use mod_misc_gw,           only : vi, atposl, pia
    use mod_kpointset

!!LOCAL VARIABLES:
    implicit none
    
    type(G_set) :: Gset
    
    integer :: ngrid(3), grid(3,2)
    integer :: ng, ngl, igl, ig
    integer :: ia, is, ias, ja, js, jas
    integer :: nr, ir, irm, jrm
    integer :: l1, m1, l2, m2
    integer :: imix, jmix
    
    integer, allocatable :: igl2ig(:), ig2igl(:)
    
    
    real(8) :: gmax, gl0
    real(8) :: dpos(3), gvec(3)
    real(8) :: x, gpr
    
    real(8), allocatable :: sinf(:), sing(:,:,:)
    real(8), allocatable :: fr(:), gr(:), cf(:,:)
    
    complex(8) ::  expg, sum
    complex(8), allocatable :: phase(:,:,:)

    
 
! !REVISION HISTORY:
! 
! Created 16th. March 2004 by RGA
! Last modified 31. March 2005 by RGA
! Revisited: June 2011 by DIN
!
!EOP
!BOC

    ! set the increased Gmax cutoff
    !gmax = 8*input%groundstate%gmaxvr
    gmax = 10*input%gw%MixBasis%gmb*gkmax
    
    ngrid(1) = idint(gmax*pia(1))+1
    ngrid(2) = idint(gmax*pia(2))+1
    ngrid(3) = idint(gmax*pia(3))+1
    do ir = 1, 3
      grid(ir,1) = -ngrid(ir)
      grid(ir,2) =  ngrid(ir)
    end do

    call generate_G_vectors(Gset,bvec,grid,gmax)
    ! total number of G vectors
    ng = Gset%ngvec
    
    !----------------------------------------------------
    ! find the number of different length G-vectors and 
    ! their position in Gset
    !----------------------------------------------------
    ! mapping igl -> ig
    allocate(igl2ig(ng))
    igl2ig(:) = 0
    ! mapping ig -> igl
    allocate(ig2igl(ng))
    ig2igl(:) = 0
        
    igl = 1
    gl0 = Gset%gc(1)
    do ig = 2, ng
      if (dabs(Gset%gc(ig)-gl0)>1.d-8) then
        igl = igl+1
        gl0 = Gset%gc(ig)
        ! store the index of this element
        igl2ig(igl) = ig
      end if
      ig2igl(ig) = igl
    end do
    ngl = igl
    
    !-------------------------
    ! Calculate phase factors
    !-------------------------
    allocate(phase(ngl,natmtot,natmtot))
    phase(:,:,:) = zzero
    
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        
        do js = 1, nspecies
          do ja = 1, natoms(js)
            jas = idxas(ja,js)
            dpos(1:3) = atposl(1:3,ia,is)-atposl(1:3,ja,js)

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ig,igl,gvec,gpr,expg)
!$OMP DO
#endif            
            do ig = 2, ng
              igl = ig2igl(ig)
              gvec(1:3) = dble(Gset%ivg(1:3,ig))
              gpr = gvec(1)*dpos(1)+ &
              &     gvec(2)*dpos(2)+ &
              &     gvec(3)*dpos(3)
              expg = cmplx(cos(2.0d0*pi*gpr),-sin(2.0d0*pi*gpr),8)
              phase(igl,ias,jas) = phase(igl,ias,jas)+expg
            enddo ! igl
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif            
            
          enddo ! ja
        enddo ! js
        
      enddo ! ia
    enddo ! is

    !---------------------------
    ! Calculate singular terms
    !---------------------------
    allocate(sing(ngl,maxnmix,natmtot)) 
    sing(:,:,:) = 0.0d0

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(igl,ig,gl0,is,nr,sinf,fr,gr,cf,ir,x,ia,ias,irm,l1)
!$OMP DO
#endif      
    do igl = 2, ngl
      ig = igl2ig(igl)
      gl0 = Gset%gc(ig)
        
      do is = 1, nspecies
        nr = nrmt(is)
        allocate(sinf(nr))          
        allocate(fr(nr),gr(nr),cf(3,nr))
    
        ! j_{l=0}.r
        do ir = 1, nr
          x = spr(ir,is)*gl0
          sinf(ir) = sin(x)
        end do
        
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do irm = 1, nmix(ias)
            l1 = bigl(irm,ias)
            if (l1==0) then
              do ir = 1, nr
                fr(ir) = umix(ir,irm,ias)*sinf(ir)
              end do
              call fderiv(-1,nr,spr(:,is),fr,gr,cf)
              sing(igl,irm,ias) = gr(nr)
            end if
          end do ! irm
        enddo ! ia
      
        deallocate(fr,gr,cf)
        deallocate(sinf)
      enddo ! is
        
    enddo ! igl
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif       

    !--------------------------------------------------------
    ! Calculate the matrix elements of the Coulomb potential 
    !--------------------------------------------------------
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(imix,is,ia,ias,irm,l1,m1,jmix,js,ja,jas,jrm,l2,m2,sum,igl,ig,gl0)
!$OMP DO
#endif     
    do imix = 1, locmatsiz
      is  = mbindex(imix,1)
      ia  = mbindex(imix,2)
      ias = idxas(ia,is)
      irm = mbindex(imix,3)
      l1  = mbindex(imix,4)
      m1  = mbindex(imix,5)
      if (l1==0) then
      
        do jmix = 1, locmatsiz
          js  = mbindex(jmix,1)
          ja  = mbindex(jmix,2)
          jas = idxas(ja,js)
          jrm = mbindex(jmix,3)
          l2  = mbindex(jmix,4)
          m2  = mbindex(jmix,5)
          if (l2==0) then
            
            ! sum over G
            sum = zzero
            do igl = 2, ngl
              ig = igl2ig(igl)
              gl0 = Gset%gc(ig)
              sum = sum + &
              &     1.d0/(gl0**4)     * &
              &     phase(igl,ias,jas)* &
              &     sing(igl,irm,ias) * &
              &     sing(igl,jrm,jas)
            end do ! igl
            
            barc(imix,jmix) = 16.d0*pi*pi*vi*sum
            
          end if ! l2==0
        enddo ! jmix
      
      end if ! l1==0
    enddo ! imix
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

    deallocate(sing)
    deallocate(phase)
    call delete_G_vectors(Gset)

    return
end subroutine
!EOC
