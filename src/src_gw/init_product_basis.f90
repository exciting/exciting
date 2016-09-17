
subroutine init_product_basis()

    use modinput
    use modmain,               only: nlomax, natmtot, nrmtmax, nspecies, &
    &                                natoms, idxas, spnstmax, apwordmax, &
    &                                chgcr, idxlm
    use modgw,                 only: Gqset, fgw
    use mod_core_states,       only: ncmax, init_core_states
    use mod_product_basis,     only: maxnup, nmix, umix, bigl, mbl, maxbigl, &
    &                                maxnmix, matsizmax, &
    &                                bradketc, bradketa, bradketlo, &
    &                                lmixmax, locmatsiz, locmixind, mbindex, &
    &                                rtl, rrint
    use mod_mpi_gw,            only: myrank, nproc_tot
    use modmpi,                only: mpi_allgatherv_ifc, barrier
    use reallocate
    implicit none

    integer(4) :: is, ia, ias
    integer(4) :: irm, im, imix
    integer(4) :: l, m, lm
    integer(4) :: ndim, nrwf
    
!_______________________________________________________________________________    
! Generate Kohn-Sham radial functions

    ! generate the core wavefunctions and densities
    call gencore
    ! find the linearisation energies
    call linengy
    ! generate APW radial functions
    call genapwfr
    ! generate local-orbital radial functions
    call genlofr
    
    ! initialize core states
    if (dabs(chgcr)>1.d-6) then
      call init_core_states
    else
      input%gw%coreflag = "vab"
    end if

!_______________________________________________________________________________    
! Generate mixed product basis functions

    ! Estimate the maximum possible number of product functions
    maxnup = 2*(ncmax+input%gw%MixBasis%lmaxmb+nlomax+1)* &
    &          (input%gw%MixBasis%lmaxmb+nlomax+1)
    
    if (allocated(nmix)) deallocate(nmix)
    allocate(nmix(natmtot))
    nmix(:) = 0
    
    if (associated(umix)) deallocate(umix)
    allocate(umix(nrmtmax,maxnup,natmtot))
    umix(:,:,:) = 0.d0
    
    if (associated(bigl)) deallocate(bigl)
    allocate(bigl(maxnup,natmtot))
    bigl(:,:) = 0
    
    if (allocated(mbl)) deallocate(mbl)
    allocate(mbl(natmtot))
    mbl(:) = 0
    
    ! loop over atoms
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        ! Generate all possible radial function products
        call setuprod(ia,is)
        ! Calculate the radial part of the mixed basis functions
        call setumix(ia,is)
      end do ! ia
    end do ! is
    
    ! Reallocate the arrays to their actual size
    maxnmix = maxval(nmix)
    maxbigl = maxval(mbl)
    call doreallocate(umix,nrmtmax,maxnmix,natmtot)
    call doreallocate(bigl,maxnmix,natmtot)
    if (maxbigl>input%groundstate%lmaxapw) then
      if (allocated(idxlm)) deallocate (idxlm)
      allocate(idxlm(0:maxbigl,-maxbigl:maxbigl))
      lm = 0
      do l = 0, maxbigl
        do m = -l, l
          lm = lm + 1
          idxlm(l,m) = lm
        end do
      end do
    end if
    
    !--------------------------------
    ! pre-calculate radial integrals
    !--------------------------------
    if (allocated(rtl)) deallocate(rtl)
    allocate(rtl(maxnmix,natmtot))
    rtl(:,:) = 0.0d0
    
    if (allocated(rrint)) deallocate(rrint)
    allocate(rrint(maxnmix*(maxnmix+1)/2,natmtot))
    rrint(:,:) = 0.0d0

    nrwf = max(spnstmax,input%groundstate%lmaxapw,nlomax)
    
    if (allocated(bradketc)) deallocate(bradketc)
    allocate(bradketc(3,maxnmix,spnstmax,0:nrwf,apwordmax,natmtot))
    bradketc = 0.d0
    
    if (allocated(bradketa)) deallocate(bradketa)
    allocate(bradketa(3,maxnmix,0:input%groundstate%lmaxapw, &
    &                 apwordmax,0:nrwf,apwordmax,natmtot))
    bradketa = 0.d0
    
    if (allocated(bradketlo)) deallocate(bradketlo)
    allocate(bradketlo(3,maxnmix,nlomax,0:nrwf,apwordmax,natmtot))
    bradketlo = 0.d0
    
    ! loop over atoms
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        ! Calculate the matrix elements <r^L>_{aNL} and <r_<^L/r_>^{L+1}> _{aNL,N'L}:
        call calcrlamint(ia,is)
        ! Calculate the matrix elements <NL,lambda|lambda'>:
        call calcbradket(ia,is)
      end do ! ia 
    end do ! is
    
    !------------------------------------------------------------------
    ! Calculate the total number of mixed wave functions (including M)    
    ! = size of the MT product basis functions
    !------------------------------------------------------------------
    lmixmax = 0
    locmatsiz = 0
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        imix = 0
        do irm = 1, nmix(ias)
          l = bigl(irm,ias)
          imix = imix+(2*bigl(irm,ias)+1)
        end do  
        if (imix>lmixmax) lmixmax = imix
        locmatsiz = locmatsiz+imix
      end do ! ia
    end do ! ias
    
    !-------------------------------------------------------------------
    ! Set an array that stores the general index of the mixed function 
    ! for a given mixed function of a given atom
    !-------------------------------------------------------------------
    if (allocated(locmixind)) deallocate(locmixind)
    allocate(locmixind(natmtot,lmixmax))
    locmixind(:,:) = 0
    
    im = 0
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        imix = 0
        do irm = 1, nmix(ias)
          l = bigl(irm,ias)
          do m = -l, l
            imix = imix+1
            im = im+1
            locmixind(ias,imix) = im
          end do ! m
        end do  
      end do ! ia
    end do ! ias
    
    ! The maximum size of MB basis
    matsizmax = locmatsiz+Gqset%ngkmax
    
    !-------------------------------------------------------------------
    ! mapping: MB function index -> (aNLM)
    !-------------------------------------------------------------------
    if (allocated(mbindex)) deallocate(mbindex)
    allocate(mbindex(locmatsiz,5))
    mbindex(:,:) = 0
    im = 0
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        do irm = 1, nmix(ias)
          l = bigl(irm,ias)
          do m = -l, l
            im = im+1
            mbindex(im,1) = is
            mbindex(im,2) = ia
            mbindex(im,3) = irm
            mbindex(im,4) = l
            mbindex(im,5) = m
          end do ! m
        end do ! irm
      end do ! ia
    end do ! is
    
    !----------------------------------------------------------------------
    ! Calculate the coefficients tildeg needed for the structure constants
    !----------------------------------------------------------------------
    call calctildeg(2*(input%gw%MixBasis%lmaxmb+1))

    ! Print info
    if (myrank==0) then
      call boxmsg(fgw,'-',"Mixed product WF info")
      write(fgw,*) ' Maximal number of MT wavefunctions per atom: ', lmixmax 
      write(fgw,*) ' Total number of MT wavefunctions:            ', locmatsiz
      write(fgw,*) ' Maximal number of PW wavefunctions:          ', Gqset%ngkmax 
      write(fgw,*) ' Total number of mixed-product wavefunctions: ', matsizmax
      write(fgw,*)
    end if
    
    return
end subroutine
