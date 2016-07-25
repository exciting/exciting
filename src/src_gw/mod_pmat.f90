
module mod_pmat

  implicit none
    
  !---------------------------------
  ! C_{G+k}^n \cdot A_{alm}^{\zeta}
  !---------------------------------
  complex(8), allocatable :: apwcmt(:,:,:,:)
  complex(8), allocatable :: locmt(:,:,:,:)
    
  !--------------------------------------------------
  ! radial matrix elements:  < u_2 | \nabla | u_1 >  
  !--------------------------------------------------
  real(8), allocatable :: uadua(:,:,:,:,:,:)
    
  real(8), allocatable :: uadulo(:,:,:,:,:,:)
  real(8), allocatable :: ulodua(:,:,:,:,:,:)
  real(8), allocatable :: ulodulo(:,:,:,:,:,:)
    
  real(8), allocatable :: uaduc(:,:,:,:,:,:)
  real(8), allocatable :: uloduc(:,:,:,:,:,:)

  integer, allocatable :: lmmaxlo(:)
  integer, allocatable :: lmapwidx(:,:)
  integer, allocatable :: lmloidx(:,:,:)
    
  private :: apwcmt, locmt, uadua, uadulo, ulodua, ulodulo, uaduc, uloduc
  private :: lmmaxlo, lmapwidx, lmloidx
    
contains

!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
  subroutine init_pmat(lcore)
    use modinput
    use modmain, only : nstsv, apwordmax, lmmaxapw, natmtot, &
    &                   nlotot, nlomax, lolmax, nspecies, &
    &                   nlorb, lorbl, apword
    implicit none
    logical, intent(in) :: lcore
    integer :: is, ilo, l, m, lm, io
    ! LM-mapping for APW
    if (allocated(lmapwidx)) deallocate(lmapwidx)
    allocate(lmapwidx(2,lmmaxapw))
    lm = 0
    do l = 0, input%groundstate%lmaxapw
    do m = -l, l
      lm = lm+1
      lmapwidx(1,lm) = l
      lmapwidx(2,lm) = m
    end do
    end do
    ! LM-mapping for LO
    if (allocated(lmmaxlo)) deallocate(lmmaxlo)
    allocate(lmmaxlo(nspecies))
    lmmaxlo(:) = 0
    do is = 1, nspecies
      lm = 0
      do ilo = 1, nlorb(is)
        l = lorbl(ilo,is)
        do m = -l, l
          lm = lm+1
        end do
      end do
      lmmaxlo(is) = lm
    end do
    ! Note that in case of LO lm .ne. idxlm(l,m)
    if (allocated(lmloidx)) deallocate(lmloidx)
    allocate(lmloidx(3,maxval(lmmaxlo),nspecies))
    lmloidx(:,:,:) = 0
    do is = 1, nspecies
      lm = 0
      do ilo = 1, nlorb(is)
        l = lorbl(ilo,is)
        do m = -l, l
          lm = lm+1
          lmloidx(1,lm,is) = ilo
          lmloidx(2,lm,is) = l
          lmloidx(3,lm,is) = m
        end do
      end do
    end do
    ! Precalculated products
    if (allocated(apwcmt)) deallocate(apwcmt)
    allocate(apwcmt(nstsv,apwordmax,lmmaxapw,natmtot))
    if (nlotot>0) then
      if (allocated(locmt)) deallocate(locmt)
      allocate(locmt(nstsv,nlomax,-lolmax:lolmax,natmtot))
    end if
    ! radial integrals
    call genpmatvv_radial
    if (lcore) call genpmatcv_radial
    return
  end subroutine

!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
  subroutine clear_pmat()
    implicit none
    if (allocated(apwcmt)) deallocate(apwcmt)
    if (allocated(locmt))  deallocate(locmt)
    if (allocated(uadua))  deallocate(uadua)
    if (allocated(uadulo)) deallocate(uadulo)
    if (allocated(ulodua)) deallocate(ulodua)
    if (allocated(ulodulo))deallocate(ulodulo)
    if (allocated(uaduc))  deallocate(uaduc)
    if (allocated(uloduc)) deallocate(uloduc)
    if (allocated(lmmaxlo)) deallocate(lmmaxlo)
    if (allocated(lmapwidx)) deallocate(lmapwidx)
    if (allocated(lmloidx)) deallocate(lmloidx)
  end subroutine

!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------  
  subroutine genpmatvv_radial()
    use modinput
    use modmain, only : nrmtmax, apwordmax, lmmaxapw, &
    &                   nlotot, nlomax, lolmax, &
    &                   natmtot, nspecies, natoms, spr, idxas, idxlm, &
    &                   apword, nlorb, lorbl, apwfr, lofr, nrmt
    implicit none
    ! local variables
    integer :: is, ia, ias, nr, ir
    integer :: l1, m1, lm1, l2, m2, lm2
    integer :: ilo1, ilo2, io1, io2, i
    ! automatic arrays
    real(8) :: r2(nrmtmax), fr(nrmtmax), gr(nrmtmax), cf(3,nrmtmax)
    ! allocatable arrays
    real(8), allocatable :: dapwfr(:,:,:,:,:), dlofr(:,:,:,:,:)
    !----------------------------------------------
    ! allocate local arrays for radial derivatives
    !----------------------------------------------
    allocate(dapwfr(lmmaxapw,nrmtmax,3,apwordmax,lmmaxapw))
    dapwfr(:,:,:,:,:) = 0.d0
    if (nlotot>0) then
     allocate(dlofr(lmmaxapw,nrmtmax,3,nlomax,-lolmax:lolmax))
     dlofr(:,:,:,:,:) = 0.d0
    end if
    !------------------------------
    ! initialize radial integrals
    !------------------------------
    if (allocated(uadua)) deallocate(uadua)
    allocate(uadua(apwordmax,lmmaxapw,apwordmax,lmmaxapw,natmtot,3))
    uadua(:,:,:,:,:,:) = 0.d0
    if (nlotot>0) then
      if (allocated(uadulo)) deallocate(uadulo)
      allocate(uadulo(apwordmax,lmmaxapw,nlomax,-lolmax:lolmax,natmtot,3))
      uadulo(:,:,:,:,:,:) = 0.d0
      if (allocated(ulodua)) deallocate(ulodua)
      allocate(ulodua(nlomax,-lolmax:lolmax,apwordmax,lmmaxapw,natmtot,3))
      ulodua(:,:,:,:,:,:) = 0.d0
      if (allocated(ulodulo)) deallocate(ulodulo)
      allocate(ulodulo(nlomax,-lolmax:lolmax,nlomax,-lolmax:lolmax,natmtot,3))
      ulodulo(:,:,:,:,:,:) = 0.d0
    end if
    !------------------------
    ! begin loop over atoms
    !------------------------
    do is = 1, nspecies
      nr = nrmt(is)
      ! calculate r^2  
      do ir = 1, nr
        r2(ir) = spr(ir,is)*spr(ir,is)
      end do
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        !--------------------!
        !     derivatives    !
        !--------------------!
        ! APW functions
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(lm1,l1,m1,io1,fr)
!$OMP DO
#endif
        do lm1 = 1, lmmaxapw
          l1  = lmapwidx(1,lm1)
          m1  = lmapwidx(2,lm1)
          do io1 = 1, apword(l1,is)
            fr(:) = apwfr(:,1,io1,l1,ias)
            call gradzfmtr(input%groundstate%lmaxapw, nr, &
            &              spr(1,is), l1, m1, lmmaxapw, nrmtmax, &
            &              fr, dapwfr(1,1,1,io1,lm1))
          end do
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
        if (nlotot>0) then
          ! local orbital functions
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(lm1,ilo1,l1,m1,fr)
!$OMP DO
#endif              
          do lm1 = 1, lmmaxlo(is)
            ilo1  = lmloidx(1,lm1,is)
            l1    = lmloidx(2,lm1,is)
            m1    = lmloidx(3,lm1,is)
            fr(:) = lofr(:,1,ilo1,ias)
            call gradzfmtr(input%groundstate%lmaxapw, nr, &
            &              spr(1,is), l1, m1, lmmaxapw, nrmtmax, &
            &              fr, dlofr(1,1,1,ilo1,m1))
          end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif              
        end if
        !----------------!
        !     APW-APW    !
        !----------------!
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(lm1,l1,io1,lm2,l2,io2,i,fr,gr,cf)
!$OMP DO
#endif
        do lm1 = 1, lmmaxapw
          l1  = lmapwidx(1,lm1)
          do io1 = 1, apword(l1,is)
            do lm2 = 1, lmmaxapw
              l2  = lmapwidx(1,lm2)
              do io2 = 1, apword(l2,is)
                do i = 1, 3
                  fr(1:nr) = apwfr(1:nr,1,io1,l1,ias)* &
                  &          dapwfr(lm1,1:nr,i,io2,lm2)*r2(1:nr)
                  call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                  uadua(io1,lm1,io2,lm2,ias,i) = gr(nr)
                end do
              end do ! io2
            end do ! lm2
          end do ! io1
        end do ! lm1
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif             
        if (nlotot>0) then
        !----------------------------!
        !     APW-local-orbital      !
        !----------------------------!
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(lm1,l1,io1,lm2,ilo2,m2,i,fr,gr,cf)
!$OMP DO
#endif
          do lm1 = 1, lmmaxapw
            l1  = lmapwidx(1,lm1)
            do io1 = 1, apword(l1,is)
              do lm2 = 1, lmmaxlo(is)
                ilo2 = lmloidx(1,lm2,is)
                m2   = lmloidx(3,lm2,is)
                do i = 1, 3
                  fr(1:nr) = apwfr(1:nr,1,io1,l1,ias)* &
                  &          dlofr(lm1,1:nr,i,ilo2,m2)*r2(1:nr)
                  call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                  uadulo(io1,lm1,ilo2,m2,ias,i) = gr(nr)
                end do
              end do ! lm2
            end do ! io1
          end do ! lm1
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif              
        !----------------------------!
        !     local-orbital-APW      !
        !----------------------------!
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(lm1,ilo1,l1,m1,lm2,l2,io2,i,fr,gr,cf)
!$OMP DO
#endif
          do lm1 = 1, lmmaxlo(is)
            ilo1 = lmloidx(1,lm1,is)
            l1   = lmloidx(2,lm1,is)
            m1   = lmloidx(3,lm1,is)
            do lm2 = 1, lmmaxapw
              l2  = lmapwidx(1,lm2)
              do io2 = 1, apword(l2,is)
                do i = 1, 3
                  fr(1:nr) = lofr(1:nr,1,ilo1,ias)* &
                  &          dapwfr(idxlm(l1,m1),1:nr,i,io2,lm2)*r2(1:nr)
                  call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                  ulodua(ilo1,m1,io2,lm2,ias,i) = gr(nr)
                end do
              end do ! io2
            end do ! lm2
          end do ! lm1
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif              
          !------------------------------------!
          !     local-orbital-local-orbital    !
          !------------------------------------!
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(lm1,ilo1,l1,m1,lm2,ilo2,m2,i,fr,gr,cf)
!$OMP DO
#endif              
          do lm1 = 1, lmmaxlo(is)
            ilo1 = lmloidx(1,lm1,is)
            l1   = lmloidx(2,lm1,is)
            m1   = lmloidx(3,lm1,is)
            do lm2 = 1, lmmaxlo(is)
              ilo2 = lmloidx(1,lm2,is)
              m2   = lmloidx(3,lm2,is)
              do i = 1, 3
                fr(1:nr) = lofr(1:nr,1,ilo1,ias)* &
                &          dlofr(idxlm(l1,m1),1:nr,i,ilo2,m2)*r2(1:nr)
                call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                ulodulo(ilo1,m1,ilo2,m2,ias,i) = gr(nr)
              end do ! i
            end do ! lm2
          end do ! lm1
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif                
        end if ! nlotot>0
      ! end loops over atoms and species
      end do
    end do
    ! deallocate
    deallocate(dapwfr)
    if (nlotot>0) deallocate(dlofr)
    return
  end subroutine

!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
  subroutine genpmatcv_radial()
    use modinput
    use modmain, only : nrmtmax, apwordmax, lmmaxapw, &
    &                   nlotot, nlomax, lolmax, &
    &                   natmtot, nspecies, natoms, spr, idxas, idxlm, &
    &                   apword, nlorb, lorbl, apwfr, lofr, rwfcr, &
    &                   spocc, spl, nrmt
    use mod_core_states, only : ncmax, lcoremax, ncore, ucore
    implicit none
    ! local variables
    integer :: is, ia, ias, nr, ir, ic, ilo, io1, i
    integer :: l1, m1, lm1, l2, m2, lm2
    ! automatic arrays
    real(8) :: r2(nrmtmax)
    real(8) :: fr(nrmtmax), gr(nrmtmax), cf(3,nrmtmax)
    ! allocatable arrays
    real(8), allocatable :: duc(:,:,:,:,:)
    !----------------------------------------------
    ! allocate local arrays for radial derivatives
    !----------------------------------------------
    allocate(duc(lmmaxapw,nrmtmax,3,ncmax,-lcoremax:lcoremax))
    duc(:,:,:,:,:) = 0.d0
    !-------------------
    ! radial integrals
    !-------------------
    allocate(uaduc(apwordmax,lmmaxapw,ncmax,-lcoremax:lcoremax,natmtot,3))
    uaduc(:,:,:,:,:,:) = 0.d0
    if (nlotot>0) then
      allocate(uloduc(nlomax,-lolmax:lolmax,ncmax,-lcoremax:lcoremax,natmtot,3))
      uloduc(:,:,:,:,:,:) = 0.d0
    end if
    !------------------
    ! loop over atoms
    !------------------
    do is = 1, nspecies
      nr = nrmt(is)
      ! calculate r^2
      do ir = 1, nr
        r2(ir) = spr(ir,is)*spr(ir,is)
      end do
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        !--------------------!
        !     derivatives    !
        !--------------------!
        do ic = 1, ncore(is)
          ! QUESTION?: Which core radial function is correct to use:
          ! rwfcr or ucore (renormalized as in GW)?
          fr(:) = ucore(:,1,ic,ias)
          l1 = spl(ic,is)
          do m1 = -l1, l1
            lm1 = idxlm(l1,m1)
            call gradzfmtr(input%groundstate%lmaxapw, nr, &
            &              spr(1,is), l1, m1, lmmaxapw, nrmtmax, &
            &              fr, duc(1,1,1,ic,m1))
          end do ! m1
        end do ! ic
        !---------------------------!
        !     APW-core-orbital      !
        !---------------------------!
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(lm1,l1,io1,ic,l2,m2,i,fr,gr,cf)
!$OMP DO
#endif
        do lm1 = 1, lmmaxapw
          l1  = lmapwidx(1,lm1)
          do io1 = 1, apword(l1,is)
            do ic = 1, ncore(is)
              l2 = spl(ic,is)
              do m2 = -l2, l2
                do i = 1, 3
                  fr(1:nr) = apwfr(1:nr,1,io1,l1,ias)* &
                  &          duc(lm1,1:nr,i,ic,m2)* &
                  &          r2(1:nr)
                  call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                  uaduc(io1,lm1,ic,m2,ias,i) = gr(nr)
                end do ! j
              end do ! m2
            end do ! ic
          end do ! io1
        end do ! lm1
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif            
        if (nlotot>0) then
        !-----------------------------------!
        !     local-orbital-core-orbital    !
        !-----------------------------------!
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(lm1,ilo,l1,m1,ic,l2,m2,i,fr,gr,cf)
!$OMP DO
#endif                
          do lm1 = 1, lmmaxlo(is)
            ilo = lmloidx(1,lm1,is)
            l1  = lmloidx(2,lm1,is)
            m1  = lmloidx(3,lm1,is)            
            do ic = 1, ncore(is)
              l2 = spl(ic,is)
              do m2 = -l2, l2
                do i = 1, 3
                  fr(1:nr) = lofr(1:nr,1,ilo,ias)* &
                  &          duc(idxlm(l1,m1),1:nr,i,ic,m2)* &
                  &          r2(1:nr)
                  call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                  uloduc(ilo,m1,ic,m2,ias,i) = gr(nr)
                end do
              end do
            end do
          end do ! lm1
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif               
        end if ! nlotot>0
      end do ! ia
    end do ! is
    ! deallocate
    deallocate(duc)
    return
  end subroutine

!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
  subroutine genevecalm(ngp,nst,evec,alm)
    use modinput
    use modmain, only : nmatmax, ngkmax, apwordmax, lmmaxapw, natmtot, &
    &                   nspecies, natoms, idxas, idxlm, apword, &
    &                   nlotot, nlomax, lolmax, nlorb, lorbl, idxlo
    implicit none
    ! input parameters
    integer, intent(in) :: ngp
    integer, intent(in) :: nst
    complex(8), intent(in) :: evec(nmatmax,nst)
    complex(8), intent(in) :: alm(ngkmax,apwordmax,lmmaxapw,natmtot)
    ! local variables
    integer :: is, ia, ias, l, m, lm, io, ilo, ist, i
    complex(8), external :: zdotu
    !----------------------------------------------------
    ! generate APW expansion coefficients for muffin-tin 
    !----------------------------------------------------
    apwcmt(:,:,:,:) = 0.d0
    do is = 1, nspecies
    do ia = 1, natoms(is)
      ias = idxas(ia,is)
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(lm,l,io,ist)
!$OMP DO
#endif
      do lm = 1, lmmaxapw
        l = lmapwidx(1,lm)
        do io = 1, apword(l,is)
          do ist = 1, nst
            apwcmt(ist,io,lm,ias) = &
            &  zdotu(ngp,evec(1:ngp,ist),1,alm(1:ngp,io,lm,ias),1)
          end do
        end do
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif      
    end do
    end do
    !----------------------------------------------------
    ! generate local orbital expansion coefficients for muffin-tin
    !----------------------------------------------------
    if (nlotot>0) then
      locmt(:,:,:,:) = 0.d0
      do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(lm,ilo,l,m,i,ist)
!$OMP DO
#endif        
        do lm = 1, lmmaxlo(is)
          ilo = lmloidx(1,lm,is)
          l   = lmloidx(2,lm,is)
          m   = lmloidx(3,lm,is)
          ! important: lm /= idxlm(l,m) (bad name for local orbital (lm) counter)
          i   = idxlo(idxlm(l,m),ilo,ias)
          do ist = 1, nst
            locmt(ist,ilo,m,ias) = evec(ngp+i,ist)
          end do
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif        
      end do
      end do
    end if
    return
  end subroutine
    
  subroutine genpmatvv_k(ngp,igpig,vgpc,nst,evec,nstart,nend,mstart,mend,pmat)
    use modinput
    use modmain, only : ngkmax, nmatmax, idxlm, apword, zzero, zone, &
    &                   nspecies, natoms, idxas, nlotot, nlorb, lorbl, zi, &
    &                   ivg, ivgig, cfunig, lmmaxapw
    implicit none
    ! input parameters
    integer, intent(in)    :: ngp
    integer, intent(in)    :: igpig(ngkmax)
    real(8), intent(in)    :: vgpc(3,ngkmax)
    integer, intent(in)    :: nst
    complex(8), intent(in) :: evec(nmatmax,nst)
    integer, intent(in)    :: nstart, nend
    integer, intent(in)    :: mstart, mend
    complex(8), intent(out):: pmat(nstart:nend,mstart:mend,3)
    ! local variables
    integer :: is, ia, ias
    integer :: l1, m1, lm1, l2, m2, lm2, io1, io2, ilo1, ilo2
    integer :: i, ndim, mdim, n, m
    integer :: igp1, igp2, ig1, ig2, ig, iv(3)
    complex(8) :: zv(nstart:nend)
    ! allocatable arrays
    complex(8), allocatable :: zm(:,:)
    complex(8), allocatable :: cfunt(:,:), pm(:,:)
        
    pmat(:,:,:) = 0.d0

    ndim = nend-nstart+1
    mdim = mend-mstart+1

    !=============================
    ! MT -- MT
    !=============================
    do is = 1, nspecies
    do ia = 1, natoms(is)
      ias = idxas(ia,is)
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,n,i,lm1,l1,m1,io1,ilo1,zv,lm2,l2,m2,io2,ilo2)
!$OMP DO
#endif          
      do m = mstart, mend
      do n = nstart, nend
      do i = 1, 3
          
        !---------------------------!
        !     APW-APW contribution  !
        !---------------------------!
        do lm1 = 1, lmmaxapw
          l1 = lmapwidx(1,lm1)
          do io1 = 1, apword(l1,is)
            do lm2 = 1, lmmaxapw
              l2 = lmapwidx(1,lm2)
              do io2 = 1, apword(l2,is)
                pmat(n,m,i) = pmat(n,m,i)+ &
                &             conjg(apwcmt(m,io1,lm1,ias))* &
                &             uadua(io1,lm1,io2,lm2,ias,i)* &
                &             apwcmt(n,io2,lm2,ias)
              end do ! io2
            end do ! lm2
          end do ! io1
        end do ! lm1
            
        if (nlotot>0) Then
          
          !--------------------------------------!
          !     APW-local-orbital contribution   !
          !--------------------------------------!
          do lm1 = 1, lmmaxapw
            l1 = lmapwidx(1,lm1)
            do io1 = 1, apword(l1,is)
              do lm2 = 1, lmmaxlo(is)
                ilo2 = lmloidx(1,lm2,is)
                m2   = lmloidx(3,lm2,is)
                pmat(n,m,i) = pmat(n,m,i)+ &
                &             conjg(apwcmt(m,io1,lm1,ias))* &
                &             uadulo(io1,lm1,ilo2,m2,ias,i)* &
                &             locmt(n,ilo2,m2,ias)
              end do ! lm2
            end do ! io1 
          end do ! lm1
            
          !--------------------------------------!
          !     local-orbital-APW contribution   !
          !--------------------------------------!
          do lm1 = 1, lmmaxlo(is)
            ilo1 = lmloidx(1,lm1,is)
            m1   = lmloidx(3,lm1,is)
            do lm2 = 1, lmmaxapw
              l2 = lmapwidx(1,lm2)
              do io2 = 1, apword(l2,is)
                pmat(n,m,i) = pmat(n,m,i)+ &
                &             conjg(locmt(m,ilo1,m1,ias))* &
                &             ulodua(ilo1,m1,io2,lm2,ias,i)* &
                &             apwcmt(n,io2,lm2,ias)
              end do ! io2
            end do ! lm2
          end do ! lm1
            
          !------------------------------------------------!
          !     local-orbital-local-orbital contribution   !
          !------------------------------------------------!
          do lm1 = 1, lmmaxlo(is)
            ilo1 = lmloidx(1,lm1,is)
            m1   = lmloidx(3,lm1,is)            
            do lm2 = 1, lmmaxlo(is)
              ilo2 = lmloidx(1,lm2,is)
              m2   = lmloidx(3,lm2,is)
              pmat(n,m,i) = pmat(n,m,i)+ &
              &             conjg(locmt(m,ilo1,m1,ias))* &
              &             ulodulo(ilo1,m1,ilo2,m2,ias,i)* &
              &             locmt(n,ilo2,m2,ias)
            end do ! lm2
          end do ! lm1
            
        ! end case of local orbitals
        end if
            
      end do ! i
      end do ! n
      end do ! m
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif            
    end do ! ia
    end do ! is
    
    ! multiply y-component with imaginary unit
    pmat(:,:,2) = zi*pmat(:,:,2)
    
    !=============================
    ! PW -- PW
    !=============================
    allocate(cfunt(ngp,ngp))
    allocate(zm(ngp,ndim))
    allocate(pm(mdim,ndim))
    do i = 1, 3
   
      do igp1 = 1, ngp
        ig1 = igpig(igp1)
        do igp2 = 1, ngp
          ig2 = igpig(igp2)
          iv(:) = ivg(:,ig1)-ivg(:,ig2)
          ig = ivgig(iv(1),iv(2),iv(3))
          cfunt(igp1,igp2) = zi*vgpc(i,igp2)*cfunig(ig)
        end do
      end do
      
      ! M_{G'n} = I_{G'G} \times C_{Gn}
      call zgemm ('n', 'n', ngp, ndim, ngp, &
      &           zone, cfunt, ngp, &
      &           evec(1:ngp,nstart:nend), ngp, &
      &           zzero, zm, ngp)
      
      ! [C_{G'm}]^* \times M_{G'n}
      call zgemm ('c', 'n', mdim, ndim, ngp, &
      &           zone, evec(1:ngp,mstart:mend), ngp, &
      &           zm, ngp, &
      &           zzero, pm, mdim)
      
      ! note that pm is transposed
      do m = mstart, mend
      do n = nstart, nend
        pmat(n,m,i) = pmat(n,m,i)+pm(m-mstart+1,n-nstart+1)
      end do
      end do
      
    end do ! i
    deallocate(cfunt,zm,pm)
        
    !================
    ! multiply by -i
    !================
    pmat(:,:,:) = -zi*pmat(:,:,:)
            
    return
  end subroutine

!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
  subroutine genpmatcv_k(vpl,cstart,cend,mstart,mend,pmat)
    use modinput
    use modmain, only : nspecies, natoms, idxas, pi, spl, idxlm, lorbl, &
    &                   zi, nlorb, apword, nlotot, lmmaxapw
    use mod_core_states, only : corind
    implicit none
    ! arguments
    real(8), intent(in) :: vpl(3)
    integer, intent(in) :: cstart, cend
    integer, intent(in) :: mstart, mend
    complex(8), intent(out) :: pmat(cstart:cend,mstart:mend,3)
    ! local variables
    integer :: is, ia, ias
    integer :: l1, m1, lm1, l2, m2, lm2, io2, ilo2, ic
    integer :: i, icg, m
    real(8) :: arg
    complex(8) :: phs
        
    pmat(:,:,:) = 0.d0
        
    !==================
    ! MT -- MT
    !==================
    do icg = cstart, cend
      is = corind(icg,1)
      ia = corind(icg,2)
      ias = idxas(ia,is)
      arg = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(1)*vpl(1)+ &
      &     input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(2)*vpl(2)+ &
      &     input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(3)*vpl(3)
      phs = cmplx(cos(2.0d0*pi*arg),sin(2.0d0*pi*arg),8)
      ic = corind(icg,3)
      l1 = corind(icg,4)
      m1 = corind(icg,5)
      lm1 = idxlm(l1,m1)
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,i,lm2,l2,io2,ilo2,m2)
!$OMP DO
#endif   
      do m = mstart, mend
      do i = 1, 3
        !-----------------------------------!
        !     APW-core-orbital contribution !
        !-----------------------------------!
        do lm2 = 1, lmmaxapw 
          l2 = lmapwidx(1,lm2)
          do io2 = 1, apword(l2,is)
            pmat(icg,m,i) = pmat(icg,m,i)+ &
            &                 conjg(apwcmt(m,io2,lm2,ias))* &
            &                 uaduc(io2,lm2,ic,m1,ias,i)*phs
          end do ! io2
        end do ! lm2
        if (nlotot>0) then
          !---------------------------------------------!
          !    local-orbital-core-orbital contribution  !
          !---------------------------------------------!
          do lm2 = 1, lmmaxlo(is)
            ilo2 = lmloidx(1,lm2,is)
            m2  = lmloidx(3,lm2,is)
            pmat(icg,m,i) = pmat(icg,m,i)+ &
            &                 conjg(locmt(m,ilo2,m2,ias))* &
            &                 uloduc(ilo2,m2,ic,m1,ias,i)*phs
          end do ! lm2
        end if
      end do ! i
      end do ! m
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif          
    end do ! icg
        
    !------------------------------------------
    ! multiply y-component with imaginary unit
    !------------------------------------------
    pmat(:,:,2) = zi*pmat(:,:,2)
        
    !------------------
    ! multiply by -i
    !------------------
    pmat(:,:,:) = -zi*pmat(:,:,:)
        
    return
  end subroutine

!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
  subroutine symmetrize(nstart,nend,mstart,mend,pmat)
    use modmain, only : nsymcrys, lsplsymc, symlatc, symlatd
    implicit none
    integer,    intent(in)    :: nstart, nend
    integer,    intent(in)    :: mstart, mend
    complex(8), intent(inout) :: pmat(nstart:nend,mstart:mend,3)
    ! local
    integer    :: n, m, isym, lspl
    real(8)    :: v1(3), v2(3), sc(3,3)
    complex(8) :: p(3), o(6)

    ! symmetrize pmat
    do n = nstart, nend
    do m = mstart, mend

      v1(:) = dble(pmat(n,m,:))
      call symvect(.false., v1)
      v2(:) = aimag(pmat(n,m,:))
      call symvect(.false., v2)

      pmat(n,m,:) = cmplx(v1,v2,8)

    end do
    end do

    return
  end subroutine

end module
