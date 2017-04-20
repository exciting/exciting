
MODULE mod_kpointset
    
    implicit none

!-------------------------------------------------------------------------------    
    type k_set
        ! k-points
        integer :: nkpt                     ! number of k-points
        integer :: nkptnr                   ! total (non-reduced) number of k-points (in case symmetry is used)
        integer :: ngridk(3)
        real(8), allocatable :: vkl(:,:)    ! lattice coordinates
        real(8), allocatable :: vkc(:,:)    ! cartesian coordinates
        integer, allocatable :: ik2ikp(:)   ! index mapping: non-reduced->reduced k-point
        integer, allocatable :: ikp2ik(:)   ! index mapping: reduced->non-reduced k-point
        real(8), allocatable :: wkpt(:)     ! weight of k-point
        ! tetrahedron integration method related data
        integer :: ntet                     ! number of tetrahedra
        integer, allocatable :: tnodes(:,:) ! coordinates of tetrahedron
        integer, allocatable :: wtet(:)     ! weight of each tetrahedron
        real(8) :: tvol                     ! volume of the tetrahedra relative to the BZ volume
    end type k_set

!-------------------------------------------------------------------------------
    type G_set
        ! grid parameters
        integer :: ngrtot     ! total number of grid points
        integer :: intgv(3,2) ! integer grid size
        ! G-points
        integer :: ngvec                ! number of G-points
        real(8), allocatable :: vgc(:,:)
        real(8), allocatable :: gc(:)
        integer, allocatable :: ivg(:,:)  
        integer, allocatable :: ivgig(:,:,:)
    end type G_set

!-------------------------------------------------------------------------------
    type Gk_set
        ! G+k points
        integer :: ngkmax
        integer, allocatable :: ngk(:,:)      ! number of G+k-vectors for augmented plane waves
        integer, allocatable :: igkig(:,:,:)  ! index from G+k-vectors to G-vectors
        integer, allocatable :: igigk(:,:,:)  ! index from G-vectors to G+k-vectors
        real(8), allocatable :: vgkl(:,:,:,:) ! G+k-vectors in lattice coordinates
        real(8), allocatable :: vgkc(:,:,:,:) ! G+k-vectors in Cartesian coordinates
        real(8), allocatable :: gkc(:,:,:)    ! length of G+k-vectors
        real(8), allocatable :: tpgkc(:,:,:,:)! (theta, phi) coordinates of G+k-vectors
        complex(8), allocatable :: sfacgk(:,:,:,:) ! structure factor for the G+k-vectors
    end type Gk_set
    
!-------------------------------------------------------------------------------    
    type kq_set
        ! k-points
        integer :: nkpt                     ! number of k/q-points
        real(8), allocatable :: vkl(:,:)    ! lattice coordinates
        real(8), allocatable :: vkc(:,:)    ! cartesian coordinates
        real(8), allocatable :: vql(:,:)    ! lattice coordinates
        real(8), allocatable :: vqc(:,:)    ! cartesian coordinates
        ! tetrahedron integration method related data
        integer :: ntet                     ! number of tetrahedra
        integer, allocatable :: tnodes(:,:) ! coordinates of tetrahedron
        integer, allocatable :: wtet(:)     ! weight of each tetrahedron
        integer, allocatable :: kqid(:,:)   ! k-dependent weight of each q-point
        integer, allocatable :: linkq(:,:)  ! number of the tetrahedra linked to by the corresponding q vector
        real(8) :: tvol                     ! volume of the tetrahedra relative to the BZ volume
        ! small group of q-vector (symmetry feature)
        integer, allocatable :: nsymq(:)    ! number of the symmetry operations in the small group of q
        integer, allocatable :: nkptq(:)    ! number of k-points in IBZ(q)
        real(8), allocatable :: wkptq(:,:)  ! q-dependent k-point weight
        integer, allocatable :: iksymq(:,:) ! index of the symmetry operation which rotates the k-point into equivalent one
        integer, allocatable :: ik2ikpq(:,:) ! map the k-point index to the corresponding irreducible one
        integer, allocatable :: ikp2ikq(:,:) ! map the irreducible k-point index to the corresponding from the non-reduced set
        integer, allocatable :: nsymkstar(:,:)   ! number of symmetry operations which form the star
        integer, allocatable :: isymkstar(:,:,:) ! index of the symmetry operations
    end type kq_set

    ! internal routine, should not be visible outside
    private  rtorat
    
CONTAINS

!-------------------------------------------------------------------------------
    subroutine rtorat(x, k, div)
        !-------------------------------------------------------
        ! Factorizes the real coordinates of a vector {\bf x}.
        ! The output is an integer vector {\bf k}, such that
        ! $$ |x(i)-k(i)/{\rm div}| < {\rm eps} $$
        ! for all $i=1,\ldots,n$.
        !   
        ! Created July 2008 by Sagmeister
        !-------------------------------------------------------    
        implicit none
        real(8),    intent(In)  :: x(3)
        integer(4), intent(Out) :: k(3)
        integer(4), intent(Out) :: div
        integer :: maxint
        real(8) :: dx
        real(8) :: eps = 1.0d-5
        maxint = nint(1.d0/eps)/10
        do div = 1, maxint
          k(:) = nint(dble(div)*x(:))
          dx = maxval(Abs(dble(k)/dble(div)-x))
          if (dx < eps) exit
        end do
        if (dx >= eps) then
          write(*,*)
          write(*, '("Error(rtorat): factorization failed")')
          write(*, '(" maximum integer :",i12)') maxint
          write(*, '(" tolerance       :",g18.10)') eps
          write(*, '(" deviation       :",g18.10)') dx
          write(*,*)
          stop
        end if
        if (dx > 1.d-12) then
          write(*,*)
          write(*, '("Warning(rtorat): small deviation in factorization")')
          write(*, '(" maximum deviation :",g18.10)') dx
          write(*,*)
        end if
    end subroutine rtorat

!-------------------------------------------------------------------------------
    subroutine generate_k_vectors(self,bvec,ngridk,vkloff,reduce)
        use modmain, only: nsymcrys, symlat, lsplsymc
        implicit none
        type(k_set), intent(OUT) :: self
        real(8), intent(IN) :: bvec(3,3)
        integer, intent(IN) :: ngridk(3)
        real(8), intent(IN) :: vkloff(3)
        logical, intent(IN) :: reduce ! apply symmetry to reduce k-set     
        ! local variables
        integer(4) :: dvk
        integer(4) :: mnd
        integer(4) :: ikloff(3)
        integer(4) :: dkloff
        integer :: i1, i2, ik, nsym, isym, lspl
        integer(4), allocatable :: symmat(:,:,:)
        integer(4), allocatable :: ivk(:,:)
        integer(4), allocatable :: iwkp(:)
        !
        ! initialize k-set
        self%nkptnr = ngridk(1)*ngridk(2)*ngridk(3)
        self%nkpt = self%nkptnr
        self%ngridk = ngridk
        if (allocated(self%vkl)) deallocate(self%vkl)
        allocate(self%vkl(3,self%nkpt))
        if (allocated(self%vkc)) deallocate(self%vkc)
        allocate(self%vkc(3,self%nkpt))
        if (allocated(self%ik2ikp)) deallocate(self%ik2ikp)
        allocate(self%ik2ikp(self%nkpt))
        if (allocated(self%ikp2ik)) deallocate(self%ikp2ik)
        allocate(self%ikp2ik(self%nkpt))
        if (allocated(self%wkpt)) deallocate(self%wkpt)
        allocate(self%wkpt(self%nkpt))
        ! tetrahedron integration method related data
        self%ntet = 6*self%nkpt
        if (allocated(self%tnodes)) deallocate(self%tnodes)
        allocate(self%tnodes(4,self%ntet))
        if (allocated(self%wtet)) deallocate(self%wtet)
        allocate(self%wtet(self%ntet))
        self%tvol = 0.d0
        
        ! explore symmetry
        nsym = 1
        if (reduce) nsym = nsymcrys
        allocate(symmat(3,3,nsym))
        do isym = 1, nsym
          lspl = lsplsymc(isym)
          ! the library used transposed rotational matrices
          do i1 = 1, 3
            do i2 = 1, 3
              symmat(i1,i2,isym) = symlat(i2,i1,lspl)
            end do
          end do
        end do
        
        ! k-mesh shift
        !call factorize(3,vkloff,ikloff,dkloff) !<-- Libbzint routine
        call rtorat(vkloff,ikloff,dkloff)
        
        ! call LibBZint library
        allocate(ivk(3,self%nkpt))
        allocate(iwkp(self%nkpt))
        
        ! modified (by DIN) version of libbzint routine kgen
        ! I have added an additional output parameter mapping ikp -> ik 
        call kgen_exciting(bvec,nsym,symmat,ngridk,ikloff,dkloff, &
        &         self%nkpt,ivk,dvk,self%ik2ikp,self%ikp2ik,iwkp, &
        &         self%ntet,self%tnodes,self%wtet,self%tvol,mnd)
        
        ! fractional and cartesian coordinates, and k-point weight
        do ik = 1, self%nkpt
            self%vkl(:,ik) = dble(ivk(:,ik))/dble(dvk)
            call r3mv(bvec,self%vkl(:,ik),self%vkc(:,ik))
            self%wkpt(ik) = dble(iwkp(ik))/dble(self%nkptnr)
        enddo ! ik
        
        deallocate(symmat,ivk,iwkp)
        return
    end subroutine
    
!-------------------------------------------------------------------------------
    subroutine delete_k_vectors(self)
        type(k_set), intent(INOUT) :: self
        if (allocated(self%vkl)) deallocate(self%vkl)
        if (allocated(self%vkc)) deallocate(self%vkc)
        if (allocated(self%ik2ikp)) deallocate(self%ik2ikp)
        if (allocated(self%ikp2ik)) deallocate(self%ikp2ik)
        if (allocated(self%wkpt)) deallocate(self%wkpt)
        if (allocated(self%tnodes)) deallocate(self%tnodes)
        if (allocated(self%wtet)) deallocate(self%wtet)
    end subroutine

!-------------------------------------------------------------------------------
    subroutine print_k_vectors(self,funit)
        implicit none
        type(k_set), intent(IN) :: self
        integer,     intent(IN) :: funit
        integer :: ik, i, j
        !
        call boxmsg(funit,'-','k-vectors')
        write(funit,*) 'Total number of k-points: < nkptnr >', self%nkptnr
        write(funit,*) 'Symmetry reduced number of k-points: < nkpt >', self%nkpt
        write(funit,*) 'k-vectors list: < ik    vkl    vkc    weight >'
        do ik = 1, self%nkpt
          write(funit,100) ik, self%vkl(1:3,ik), self%vkc(1:3,ik), self%wkpt(ik)
        enddo
        100 format(i4,4x,3f8.4,4x,3f8.4,4x,f8.4)
        write(funit,*) 'Mapping from non-reduced to reduced grid: < ik2ikp >'
        write(funit,*) (self%ik2ikp(ik), ik=1,self%nkptnr)
        write(funit,*) 'Mapping from reduced to non-reduced grid: < ikp2ik >'
        write(funit,*) (self%ikp2ik(ik), ik=1,self%nkpt)
        write(funit,*)
        write(funit,*) 'Tetrahedron method info: < ntet    tvol >'
        write(funit,101) self%ntet, self%tvol
        write(funit,*) 'Nodal points of tetrahedron: < itet    tnodes    wtet >'
        do i = 1, self%ntet
          write(funit,*) i, (self%tnodes(j,i),j=1,4), self%wtet(i)
        enddo 
        101 format(i6,e16.8)
        return
    end subroutine

!-------------------------------------------------------------------------------
    subroutine generate_G_vectors(self,bvec,intgv,gmaxvr)
        implicit none
        type(G_set), intent(OUT) :: self
        real(8), intent(IN) :: bvec(3,3)
        integer, intent(IN) :: intgv(3,2) ! integer ranges for G-grid
        real(8), intent(IN) :: gmaxvr
        ! local variables
        integer :: ig, i1, i2, i3, k
        real(8) :: v(3), t1
        ! allocatable arrays
        integer, allocatable :: idx(:)
        integer, allocatable :: iar(:)
        real(8), allocatable :: rar(:)
        ! initialize G-set
        self%intgv(:,:) = intgv(:,:)
        self%ngrtot = (self%intgv(1,2)-self%intgv(1,1)+1)* &
        &             (self%intgv(2,2)-self%intgv(2,1)+1)* &
        &             (self%intgv(3,2)-self%intgv(3,1)+1)
        self%ngvec = 0
        if (allocated(self%ivg)) deallocate(self%ivg)
        allocate(self%ivg(3,self%ngrtot))
        if (allocated(self%ivgig)) deallocate(self%ivgig)
        allocate(self%ivgig(intgv(1,1):intgv(1,2), &
        &                   intgv(2,1):intgv(2,2), &
        &                   intgv(3,1):intgv(3,2)))
        if (allocated(self%vgc)) deallocate(self%vgc)
        allocate(self%vgc(3,self%ngrtot))
        if (allocated(self%gc)) deallocate(self%gc)
        allocate(self%gc(self%ngrtot))
        ! allocate local arrays
        allocate(idx(self%ngrtot))
        allocate(iar(self%ngrtot))
        allocate(rar(self%ngrtot))
        ig = 0
        do i1 = intgv(1,1), intgv(1,2)
          do i2 = intgv(2,1), intgv(2,2)
            do i3 = intgv(3,1), intgv(3,2)
              v(:) = dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2)+dble(i3)*bvec(:,3)
              t1 = v(1)**2+v(2)**2+v(3)**2
              ig = ig+1
              ! map from G-vector to (i1,i2,i3) index
              self%ivg(1,ig) = i1
              self%ivg(2,ig) = i2
              self%ivg(3,ig) = i3
              ! length of each G-vector
              self%gc(ig) = sqrt(t1)
            end do
          end do
        end do
        ! sort by vector length
        call sortidx(self%ngrtot,self%gc,idx)
        ! re-order arrays
        do ig = 1, self%ngrtot
          rar(ig) = self%gc(ig)
        end do
        do ig = 1, self%ngrtot
         self%gc(ig) = rar(idx(ig))
        end do
        do k = 1, 3
          do ig = 1, self%ngrtot
            iar(ig) = self%ivg(k,ig)
          end do
          do ig = 1, self%ngrtot
            self%ivg(k,ig) = iar(idx(ig))
          end do
        end do
        self%ivgig(:,:,:) = 0
        do ig = 1, self%ngrtot
          i1 = self%ivg(1,ig)
          i2 = self%ivg(2,ig)
          i3 = self%ivg(3,ig)
          ! map from (i1,i2,i3) index to G-vector
          self%ivgig(i1,i2,i3) = ig
          ! assign G-vectors to global array
          self%vgc(:,ig) = dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2)+dble(i3)*bvec(:,3)
        end do
        ! find the number of vectors with G < gmaxvr
        self%ngvec = 1
        do ig = self%ngrtot, 1, -1
          if (self%gc(ig) < gmaxvr) then
            self%ngvec = ig
            exit
          end if
        end do
        deallocate(idx,iar,rar)
        return
    end subroutine
    
!-------------------------------------------------------------------------------
    subroutine delete_G_vectors(self)
        type(G_set), intent(INOUT) :: self
        if (allocated(self%vgc)) deallocate(self%vgc)
        if (allocated(self%gc)) deallocate(self%gc)
        if (allocated(self%ivg)) deallocate(self%ivg)
        if (allocated(self%ivgig)) deallocate(self%ivgig)
    end subroutine

!-------------------------------------------------------------------------------
    subroutine print_G_vectors(self,funit)
        implicit none
        type(G_set), intent(IN) :: self
        integer,     intent(IN) :: funit
        integer :: ig, i1, i2, i3
        call boxmsg(funit,'-','G-vectors')
        write(funit,*) 'Grid size: < b-axis  intgv(start) intgv(end) >'
        do i1 = 1, 3
          write(funit,*) i1, self%intgv(i1,:)
        end do ! i1
        write(funit,*) 'Total number of grid points: < ngrtot >', self%ngrtot
        write(funit,*) 'Number of G-vectors: < ngvec >', self%ngvec
        write(funit,*) 'List of G-vectors:   < ig    i1,i2,i3    ivgig    vgc    gc >'
        do ig = 1, self%ngvec
          i1 = self%ivg(1,ig)
          i2 = self%ivg(2,ig)
          i3 = self%ivg(3,ig)
          write(funit,102) ig, i1, i2, i3, self%ivgig(i1,i2,i3), &
          &                self%vgc(1:3,ig), self%gc(ig)
        end do
        102 format(i6,4x,3i4,4x,i6,4x,4f10.4)
        return
    end subroutine

!-------------------------------------------------------------------------------
    subroutine generate_Gk_vectors(self,kset,Gset,gkmax)
        use modmain, only: nspnfv, natmtot
        implicit none
        type(Gk_set), intent(OUT) :: self
        type(k_set),  intent(IN)  :: kset
        type(G_set),  intent(IN)  :: Gset
        real(8),      intent(IN)  :: gkmax
        ! local variables
        integer :: ispn, ik, ig, igp
        real(8) :: v(3), t1
        integer, allocatable :: igk2ig(:,:,:)
        ! 
        allocate(igk2ig(Gset%ngrtot,kset%nkpt,nspnfv))
        igk2ig(:,:,:) = 0
        !
        ! determine the number of G+k combinations which satisfy |G+k|<gkmax
        if (allocated(self%ngk)) deallocate(self%ngk)
        allocate(self%ngk(nspnfv,kset%nkpt))
        do ispn = 1, nspnfv
          do ik = 1, kset%nkpt
            igp = 0
            do ig = 1, Gset%ngrtot
              v(:) = Gset%vgc(:,ig)+kset%vkc(:,ik)
              t1 = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
              if (t1 < gkmax) then
                igp = igp+1
                igk2ig(igp,ik,ispn) = ig
              end if
            end do ! ig
            self%ngk(ispn,ik) = igp
          end do ! ik
        end do ! ispn
        !
        ! maximum number of G+k vectors
        self%ngkmax = maxval(self%ngk)
        if (self%ngkmax > Gset%ngrtot) then
          write(*,*) 'ERROR(mod_kpoints::generate_Gk_vectors) ngkmax > ngrtot'
          stop
        end if
        !
        ! generate G+k data set
        if (allocated(self%igkig)) deallocate(self%igkig)
        allocate(self%igkig(self%ngkmax,nspnfv,kset%nkpt))
        if (allocated(self%vgkl)) deallocate(self%vgkl)
        allocate(self%vgkl(3,self%ngkmax,nspnfv,kset%nkpt))
        if (allocated(self%vgkc)) deallocate(self%vgkc)
        allocate(self%vgkc(3,self%ngkmax,nspnfv,kset%nkpt))
        if (allocated(self%gkc)) deallocate(self%gkc)
        allocate(self%gkc(self%ngkmax,nspnfv,kset%nkpt))
        if (allocated(self%tpgkc)) deallocate(self%tpgkc)
        allocate(self%tpgkc(2,self%ngkmax,nspnfv,kset%nkpt))
        if (allocated(self%sfacgk)) deallocate(self%sfacgk)
        allocate(self%sfacgk(self%ngkmax,natmtot,nspnfv,kset%nkpt))
        do ispn = 1, nspnfv
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,igp,ig)
!$OMP DO
#endif    
          do ik = 1, kset%nkpt
            do igp = 1, self%ngk(ispn,ik)
              ig = igk2ig(igp,ik,ispn)
              ! index to G-vector
              self%igkig(igp,ispn,ik) = ig
              ! G+k-vector in lattice coordinates
              self%vgkl(:,igp,ispn,ik) = dble(Gset%ivg(:,ig))+kset%vkl(:,ik)
              ! G+k-vector in Cartesian coordinates
              self%vgkc(:,igp,ispn,ik) = Gset%vgc(:,ig)+kset%vkc(:,ik)
              ! G+k-vector length and (theta, phi) coordinates
              call sphcrd(self%vgkc(:,igp,ispn,ik),self%gkc(igp,ispn,ik),&
              &           self%tpgkc(:,igp,ispn,ik))
              ! generate structure factors for G+k-vectors
              call gensfacgp(self%ngk(ispn,ik),self%vgkc(:,:,ispn,ik), &
              &              self%ngkmax,self%sfacgk(:,:,ispn,ik))
            end do
          end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
        end do
        deallocate(igk2ig)
        return
    end subroutine
    
!-------------------------------------------------------------------------------
    subroutine delete_Gk_vectors(self)
        type(Gk_set), intent(INOUT) :: self
        if (allocated(self%ngk)) deallocate(self%ngk)
        if (allocated(self%igkig)) deallocate(self%igkig)
        if (allocated(self%vgkl)) deallocate(self%vgkl)
        if (allocated(self%vgkc)) deallocate(self%vgkc)
        if (allocated(self%gkc)) deallocate(self%gkc)
        if (allocated(self%tpgkc)) deallocate(self%tpgkc)
        if (allocated(self%sfacgk)) deallocate(self%sfacgk)
    end subroutine

!-------------------------------------------------------------------------------
    subroutine print_Gk_vectors(self,ik,funit)
        implicit none
        type(Gk_set), intent(IN) :: self
        integer,      intent(IN) :: ik
        integer,      intent(IN) :: funit
        integer :: igp
        call boxmsg(funit,'-','G+k-vectors')
        write(funit,*) 'Maximal number of G+k-vectors: < ngkmax >', self%ngkmax
        write(funit,*) 'List of G+k-vectors for k-point:', ik
        write(funit,*) ' * Number of G+k combinations: < ngk >', self%ngk(1,ik)
        write(funit,*) ' * < igp    igkig    vgkl    vgkc    gkc    tpgkc    sfacgk >'
        do igp = 1, self%ngk(1,ik)
           write(funit,103) igp, self%igkig(igp,1,ik), &
           &                self%vgkl(1:3,igp,1,ik), self%vgkc(1:3,igp,1,ik), &
           &                self%gkc(igp,1,ik), self%tpgkc(1:2,igp,1,ik), &
           &                self%sfacgk(igp,1,1,ik)
        end do
        103 format(i6,2x,i6,4x,3f8.4,4x,3f8.4,4x,f8.4,4x,2f8.4,4x,2f8.4)
        return
    end subroutine

!-------------------------------------------------------------------------------
    subroutine generate_kq_vectors(self,bvec,ngridk,vkloff,reduce)
        implicit none
        type(kq_set), intent(OUT) :: self
        real(8), intent(IN) :: bvec(3,3)
        integer, intent(IN) :: ngridk(3)
        real(8), intent(IN) :: vkloff(3)
        logical, intent(IN) :: reduce
        ! local variables        
        integer(4) :: ikloff(3)
        integer(4) :: dkloff
        integer(4) :: ik
        integer(4) :: dvk, dvq
        integer(4), allocatable :: ivk(:,:), ivq(:,:)
        !
        ! initialize k/q-sets
        self%nkpt = ngridk(1)*ngridk(2)*ngridk(3)
        if (allocated(self%vkl)) deallocate(self%vkl)
        allocate(self%vkl(3,self%nkpt))
        if (allocated(self%vkc)) deallocate(self%vkc)
        allocate(self%vkc(3,self%nkpt))
        if (allocated(self%vql)) deallocate(self%vql)
        allocate(self%vql(3,self%nkpt))
        if (allocated(self%vqc)) deallocate(self%vqc)
        allocate(self%vqc(3,self%nkpt))
        
        ! tetrahedron integration method related data
        self%ntet = 6*self%nkpt
        if (allocated(self%linkq)) deallocate(self%linkq)
        allocate(self%linkq(self%ntet,self%nkpt))
        if (allocated(self%kqid)) deallocate(self%kqid)
        allocate(self%kqid(self%nkpt,self%nkpt))
        if (allocated(self%tnodes)) deallocate(self%tnodes)
        allocate(self%tnodes(4,self%ntet))
        if (allocated(self%wtet)) deallocate(self%wtet)
        allocate(self%wtet(self%ntet))
        self%tvol = 0.d0
        
        ! k-mesh shift
        !call factorize(3,vkloff,ikloff,dkloff) !<-- Libbzint routine
        call rtorat(vkloff,ikloff,dkloff)
        
        ! call LibBZint library
        allocate(ivk(3,self%nkpt))
        allocate(ivq(3,self%nkpt))
        !
        
        ! suppress debug output in tetrahedron integration library (0)
        !call tetrasetdbglv(0)
        ! Generate the k- and q-points meshes
        !call kqgen(bvec,ngridk,ikloff,dkloff, &
        !&          self%nkpt,ivk,ivq,dvk,dvq,self%kqid,self%ntet, &
        !&          self%tnodes,self%wtet,self%linkq,self%tvol)
        
        call kqgen(bvec,ngridk,ikloff,dkloff, &
        &          self%nkpt,ivk,ivq,dvk,dvq,self%kqid,self%ntet, &
        &          self%tnodes,self%wtet,self%linkq,self%tvol)

        !      
        ! fractional and cartesian coordinates
        do ik = 1, self%nkpt
          self%vkl(:,ik) = dble(ivk(:,ik))/dble(dvk)
          call r3mv(bvec,self%vkl(:,ik),self%vkc(:,ik))
          self%vql(:,ik) = dble(ivq(:,ik))/dble(dvq)
          call r3mv(bvec,self%vql(:,ik),self%vqc(:,ik))
        end do ! ik
        deallocate(ivk,ivq)
        !
        if (reduce) call generate_small_group_q(self,bvec,ngridk,vkloff,reduce)
        !
        return
    end subroutine

!-------------------------------------------------------------------------------
    subroutine generate_small_group_q(self,bvec,ngridk,vkloff,reduce)
        use modmain, only: nsymcrys, symlat, lsplsymc
        implicit none
        type(kq_set), intent(InOut) :: self
        real(8), intent(IN) :: bvec(3,3)
        integer, intent(IN) :: ngridk(3)
        real(8), intent(IN) :: vkloff(3)        
        logical,      intent(In)    :: reduce
        ! local variables    
        integer :: iqp, ikp  ! q-, k-point indexes for IBZ
        integer :: ik, nqpt, nkptnr
        integer :: ip, jp, iv(3)
        integer :: isym, lspl, nsym
        real(8) :: s(3,3), v1(3), v2(3), t1
        real(8), allocatable :: vklq(:,:)
        integer, allocatable :: scmapq(:)    
        integer, allocatable :: ivwrapq(:,:)
        integer, allocatable :: iwkpq(:)
        type(k_set) :: qpt
        real(8), external :: r3taxi
        real(8), parameter :: epslat=1.d-6
        ! total non-reduced number of k-points (shortcut)
        nkptnr = self%nkpt
        ! generate q-point set
        call generate_k_vectors(qpt,bvec,ngridk,vkloff,reduce)        
        ! if symmetry is used        
        if (reduce) then
          nsym = nsymcrys
        else
          nsym = 1
        end if
        ! resulting number of q-points
        nqpt = qpt%nkpt
        ! allocate global data arrays
        if (allocated(self%nsymq)) deallocate(self%nsymq)
        allocate(self%nsymq(nqpt))
        self%nsymq(:) = 0
        if (allocated(self%nkptq)) deallocate(self%nkptq)
        allocate(self%nkptq(nqpt))
        self%nkptq(:) = 0
        if (allocated(self%wkptq)) deallocate(self%wkptq)
        allocate(self%wkptq(nkptnr,nqpt))
        self%wkptq(:,:) = 0.d0
        if (allocated(self%ik2ikpq)) deallocate(self%ik2ikpq)
        allocate(self%ik2ikpq(nkptnr,nqpt))
        self%ik2ikpq(:,:) = 0
        if (allocated(self%ikp2ikq)) deallocate(self%ikp2ikq)
        allocate(self%ikp2ikq(nkptnr,nqpt))
        self%ikp2ikq(:,:) = 0
        if (allocated(self%iksymq)) deallocate(self%iksymq)
        allocate(self%iksymq(nkptnr,nqpt))
        self%iksymq(:,:) = 0
        if (allocated(self%nsymkstar)) deallocate(self%nsymkstar)
        allocate(self%nsymkstar(nkptnr,nqpt))
        self%nsymkstar(:,:) = 0
        if (allocated(self%isymkstar)) deallocate(self%isymkstar)
        allocate(self%isymkstar(nsym,nkptnr,nqpt))
        self%isymkstar(:,:,:) = 0
        !
        ! allocate local arrays
        allocate(vklq(3,nkptnr))
        allocate(scmapq(nsym))
        allocate(ivwrapq(3,nsym))
        allocate(iwkpq(nkptnr))
        !
        ! for each q-point
        do iqp = 1, nqpt
          !---------------------------!
          ! find the small group of q ! 
          !---------------------------!
          v1(:) = qpt%vkl(:,iqp)
          !
          call findgroupq(.false.,v1,epslat, &
          &               bvec,symlat,nsym,lsplsymc, &
          &               self%nsymq(iqp),scmapq,ivwrapq)
          !
          !--------------------------------------------------------------------------------
          ! determine the q-dependent BZ(q) and symmetry operations which regenerate
          ! the full (non-reduced) k-grid
          !--------------------------------------------------------------------------------      
          ip = 0
          do ik = 1, nkptnr
            v1(:) = self%vkl(:,ik)
            ! determine if this point is equivalent to that already in the set
            do isym = 1, self%nsymq(iqp)
              lspl = lsplsymc(scmapq(isym))
              s(:,:) = dble(symlat(:,:,lspl))
              call r3mtv(s,v1,v2)
              call r3frac(epslat,v2,iv)
              do jp = 1, ip
                t1 = r3taxi(vklq(:,jp),v2)
                if (t1 < epslat) then
                  ! equivalent k-point found so add to current weight
                  self%ik2ikpq(ik,iqp) = jp 
                  iwkpq(jp) = iwkpq(jp) + 1
                  self%iksymq(ik,iqp) = scmapq(isym)
                  goto 10
                end if
              end do !jp
            end do ! isym
            ! add new point to set
            ip = ip + 1
            vklq(:,ip) = v1(:)
            self%ik2ikpq(ik,iqp) = ip 
            iwkpq(ip) = 1
            self%iksymq(ik,iqp) = 1
            10 continue
          end do ! ik
          !
          ! save the number N(q) of k-points            
          self%nkptq(iqp) = ip 
          !
          ! q-dependent k-point weight      
          do ik = 1, self%nkptq(iqp)
            self%wkptq(ik,iqp)=dble(iwkpq(ik))/dble(nkptnr)
          end do
          !
          ! Determine the index of the irreducible point in the non-reduced set
          do ik = 1, nkptnr
            v1(:) = self%vkl(:,ik)
            v2(:) = vklq(:,self%ik2ikpq(ik,iqp))
            t1 = r3taxi(v1,v2)
            if (t1 < epslat) then
              self%ikp2ikq(self%ik2ikpq(ik,iqp),iqp) = ik
            end if
          end do
          !
          ! Determine the star G_k
          do ik = 1, nkptnr
            ikp = self%ik2ikpq(ik,iqp)
            self%nsymkstar(ikp,iqp) = self%nsymkstar(ikp,iqp)+1
            self%isymkstar(self%nsymkstar(ikp,iqp),ikp,iqp) = self%iksymq(ik,iqp)
          end do ! ik
          !
        end do ! iqp
        !
        call delete_k_vectors(qpt)
        deallocate(iwkpq)
        deallocate(vklq)
        deallocate(scmapq)
        deallocate(ivwrapq)
        !
        return
    end subroutine
    
!-------------------------------------------------------------------------------
    subroutine delete_kq_vectors(self)
        type(kq_set), intent(INOUT) :: self
        if (allocated(self%vkl)) deallocate(self%vkl)
        if (allocated(self%vkc)) deallocate(self%vkc)
        if (allocated(self%vql)) deallocate(self%vql)
        if (allocated(self%vqc)) deallocate(self%vqc)
        if (allocated(self%tnodes)) deallocate(self%tnodes)
        if (allocated(self%wtet)) deallocate(self%wtet)
        if (allocated(self%kqid)) deallocate(self%kqid)
        if (allocated(self%linkq)) deallocate(self%linkq)
        if (allocated(self%nsymq)) deallocate(self%nsymq)
        if (allocated(self%nkptq)) deallocate(self%nkptq)
        if (allocated(self%wkptq)) deallocate(self%wkptq)
        if (allocated(self%ik2ikpq)) deallocate(self%ik2ikpq)
        if (allocated(self%ikp2ikq)) deallocate(self%ikp2ikq)
        if (allocated(self%nsymkstar)) deallocate(self%nsymkstar)
        if (allocated(self%isymkstar)) deallocate(self%isymkstar)
    end subroutine

!-----------------------------------------------------------------------------
    subroutine print_kq_vectors(self,funit)
        implicit none
        type(kq_set), intent(IN) :: self
        integer,      intent(IN) :: funit
        integer :: ik, iq, i, j, k, kmax
        call boxmsg(funit,'-','k/q-vectors')
        write(funit,*)
        write(funit,*) 'Total number of k/q-points:', self%nkpt
        write(funit,*) 'k-vectors list:  < ik    vkl    vkc    weight >'
        do ik = 1, self%nkpt
          write(funit,104) ik, self%vkl(1:3,ik), self%vkc(1:3,ik), 1.d0/dble(self%nkpt)
        enddo
        write(funit,*) 'q-vectors list:  < iq    vql    vqc    weight >'
        do ik = 1, self%nkpt
          write(funit,104) ik, self%vql(1:3,ik), self%vqc(1:3,ik), 1.d0/dble(self%nkpt)
        enddo        
        104 format(i4,4x,3f8.4,4x,3f8.4,4x,f8.4)
        write(funit,*)
        write(funit,*) 'Tetrahedron method info: < ntet    tvol >'
        write(funit,105) self%ntet, self%tvol
        write(funit,*) 'Nodal points of tetrahedron: < ik    tnodes    wtet >'
        do i = 1, self%ntet
          write(funit,*) i, (self%tnodes(j,i),j=1,4), self%wtet(i)
        enddo 
        105 format(i6,e16.8)
        write(funit,*)
        write(funit,*) 'Id of tetrahedra linked by the corresponding q vector: < iq    linkq >'
        do i = 1, self%ntet
          write(funit,*) i, (self%linkq(i,j), j=1,self%nkpt)
        enddo 
        write(funit,*)
        write(funit,*) 'Link between k and k-q points: < kqid >'
        do i = 1, self%nkpt, 8
          kmax = 7
          if ((self%nkpt-i) < kmax) kmax = self%nkpt-i
          write(funit,106) (i+k, k=0,kmax)
          do j = 1, self%nkpt
            write(funit,107) j,(self%kqid(j,i+k), k=0,kmax)
          enddo
        enddo
        106 format(' ik/iq |',8i4)      
        107 format(i4,'   |',8i4)      
        !
        ! Print small group of q-vectors (when symmetry is used)
        if (allocated(self%nkptq)) then
          call boxmsg(funit,'-','Generate a small group of q-vectors')
          do iq = 1, size(self%nkptq)
            write(funit,*)
            write(funit,'("q-point: ",i4,4x,3f8.4)') iq, self%vkl(:,iq)
            write(funit,*)
            write(funit,*) '  nsymq = ', self%nsymq(iq)
            write(funit,*) '  nkptq = ', self%nkptq(iq)
            write(funit,*) '  Small group of k-points: < ik    vklq    wkptq >'
            do ik = 1, self%nkptq(iq)
              write(funit,'(i4,4f12.4)') ik, self%vkl(:,self%ikp2ikq(ik,iq)), self%wkptq(ik,iq)
            end do
            write(funit,*)  
            write(funit,*) '  Mapping arrays:'
            write(funit,*) '    ik2ikpq: ', self%ik2ikpq(:,iq)
            write(funit,*) '    ikp2ikq: ', self%ikp2ikq(1:self%nkptq(iq),iq)      
            write(funit,*) 
            write(funit,*) '  Star operations:'
            write(funit,*) '    iksymq: ', self%iksymq(:,iq)
            do ik = 1, self%nkptq(iq)
              write(funit,*) '    ik =', ik, '    nsymkstar =', self%nsymkstar(ik,iq)
              write(funit,*) '    isymkstar =', self%isymkstar(1:self%nsymkstar(ik,iq),ik,iq)
            enddo
            call linmsg(funit,'-','')
          end do
        end if ! symmetry is used

        return
    end subroutine

END MODULE

