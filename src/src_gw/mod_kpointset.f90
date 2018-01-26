
MODULE mod_kpointset
    
    implicit none

!-------------------------------------------------------------------------------    
    type k_set
        ! k-points
        integer(4) :: nkptnr                ! total (non-reduced) number of k-points (in case symmetry is used)
        integer(4) :: nkpt                  ! number of k-points
        logical :: isreduced                ! Symmetry was used in construction
        logical :: usedlibzint              ! Was build with libzint
        real(8), allocatable :: bvec(:,:)   ! Reciprocal lattice basis
        real(8), allocatable :: vkloff(:)   ! Offset of k-grid in k-coordinates
        integer(4), allocatable :: ngridk(:)! Number of k-points along [0,1) in each lattice direction

        ! Reduced quantities
        integer(4), allocatable :: ivk(:,:) ! 3d integer index of k points
        real(8), allocatable :: vkl(:,:)    ! lattice coordinates
        real(8), allocatable :: vkc(:,:)    ! cartesian coordinates
        real(8), allocatable :: wkpt(:)     ! weight of k-point

        ! Non reduced quantities (created only if no libzint is used)
        integer(4), allocatable :: ivknr(:,:)
        real(8), allocatable :: vklnr(:,:)
        real(8), allocatable :: vkcnr(:,:)
        real(8), allocatable :: wkptnr(:) 

        ! 3D index maps
        integer(4), allocatable ::  ikmap(:,:,:)   ! non-reduced 3d index -> 1d reduced index
        integer(4), allocatable ::  ikmapnr(:,:,:) ! non-reduced 3d index -> 1d non-reduced index

        ! 1D index maps
        integer(4), allocatable :: ik2ikp(:)   ! 1d non-reduced index -> 1d reduced k-point index
        integer(4), allocatable :: ikp2ik(:)   ! 1d reduced index -> 1d non-reduced k-point index

        ! tetrahedron integration method related data
        integer(4) :: ntet                     ! number of tetrahedra
        integer(4), allocatable :: tnodes(:,:) ! coordinates of tetrahedron
        integer(4), allocatable :: wtet(:)     ! weight of each tetrahedron
        real(8) :: tvol                     ! volume of the tetrahedra relative to the BZ volume

    end type k_set

!-------------------------------------------------------------------------------
    type G_set
        ! grid parameters
        real(8) :: gmaxvr     ! Maximum G vector length
        real(8), allocatable :: bvec(:,:)    ! Lattice basis vectors
        integer :: ngrtot     ! total number of grid points
        integer, allocatable :: intgv(:,:) ! integer grid size
        ! G-points
        integer :: ngvec                  ! number of G-points
        real(8), allocatable :: vgc(:,:)  ! Cartesian coodinates of lattice
        real(8), allocatable :: gc(:)     ! Length of G vector
        integer, allocatable :: ivg(:,:)  ! integer coordinates
        integer, allocatable :: ivgig(:,:,:) ! integer coordinates -> 1d index
    end type G_set

!-------------------------------------------------------------------------------
    type Gk_set

        ! Cutoff for G+k length
        real(8) :: gkmax   

        ! Reduced (potentially) quantities
        integer :: ngkmax
        integer, allocatable :: ngk(:,:)      ! number of G+k-vectors for augmented plane waves
        integer, allocatable :: igkig(:,:,:)  ! index from G+k-vectors to G-vectors
        integer, allocatable :: igigk(:,:,:)  ! index from G-vectors to G+k-vectors
        real(8), allocatable :: vgkl(:,:,:,:) ! G+k-vectors in lattice coordinates
        real(8), allocatable :: vgkc(:,:,:,:) ! G+k-vectors in Cartesian coordinates
        real(8), allocatable :: gkc(:,:,:)    ! length of G+k-vectors
        real(8), allocatable :: tpgkc(:,:,:,:)! (theta, phi) coordinates of G+k-vectors
        complex(8), allocatable :: sfacgk(:,:,:,:) ! structure factor for the G+k-vectors

        ! Non-reduced quantities (created if k set was build without libzint)
        integer :: ngknrmax
        integer, allocatable :: ngknr(:,:)      ! number of G+k-vectors for augmented plane waves
        integer, allocatable :: igknrig(:,:,:)  ! index from G+k-vectors to G-vectors
        integer, allocatable :: igigknr(:,:,:)  ! index from G-vectors to G+k-vectors
        real(8), allocatable :: vgknrl(:,:,:,:) ! G+k-vectors in lattice coordinates
        real(8), allocatable :: vgknrc(:,:,:,:) ! G+k-vectors in Cartesian coordinates
        real(8), allocatable :: gknrc(:,:,:)    ! length of G+k-vectors
        real(8), allocatable :: tpgknrc(:,:,:,:)! (theta, phi) coordinates of G+k-vectors
        complex(8), allocatable :: sfacgknr(:,:,:,:) ! structure factor for the G+k-vectors

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

!-------------------------------------------------------------------------------    
    type kkqmt_set
        ! Momentum transfer Q-vector Q = q_Q + G_Q
        real(8), allocatable :: vqmtl(:)      ! Lattice coordinates
        real(8), allocatable :: vqmtc(:)      ! Cartesian coordinates
        !   Unit cell part
        real(8), allocatable :: vqmtl_q(:)    ! Lattice coordinates 
        real(8), allocatable :: vqmtc_q(:)    ! Cartesian coordinates
        !   Lattice vector part
        integer(4), allocatable :: vqmtl_g(:) ! Lattice coordinates
        real(8), allocatable :: vqmtc_g(:)    ! Cartesian coordinates
        ! k-points
        type(k_set) :: kset
        ! k+qmt-points
        type(k_set) :: kqmtset

        !! Index mapping between k-grid and k+qmt-grid
        ! ik --> ik'+ig, where ik' is ik+qmt mapped back to unit cell
        ! and ig the index of the corresponding G vector

        ! Reduced
        !   ik --> ik'
        integer(4), allocatable :: ik2ikqmt(:)
        !   ik --> ig
        integer(4), allocatable :: ik2ig(:)
        !   ik' --> ik
        integer(4), allocatable :: ikqmt2ik(:)

        ! Non-reduced
        !   ik --> ik'
        integer(4), allocatable :: ik2ikqmt_nr(:)
        !   ik --> ig
        integer(4), allocatable :: ik2ig_nr(:)
        !   ik' --> ik
        integer(4), allocatable :: ikqmt2ik_nr(:)

    end type kkqmt_set

!-------------------------------------------------------------------------------    
    type km_set

      ! The k-grid build by inverting the input k-grid
      ! and mapping the result back to the unit cell.
      ! It differs form the input grid, if vkloff is non-zero
      type(k_set) :: kset

      ! Index mapping between k-grid and -k-grid
      integer(4), allocatable :: ik2ikm(:)
      integer(4), allocatable :: ikm2ik(:)
      integer(4), allocatable :: ik2ikm_nr(:)
      integer(4), allocatable :: ikm2ik_nr(:)

    end type km_set

!-------------------------------------------------------------------------------    
    type q_set
        ! q-grid, i.e. the differece vectors k'-k
        type(k_set) :: qset

        ! Non-reduced maps from k'-k combinations to q grid
        ! k'-k = q + G, where k',k and q are in [0,1) unit cell
        ! iknr,iknr' --> iqnr + ig
        integer(4), allocatable :: ikikp2iq_nr(:,:)
        integer(4), allocatable :: ikikp2ig_nr(:,:)
        ! iknr,iqnr --> iknr' + ig
        integer(4), allocatable :: ikiq2ikp_nr(:,:)
        integer(4), allocatable :: ikiq2ig_nr(:,:)

        ! 1d index mapping  ikkpnr --> iqnr + ig for iknr' >= iknr
        integer(4), allocatable :: ikkp2iq_nr(:)
        integer(4), allocatable :: ikkp2ig_nr(:)

        ! ikkp with ik' >= ik ordered according to corresponding
        ! iq value
        integer(4), allocatable :: ikkp_qordered(:)

    end type q_set

    type p_set
        ! q-grid, i.e. the negative sum vectors -(k'+k) mapped back to 
        ! the unit cell
        type(k_set) :: pset

        ! Non-reduced maps from -(k'+k) combinations to p grid
        ! -(k'+k) = p + G, where k',k and p are in [0,1) unit cell
        ! iknr,iknr' --> ipnr + ig
        integer(4), allocatable :: ikikp2ip_nr(:,:)
        integer(4), allocatable :: ikikp2ig_nr(:,:)
        ! iknr,ipnr --> iknr' + ig
        integer(4), allocatable :: ikip2ikp_nr(:,:)
        integer(4), allocatable :: ikip2ig_nr(:,:)

    end type p_set

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
        maxint = nint(1.d0/eps)
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
    end subroutine rtorat

!-------------------------------------------------------------------------------
    subroutine generate_k_vectors(self,bvec,ngridk,vkloff,reduce,uselibzint)
        use modmain, only: nsymcrys, symlat, lsplsymc
        implicit none
        type(k_set), intent(OUT) :: self
        real(8), intent(IN) :: bvec(3,3)
        integer, intent(IN) :: ngridk(3)
        real(8), intent(IN) :: vkloff(3)
        logical, intent(IN) :: reduce ! apply symmetry to reduce k-set     
        logical, intent(IN), optional :: uselibzint
        ! local variables
        integer(4) :: dvk
        integer(4) :: mnd
        integer(4) :: ikloff(3)
        integer(4) :: dkloff
        integer(4) :: i1, i2, ik, nsym, isym, lspl
        integer(4), allocatable :: symmat(:,:,:)
        integer(4), allocatable :: ivk(:,:)
        integer(4), allocatable :: iwkp(:)
        integer(4) :: iv(3)

        logical :: uselz
        real(8), parameter :: epslat=1.d-6
        real(8) :: boxl(3,4)

        ! Check optional input
        if(present(uselibzint)) then
          uselz = uselibzint
        else
          uselz = .true.
        end if
        
        ! initialize k-set
        self%usedlibzint = uselz
        self%isreduced = reduce

        if (allocated(self%bvec)) deallocate(self%bvec)
        allocate(self%bvec(3,3))
        self%bvec = bvec

        if (allocated(self%vkloff)) deallocate(self%vkloff)
        allocate(self%vkloff(3))
        if(any(abs(vkloff) > 1.0d0) .or. any(vkloff < 0.0d0)) then 
          write(*,*) "Warning(generate_k_vectors): vkloff mapped back to first k-parallelepiped"
          write(*,*) "vkloff",vkloff
          call r3frac(epslat, vkloff, iv)
          write(*,*) "vkloffp",vkloff, "kshift", iv
        end if
        self%vkloff = vkloff

        if (allocated(self%ngridk)) deallocate(self%ngridk)
        allocate(self%ngridk(3))
        self%ngridk = ngridk

        ! non reduced 
        self%nkptnr = ngridk(1)*ngridk(2)*ngridk(3)
        if (allocated(self%ivknr)) deallocate(self%ivknr)
        allocate(self%ivknr(3,self%nkptnr))
        if (allocated(self%vklnr)) deallocate(self%vklnr)
        allocate(self%vklnr(3,self%nkptnr))
        if (allocated(self%vkcnr)) deallocate(self%vkcnr)
        allocate(self%vkcnr(3,self%nkptnr))
        if (allocated(self%wkptnr)) deallocate(self%wkptnr)
        allocate(self%wkptnr(self%nkptnr))

        ! reduced
        self%nkpt = self%nkptnr ! Initialize arrays for worst case symmetry
        if (allocated(self%ivk)) deallocate(self%ivk)
        allocate(self%ivk(3,self%nkpt))
        if (allocated(self%vkl)) deallocate(self%vkl)
        allocate(self%vkl(3,self%nkpt))
        if (allocated(self%vkc)) deallocate(self%vkc)
        allocate(self%vkc(3,self%nkpt))
        if (allocated(self%wkpt)) deallocate(self%wkpt)
        allocate(self%wkpt(self%nkpt))

        ! 3d index maps
        if (allocated(self%ikmap)) deallocate(self%ikmap)
        allocate(self%ikmap(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
        if (allocated(self%ikmapnr)) deallocate(self%ikmapnr)
        allocate(self%ikmapnr(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))

        ! 1d index maps
        if (allocated(self%ik2ikp)) deallocate(self%ik2ikp)
        allocate(self%ik2ikp(self%nkpt))
        if (allocated(self%ikp2ik)) deallocate(self%ikp2ik)
        allocate(self%ikp2ik(self%nkpt))

        ! tetrahedron integration method related data
        self%ntet = 6*self%nkpt
        if (allocated(self%tnodes)) deallocate(self%tnodes)
        allocate(self%tnodes(4,self%ntet))
        if (allocated(self%wtet)) deallocate(self%wtet)
        allocate(self%wtet(self%ntet))
        self%tvol = 0.d0

        if(uselz) then 

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
          ! Note (Aurich): In kgen_exciting the index labeling is different than in genppts.
          !                In genppts the 1st dimension counts up the fastest, followed by 2 and 3
          !                In kgen_exciting the opposite is the case 3,2,1.
          !                Also offsets seem not to work properly in kgen_exciting.
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

        ! No libzint 
        else

          ! Set up k grid box
          boxl(:, 1) = self%vkloff(:) / dble(self%ngridk(:))
          boxl(:, 2) = boxl(:, 1)
          boxl(:, 3) = boxl(:, 1)
          boxl(:, 4) = boxl(:, 1)
          boxl(1, 2) = boxl(1, 2) + 1.d0
          boxl(2, 3) = boxl(2, 3) + 1.d0
          boxl(3, 4) = boxl(3, 4) + 1.d0

          call genppts(.false., .false., &
          &            self%ngridk, boxl, self%nkptnr, &
          &            self%ikmapnr, self%ivknr, self%vklnr, self%vkcnr, self%wkptnr)

          call genppts(reduce, .false., &
          &            self%ngridk, boxl, self%nkpt, &
          &            self%ikmap, self%ivk, self%vkl, self%vkc, self%wkpt)

          ! 1d non-reduced index -> 1d reduced index
          do ik =1, self%nkptnr
            self%ik2ikp(ik) =&
              & self%ikmap(self%ivknr(1,ik), self%ivknr(2,ik), self%ivknr(3,ik))
          end do

          ! 1d reduced index -> 1d non-reduced index
          do ik =1, self%nkpt
            self%ikp2ik(ik) =&
              & self%ikmapnr(self%ivk(1,ik), self%ivk(2,ik), self%ivk(3,ik))
          end do

        end if

        return
    end subroutine generate_k_vectors
    
!-------------------------------------------------------------------------------
    subroutine delete_k_vectors(self)
        type(k_set), intent(INOUT) :: self
        if (allocated(self%bvec)) deallocate(self%bvec)
        if (allocated(self%vkloff)) deallocate(self%vkloff)
        if (allocated(self%ngridk)) deallocate(self%ngridk)

        if (allocated(self%ivk)) deallocate(self%ivk)
        if (allocated(self%vkl)) deallocate(self%vkl)
        if (allocated(self%vkc)) deallocate(self%vkc)
        if (allocated(self%wkpt)) deallocate(self%wkpt)

        if (allocated(self%ivknr)) deallocate(self%ivknr)
        if (allocated(self%vklnr)) deallocate(self%vklnr)
        if (allocated(self%vkcnr)) deallocate(self%vkcnr)
        if (allocated(self%wkptnr)) deallocate(self%wkptnr)

        if (allocated(self%ikmap)) deallocate(self%ikmap)
        if (allocated(self%ikmapnr)) deallocate(self%ikmapnr)
        if (allocated(self%ik2ikp)) deallocate(self%ik2ikp)
        if (allocated(self%ikp2ik)) deallocate(self%ikp2ik)

        if (allocated(self%tnodes)) deallocate(self%tnodes)
        if (allocated(self%wtet)) deallocate(self%wtet)
    end subroutine delete_k_vectors

!-------------------------------------------------------------------------------
    subroutine print_k_vectors(self,funit)
        implicit none
        type(k_set), intent(IN) :: self
        integer,     intent(IN) :: funit
        integer :: ik, i, j
        
        call boxmsg(funit,'-','k-space basis vectors')
        write(funit,*) '< b1 b2 b3 >'
        write(funit,99) transpose(self%bvec)
        99 format(3f8.4)

        call boxmsg(funit,'-','k-grid basis vectors')
        write(funit,*) '< N1 N2 N3 >'
        write(funit,98) self%ngridk
        98 format(3i4)
        write(funit,*) '< b1/N1 b2/N2 b3/N3 >'
        write(funit,99) ((self%bvec(i,j)/self%ngridk(j), j=1,3), i=1,3)

        call boxmsg(funit,'-','k-vectors')
        write(funit,*) 'Total number of k-points: < nkptnr >', self%nkptnr
        write(funit,*) 'Symmetry reduced number of k-points: < nkpt >', self%nkpt
        write(funit,*) 'Lattice offset in k-grid coordinates < vkloff >'
        write(funit,99) self%vkloff
        write(funit,*) 'Lattice offset in lattice coordinates < vkloff/ngridk >'
        write(funit,99) self%vkloff/self%ngridk
        write(funit,*)
        if(self%usedlibzint) then
          write(funit,*) 'k-vectors list:'
          write(funit,*) '< ik   vkl    vkc    weight >'
          do ik = 1, self%nkpt
            write(funit,100) ik, self%vkl(1:3,ik), self%vkc(1:3,ik), self%wkpt(ik)
          enddo
        else
          write(funit,*) 'k-vectors list:'
          write(funit,*) '< ik   ivk    vkl    vkc    weight >'
          do ik = 1, self%nkpt
            write(funit,108) ik, self%ivk(1:3,ik), self%vkl(1:3,ik), self%vkc(1:3,ik), self%wkpt(ik)
          enddo
          write(funit,*)
          write(funit,*) 'k-vectors list (nr):'
          write(funit,*) '< iknr   ivknr    vklnr    vkcnr    weightnr >'
          do ik = 1, self%nkptnr
            write(funit,108) ik, self%ivknr(1:3,ik), self%vklnr(1:3,ik), self%vkcnr(1:3,ik),&
              & self%wkptnr(ik)
          enddo
        end if
        100 format(i4,4x,3f8.4,4x,3f8.4,4x,f8.4)
        108 format(i4,4x,3i4,4x,3f8.4,4x,3f8.4,4x,f8.4)
        write(funit,*) 'Mapping from non-reduced to reduced grid: < ik2ikp >'
        do ik = 1, self%nkptnr
          write(funit,110) ik, self%ik2ikp(ik)
        end do
        110 format(i4,4x,i4)
        write(funit,*) 'Mapping from reduced to non-reduced grid: < ikp2ik >'
        do ik = 1, self%nkpt
          write(funit,110) ik, self%ikp2ik(ik)
        end do
        if(self%usedlibzint) then
          write(funit,*)
          write(funit,*) 'Tetrahedron method info: < ntet    tvol >'
          write(funit,101) self%ntet, self%tvol
          write(funit,*) 'Nodal points of tetrahedron: < itet    tnodes    wtet >'
          do i = 1, self%ntet
            write(funit,*) i, (self%tnodes(j,i),j=1,4), self%wtet(i)
          enddo 
        end if
        101 format(i6,e16.8)
        return
    end subroutine print_k_vectors

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

        ! Reciprocal lattice basis
        if (allocated(self%bvec)) deallocate(self%bvec)
        allocate(self%bvec(3,3))
        self%bvec = bvec

        ! G grid dimensions
        if (allocated(self%intgv)) deallocate(self%intgv)
        allocate(self%intgv(3,2))
        self%intgv(:,:) = intgv(:,:)

        ! Number of G points on grid 
        self%ngrtot = (self%intgv(1,2)-self%intgv(1,1)+1)* &
        &             (self%intgv(2,2)-self%intgv(2,1)+1)* &
        &             (self%intgv(3,2)-self%intgv(3,1)+1)

        ! Cutoff G vector length 
        self%gmaxvr = gmaxvr

        ! Number of G vectors shorter than gmaxvr
        self%ngvec = 0

        ! 3D index of G vectors (integer coordinates)
        if (allocated(self%ivg)) deallocate(self%ivg)
        allocate(self%ivg(3,self%ngrtot))

        ! 3D index (integer coordinates) --> 1D index map
        if (allocated(self%ivgig)) deallocate(self%ivgig)
        allocate(self%ivgig(intgv(1,1):intgv(1,2), &
        &                   intgv(2,1):intgv(2,2), &
        &                   intgv(3,1):intgv(3,2)))

        ! Cartesian coordinates of G vectors
        if (allocated(self%vgc)) deallocate(self%vgc)
        allocate(self%vgc(3,self%ngrtot))

        ! Lengths of G vectors
        if (allocated(self%gc)) deallocate(self%gc)
        allocate(self%gc(self%ngrtot))

        ! Helper arrays
        allocate(idx(self%ngrtot))
        allocate(iar(self%ngrtot))
        allocate(rar(self%ngrtot))

        ! Generate 3D index and modulus of all G vectors
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

        ! Sort by vector length
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

        ! Make integer coordinates --> 1d index
        ! and calculate G vector Cartesian coordinates
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

        ! Find the number of vectors with G < gmaxvr
        self%ngvec = 1
        do ig = self%ngrtot, 1, -1
          if (self%gc(ig) < gmaxvr) then
            self%ngvec = ig
            exit
          end if
        end do

        deallocate(idx,iar,rar)
        return
    end subroutine generate_G_vectors
    
!-------------------------------------------------------------------------------
    subroutine delete_G_vectors(self)
        type(G_set), intent(INOUT) :: self
        if (allocated(self%bvec)) deallocate(self%bvec)
        if (allocated(self%intgv)) deallocate(self%intgv)
        if (allocated(self%vgc)) deallocate(self%vgc)
        if (allocated(self%gc)) deallocate(self%gc)
        if (allocated(self%ivg)) deallocate(self%ivg)
        if (allocated(self%ivgig)) deallocate(self%ivgig)
    end subroutine delete_G_vectors

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
    end subroutine print_G_vectors

!-------------------------------------------------------------------------------
    subroutine generate_Gk_vectors(self,kset,Gset,gkmax)
        use modmain, only: nspnfv, natmtot

        implicit none

        type(Gk_set), intent(OUT) :: self
        type(k_set),  intent(IN)  :: kset
        type(G_set),  intent(IN)  :: Gset
        real(8),      intent(IN)  :: gkmax

        ! local variables
        integer :: ispn, ik, ig, igp, igmax
        logical :: fg0
        real(8) :: v(3), t1
        integer, allocatable :: igk2ig(:,:,:), ig2igk(:,:,:)

        ! Reset self
        call delete_Gk_vectors(self)

        ! Save requested gkmax
        self%gkmax = gkmax

        ! If gkmax == 0, then instead of no vectors construct
        ! those for G=0 only
        fg0 = .false.
        if(gkmax <= 0.0d0) fg0 = .true.

        !!  Reduced (potentially) k-set

        ! Map (igk,ik,ispin) --> ig
        allocate(igk2ig(gset%ngrtot,kset%nkpt,nspnfv))
        igk2ig(:,:,:) = 0
        allocate(ig2igk(gset%ngrtot,kset%nkpt,nspnfv))
        ig2igk(:,:,:) = 0
        
        ! Determine the number of G+k combinations which satisfy |G+k|<gkmax
        allocate(self%ngk(nspnfv,kset%nkpt))
        igmax = 0
        do ispn = 1, nspnfv
          do ik = 1, kset%nkpt
            igp = 0
            do ig = 1, Gset%ngrtot
              v(:) = Gset%vgc(:,ig)+kset%vkc(:,ik)
              t1 = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
              if (t1 < gkmax .or. fg0 .and. ig == 1) then
                igp = igp+1
                igk2ig(igp,ik,ispn) = ig
                ig2igk(ig,ik,ispn) = igp
                igmax = max(igmax,ig)
              end if
            end do ! ig
            self%ngk(ispn,ik) = igp
          end do ! ik
        end do ! ispn
        
        ! Maximum number of G+k vectors over all k's
        self%ngkmax = maxval(self%ngk)
        if (self%ngkmax > Gset%ngrtot) then
          write(*,*) 'ERROR(mod_kpoints::generate_Gk_vectors) ngkmax > ngrtot'
          stop
        end if
        
        ! generate G+k data set

        ! Map (ig, ispin, ik) --> igk
        allocate(self%igigk(igmax,kset%nkpt,nspnfv))
        self%igigk(:,:,:) = ig2igk(1:igmax,:,:)
        deallocate(ig2igk)

        ! Map (igk, ispin, ik) --> ig
        allocate(self%igkig(self%ngkmax,nspnfv,kset%nkpt))

        ! Lattice coordinates of G+k(ig(ik), ispin, ik)
        allocate(self%vgkl(3,self%ngkmax,nspnfv,kset%nkpt))

        ! Cartesian coordinates of G+k(ig(ik), ispin, ik)
        allocate(self%vgkc(3,self%ngkmax,nspnfv,kset%nkpt))

        ! Length of G+k vector
        allocate(self%gkc(self%ngkmax,nspnfv,kset%nkpt))

        ! Theta, Phi coordinates
        allocate(self%tpgkc(2,self%ngkmax,nspnfv,kset%nkpt))

        ! Structure factor
        allocate(self%sfacgk(self%ngkmax,natmtot,nspnfv,kset%nkpt))

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ispn, ik, igp, ig)
!$OMP DO COLLAPSE(2)
#endif    
        do ispn = 1, nspnfv
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
            end do
            ! generate structure factors for G+k-vectors
            call gensfacgp(self%ngk(ispn,ik),self%vgkc(:,:,ispn,ik), &
            &              self%ngkmax,self%sfacgk(:,:,ispn,ik))
          end do
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
        deallocate(igk2ig)

        !! Also make non-reduced quantities
        if(kset%usedlibzint == .false.) then

          ! Map (igknr,iknr,ispin) --> ig
          allocate(igk2ig(gset%ngrtot,kset%nkptnr,nspnfv))
          igk2ig(:,:,:) = 0
          allocate(ig2igk(gset%ngrtot,kset%nkptnr,nspnfv))
          ig2igk(:,:,:) = 0


          ! Determine the number of G+k combinations which satisfy |G+k|<gkmax
          allocate(self%ngknr(nspnfv,kset%nkptnr))
          igmax = 0
          do ispn = 1, nspnfv
            do ik = 1, kset%nkptnr
              igp = 0
              do ig = 1, Gset%ngrtot
                v(:) = gset%vgc(:,ig)+kset%vkcnr(:,ik)
                t1 = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
                if (t1 < gkmax .or. fg0 .and. ig == 1) then
                  igp = igp+1
                  igk2ig(igp,ik,ispn) = ig
                  ig2igk(ig,ik,ispn) = igp
                  igmax = max(igmax,ig)
                end if
              end do ! ig
              self%ngknr(ispn,ik) = igp
            end do ! ik
          end do ! ispn
          
          ! Maximum number of G+k vectors over all k's
          self%ngknrmax = maxval(self%ngknr)
          if (self%ngknrmax > gset%ngrtot) then
            write(*,*) 'ERROR(mod_kpoints::generate_Gk_vectors) ngknrmax > ngrtot'
            stop
          end if
          
          ! generate G+k data set

          ! Map (ig, ispin, iknr) --> igknr
          allocate(self%igigknr(igmax,kset%nkptnr,nspnfv))
          self%igigknr(:,:,:) = ig2igk(1:igmax,:,:)
          deallocate(ig2igk)

          ! Map (igknr, ispin, iknr) --> ig
          allocate(self%igknrig(self%ngknrmax,nspnfv,kset%nkptnr))

          ! Lattice coordinates of G+k(ig(iknr), ispin, iknr)
          allocate(self%vgknrl(3,self%ngknrmax,nspnfv,kset%nkptnr))

          ! Cartesian coordinates of G+knr(ig(iknr), ispin, iknr)
          allocate(self%vgknrc(3,self%ngknrmax,nspnfv,kset%nkptnr))

          ! Length of G+k vector
          allocate(self%gknrc(self%ngknrmax,nspnfv,kset%nkptnr))

          ! Theta, Phi coordinates
          allocate(self%tpgknrc(2,self%ngknrmax,nspnfv,kset%nkptnr))

          ! Structure factor
          allocate(self%sfacgknr(self%ngknrmax,natmtot,nspnfv,kset%nkptnr))

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ispn, ik, igp, ig)
!$OMP DO COLLAPSE(2)
#endif    
          do ispn = 1, nspnfv
            do ik = 1, kset%nkptnr
              do igp = 1, self%ngknr(ispn,ik)
                ig = igk2ig(igp,ik,ispn)
                ! index to G-vector
                self%igknrig(igp,ispn,ik) = ig
                ! G+k-vector in lattice coordinates
                self%vgknrl(:,igp,ispn,ik) = dble(gset%ivg(:,ig))+kset%vklnr(:,ik)
                ! G+k-vector in Cartesian coordinates
                self%vgknrc(:,igp,ispn,ik) = gset%vgc(:,ig)+kset%vkcnr(:,ik)
                ! G+k-vector length and (theta, phi) coordinates
                call sphcrd(self%vgknrc(:,igp,ispn,ik),self%gknrc(igp,ispn,ik),&
                &           self%tpgknrc(:,igp,ispn,ik))
              end do
              ! generate structure factors for G+k-vectors
              call gensfacgp(self%ngknr(ispn,ik),self%vgknrc(:,:,ispn,ik), &
              &              self%ngknrmax,self%sfacgknr(:,:,ispn,ik))
            end do
          end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
          deallocate(igk2ig)

        ! No libzint
        end if

        return
    end subroutine generate_Gk_vectors
    
!-------------------------------------------------------------------------------
    subroutine delete_Gk_vectors(self)
        type(Gk_set), intent(INOUT) :: self

        if (allocated(self%ngk)) deallocate(self%ngk)
        if (allocated(self%igkig)) deallocate(self%igkig)
        if (allocated(self%igigk)) deallocate(self%igigk)
        if (allocated(self%vgkl)) deallocate(self%vgkl)
        if (allocated(self%vgkc)) deallocate(self%vgkc)
        if (allocated(self%gkc)) deallocate(self%gkc)
        if (allocated(self%tpgkc)) deallocate(self%tpgkc)
        if (allocated(self%sfacgk)) deallocate(self%sfacgk)

        if (allocated(self%ngknr)) deallocate(self%ngknr)
        if (allocated(self%igknrig)) deallocate(self%igknrig)
        if (allocated(self%igigknr)) deallocate(self%igigknr)
        if (allocated(self%vgknrl)) deallocate(self%vgknrl)
        if (allocated(self%vgknrc)) deallocate(self%vgknrc)
        if (allocated(self%gknrc)) deallocate(self%gknrc)
        if (allocated(self%tpgknrc)) deallocate(self%tpgknrc)
        if (allocated(self%sfacgknr)) deallocate(self%sfacgknr)
    end subroutine delete_Gk_vectors

!-------------------------------------------------------------------------------
    subroutine print_Gk_vectors(self,ik,funit)
        implicit none
        type(Gk_set), intent(IN) :: self
        integer,      intent(IN) :: ik
        integer,      intent(IN) :: funit
        integer :: igp
        call boxmsg(funit,'-','G+k-vectors (reduced k)')
        write(funit,*) 'G+k cutoff gkmax : < gkmax >', self%gkmax
        write(funit,*) 'Maximal number of G+k-vectors over all k: < ngkmax >', self%ngkmax
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
    end subroutine print_Gk_vectors

    subroutine print_Gknr_vectors(self,iknr,funit)
        implicit none
        type(Gk_set), intent(IN) :: self
        integer,      intent(IN) :: iknr
        integer,      intent(IN) :: funit
        integer :: igp, ig
        call boxmsg(funit,'-','G+k-vectors (non-reduced k)')
        write(funit,*) 'G+k cutoff gkmax : < gkmax >', self%gkmax
        write(funit,*) 'Maximal number of G+k-vectors over all k: < ngknrmax >', self%ngknrmax
        write(funit,*) 'List of G+k-vectors for k-point:', iknr
        write(funit,*) ' * Number of G+k combinations: < ngknr >', self%ngknr(1,iknr)
        write(funit,*) ' * < igpnr    igknrig    vgknrl    vgknrc    gknrc    tpgknrc    sfacgknr >'
        do igp = 1, self%ngknr(1,iknr)
           write(funit,103) igp, self%igknrig(igp,1,iknr), &
           &                self%vgknrl(1:3,igp,1,iknr), self%vgknrc(1:3,igp,1,iknr), &
           &                self%gknrc(igp,1,iknr), self%tpgknrc(1:2,igp,1,iknr), &
           &                self%sfacgknr(igp,1,1,iknr)
        end do
        103 format(i6,2x,i6,4x,3f8.4,4x,3f8.4,4x,f8.4,4x,2f8.4,4x,2f8.4)
        call boxmsg(funit,'-','G to G+k map')
        write(funit,*) ' * < ig  igpnr >'
        do ig = 1, size(self%igigknr,1)
          write(funit, 104) ig, self%igigknr(ig, iknr, 1)
        end do
        104 format(i6,2x,i6)
        return
    end subroutine print_Gknr_vectors

!-------------------------------------------------------------------------------
    subroutine generate_kq_vectors(self,bvec,ngridk,vkloff,reduce,uselibzint)
        implicit none
        type(kq_set), intent(OUT) :: self
        real(8), intent(IN) :: bvec(3,3)
        integer, intent(IN) :: ngridk(3)
        real(8), intent(IN) :: vkloff(3)
        logical, intent(IN) :: reduce
        logical, intent(IN), optional :: uselibzint

        ! local variables        
        integer(4) :: ikloff(3)
        integer(4) :: dkloff
        integer(4) :: ik
        integer(4) :: dvk, dvq
        integer(4), allocatable :: ivk(:,:), ivq(:,:)
        logical :: uselz

        ! Check optional input
        if(present(uselibzint)) then
          uselz = uselibzint
        else
          uselz = .true.
        end if

        ! Reset self
        call delete_kq_vectors(self)
        
        ! initialize k/q-sets
        self%nkpt = ngridk(1)*ngridk(2)*ngridk(3)
        allocate(self%vkl(3,self%nkpt))
        allocate(self%vkc(3,self%nkpt))
        allocate(self%vql(3,self%nkpt))
        allocate(self%vqc(3,self%nkpt))
        
        ! tetrahedron integration method related data
        self%ntet = 6*self%nkpt
        allocate(self%tnodes(4,self%ntet))
        allocate(self%wtet(self%ntet))
        allocate(self%kqid(self%nkpt,self%nkpt))
        allocate(self%linkq(self%ntet,self%nkpt))
        self%tvol = 0.d0
        
        ! k-mesh shift
        !call factorize(3,vkloff,ikloff,dkloff) !<-- Libbzint routine
        call rtorat(vkloff,ikloff,dkloff)
        
        ! call LibBZint library
        allocate(ivk(3,self%nkpt))
        allocate(ivq(3,self%nkpt))
        
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
    end subroutine generate_kq_vectors

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
    end subroutine generate_small_group_q
    
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
    end subroutine delete_kq_vectors

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
    end subroutine print_kq_vectors

!-------------------------------------------------------------------------------
    subroutine generate_kkqmt_vectors(self,gset,bvec,ngridk,vkloff,reduce,veclqmt,uselibzint)
        type(kkqmt_set), intent(OUT) :: self
        type(g_set), intent(in) :: gset
        real(8), intent(IN) :: bvec(3,3)
        integer, intent(IN) :: ngridk(3)
        real(8), intent(IN) :: vkloff(3)
        logical, intent(IN) :: reduce ! apply symmetry to reduce k-set     
        real(8), intent(IN) :: veclqmt(3)
        logical, intent(IN), optional :: uselibzint

        real(8), parameter :: epslat=1.d-6
        real(8) :: v1(3), vkloff_kqmt(3)
        integer :: iv(3), idxnr, idxr, ik
        logical :: uselz

        if(present(uselibzint)) then 
          uselz = uselibzint
          if(uselz == .true.) then
            write(*,*) "Error(generate_kkqmt_vectors): libzint not supported. Using genppts"
            uselz = .false.
          end if
        else
          uselz = .false.
        end if

        ! Clear self
        call delete_kkqmt_vectors(self)

        ! Save Q to type
        allocate(self%vqmtl(3))
        allocate(self%vqmtc(3))
        allocate(self%vqmtl_q(3))
        allocate(self%vqmtc_q(3))
        allocate(self%vqmtl_g(3))
        allocate(self%vqmtc_g(3))

        self%vqmtl = veclqmt

        ! Get Cartesian coordinates
        call r3mv(bvec,self%vqmtl,self%vqmtc)
        ! Get q_Q \in [0,1) unit cell and G_Q \el lattice
        self%vqmtl_q=self%vqmtl
        call r3frac(epslat, self%vqmtl_q, self%vqmtl_g)
        ! Get Cartesian coordinates
        call r3mv(bvec, self%vqmtl_q, self%vqmtc_q)
        call r3mv(bvec, dble(self%vqmtl_g), self%vqmtc_g)

        ! Generate k-set
        call generate_k_vectors(self%kset, bvec, ngridk, vkloff, reduce, uselibzint=uselz)

        ! Derive k+qmt-grid offset from qmt
        ! vkloff: Offset vector in coordinates "k-grid coordinates {b_i/N_i}"
        !   --> vkloff/ngridk: Offset vector in b_i coordinates
        ! ngridk: The N_i's 
        ! Check origin of shifted k-grid
        !   Shifted k-grid origin vector is outside [0,1) unit cell
        v1 = vkloff/ngridk + self%vqmtl_q
        if(any(v1 .ge. 1.d0)) then
          ! Replace v1 with corresponding vector in unit cell
          ! and discard shifting G vector
          call r3frac(epslat, v1, iv)
          ! v1 in k-grid coordinates
          vkloff_kqmt = v1*ngridk
          ! Get corresponding vector in first k-parallelepiped
          ! (The components of vkloff should be in [0,1) )
          if(any(vkloff_kqmt .ge. 1.d0)) then
            call r3frac(epslat, vkloff_kqmt, iv)
          end if
        !   Shifted k-grid origin vector is inside [0,1) unit cell
        !   but not within first k-parallelepiped
        else if(any(v1*ngridk .ge. 1.d0)) then
          vkloff_kqmt = v1*ngridk
          call r3frac(epslat, vkloff_kqmt, iv)
        !   Shifted k-grid origin vector is inside first k-parallelepiped 
        else
          vkloff_kqmt = v1*ngridk
        end if

        ! Generate k+qmt-set
        call generate_k_vectors(self%kqmtset, bvec, ngridk, vkloff_kqmt, reduce, uselibzint=uselz)

        ! Generate map between non reduced k and  non reduced k+qmt set

        ! Build map iknr --> ikqnr
        allocate(self%ik2ikqmt_nr(self%kset%nkptnr))
        allocate(self%ik2ig_nr(self%kset%nkptnr))
        do ik = 1, self%kset%nkptnr
          ! Build k+qmt vector from k grid
          v1 = self%kset%vklnr(:, ik) + self%vqmtl_q
          ! Map back to [0,1) if needed
          call r3frac(epslat, v1, iv)
          ! Save index of G shift 
          self%ik2ig_nr(ik) = gset%ivgig(iv(1),iv(2),iv(3))
          ! Get corresponding non-reduced 3d index of k+qmt grid
          iv = nint(v1*self%kqmtset%ngridk-self%kqmtset%vkloff)
          ! Get non-reduced 1d index form 3d index
          idxnr = self%kqmtset%ikmapnr(iv(1), iv(2), iv(3))
          ! Write map 
          self%ik2ikqmt_nr(ik) = idxnr
        end do

        ! Build map ikqnr --> iknr
        allocate(self%ikqmt2ik_nr(self%kqmtset%nkptnr))
        call sortidx(self%kqmtset%nkptnr, dble(self%ik2ikqmt_nr),self%ikqmt2ik_nr)

        ! Build map ik --> ikq
        allocate(self%ik2ikqmt(self%kset%nkpt))
        allocate(self%ik2ig(self%kset%nkpt))
        do ik = 1, self%kset%nkpt
          self%ik2ikqmt(ik) = &
            & self%kqmtset%ik2ikp( self%ik2ikqmt_nr( self%kset%ikp2ik(ik) ) )
          self%ik2ig(ik) = self%ik2ig_nr( self%kset%ikp2ik(ik) ) 
        end do

        ! Build map ikq --> ik 
        allocate(self%ikqmt2ik(self%kqmtset%nkpt))
        do ik = 1, self%kqmtset%nkpt
          self%ikqmt2ik(ik) = &
            & self%kset%ik2ikp( self%ikqmt2ik_nr( self%kqmtset%ikp2ik(ik) ) )
        end do

    end subroutine generate_kkqmt_vectors

!-------------------------------------------------------------------------------
    subroutine delete_kkqmt_vectors(self)
        type(kkqmt_set), intent(INOUT) :: self
        call delete_k_vectors(self%kset)
        call delete_k_vectors(self%kqmtset)
        if(allocated(self%vqmtl)) deallocate(self%vqmtl)
        if(allocated(self%vqmtc)) deallocate(self%vqmtc)
        if(allocated(self%vqmtl_q)) deallocate(self%vqmtl_q)
        if(allocated(self%vqmtc_q)) deallocate(self%vqmtc_q)
        if(allocated(self%vqmtl_g)) deallocate(self%vqmtl_g)
        if(allocated(self%vqmtc_g)) deallocate(self%vqmtc_g)
        if(allocated(self%ik2ikqmt)) deallocate(self%ik2ikqmt)
        if(allocated(self%ik2ig)) deallocate(self%ik2ig)
        if(allocated(self%ikqmt2ik)) deallocate(self%ikqmt2ik)
        if(allocated(self%ik2ikqmt_nr)) deallocate(self%ik2ikqmt_nr)
        if(allocated(self%ik2ig_nr)) deallocate(self%ik2ig_nr)
        if(allocated(self%ikqmt2ik_nr)) deallocate(self%ikqmt2ik_nr)
    end subroutine delete_kkqmt_vectors
!-------------------------------------------------------------------------------

    subroutine print_kkqmt_vectors(self,gset,funit)
        implicit none
        type(kkqmt_set), intent(IN) :: self
        type(g_set), intent(IN) :: gset
        integer, intent(IN) :: funit
        integer :: ik, i, j

        call boxmsg(funit,'-','k set')
        call print_k_vectors(self%kset,funit)

        call boxmsg(funit,'-','Momentum transfer vector')
        write(funit,*) 'Momentum tranfer vector in lattice coordinates < vqmtl >'
        write(funit,100) (self%vqmtl(i), i=1,3)
        write(funit,*) 'Lattice vector component of momentum tranfer vector < vqmtl_g >'
        write(funit,101) (self%vqmtl_g(i), i=1,3)
        write(funit,*) 'Unit cell part of momentum tranfer vector < vqmtl_q >'
        write(funit,100) (self%vqmtl_q(i), i=1,3)
        100 format(3E11.4)
        101 format(3i11)

        call boxmsg(funit,'-','k+qmt set')
        call print_k_vectors(self%kqmtset, funit)

        call boxmsg(funit,'-','Maps non-reduced')
        write(funit,*) 'Mapping from k to k+qmt grid k index: < ik2ikqmt_nr >'
        write(funit,*) 'Mapping from k to k+qmt grid g index: < ik2ig_nr >'
        write(funit,*) 'Mapping from k to k+qmt grid g vector: < ivg >'
        write(funit,*) '< iknr    ik2ikqmt_nr    ik2ig_nr    ivg(ig)>'
        do ik = 1, self%kset%nkptnr
          write(funit,'(6i11)') ik, self%ik2ikqmt_nr(ik),&
            self%ik2ig_nr(ik), gset%ivg(:,self%ik2ig_nr(ik))
        end do
        write(funit,*) 'Mapping from k+qmt to k grid: < ikqmt2ik_nr >'
        write(funit,*) '< ikqmt    ikqmt2ik_nr >'
        do ik = 1, self%kqmtset%nkptnr
          write(funit,'(2i11)') ik, self%ikqmt2ik_nr(ik)
        end do

        call boxmsg(funit,'-','Maps reduced')
        write(funit,*) 'Mapping from k to k+qmt grid k index: < ik2ikqmt >'
        write(funit,*) 'Mapping from k to k+qmt grid g index: < ik2ig >'
        write(funit,*) '< ik    ik2ikqmt    ik2ig>'
        do ik = 1, self%kset%nkpt
          write(funit,'(3i11)') ik, self%ik2ikqmt(ik), self%ik2ig(ik)
        end do
        write(funit,*) 'Mapping from k+qmt to k grid: < ikqmt2ik >'
        write(funit,*) '< ikqmt    ikqmt2ik >'
        do ik = 1, self%kqmtset%nkpt
          write(funit,'(2i11)') ik, self%ikqmt2ik(ik)
        end do

        return

    end subroutine print_kkqmt_vectors

!-------------------------------------------------------------------------------

    subroutine generate_km_vectors(self,kset)
        type(km_set), intent(OUT) :: self
        type(k_set), intent(IN) :: kset

        real(8), parameter :: epslat=1.d-6
        real(8) :: v1(3), vkloff_km(3)
        integer :: iv(3), idxnr, idxr, ik

        ! Clear self
        call delete_km_vectors(self)

        ! Derive -k-grid offset form k-grid offset
        ! vkloff: Offset vector in coordinates "k-grid coordinates {b_i/N_i}"
        !   --> vkloff/ngridk: Offset vector in b_i coordinates
        ! ngridk: The N_i's 

        ! Map -vkloff in lattice coordinates back to [0,1]
        v1 = -kset%vkloff/kset%ngridk 
        call r3frac(epslat, v1, iv)

        ! Check origin of shifted k-grid
        !   Shifted k-grid origin vector is outside [0,1) unit cell
        if(any(v1 .ge. 1.d0)) then
          ! Replace v1 with corresponding vector in unit cell
          ! and discard shifting G vector
          call r3frac(epslat, v1, iv)
          ! v1 in k-grid coordinates
          vkloff_km = v1*kset%ngridk
          ! Get corresponding vector in first k-parallelepiped
          ! (The components of vkloff should be in [0,1) )
          if(any(vkloff_km .ge. 1.d0)) then
            call r3frac(epslat, vkloff_km, iv)
          end if
        !   Shifted k-grid origin vector is inside [0,1) unit cell
        !   but not within first k-parallelepiped
        else if(any(v1*kset%ngridk .ge. 1.d0)) then
          vkloff_km = v1*kset%ngridk
          call r3frac(epslat, vkloff_km, iv)
        !   Shifted k-grid origin vector is inside first k-parallelepiped 
        else
          vkloff_km = v1*kset%ngridk
        end if

        ! Generate -k-set
        call generate_k_vectors(self%kset, kset%bvec, kset%ngridk,&
          & vkloff_km, kset%isreduced, uselibzint=kset%usedlibzint)

        ! Generate map between non reduced k and  non reduced -k set

        ! Build map iknr --> ikmnr
        allocate(self%ik2ikm_nr(kset%nkptnr))
        do ik = 1, kset%nkptnr
           ! Build -k vector from k grid
           v1 = -kset%vklnr(:, ik)
           ! Map back to [0,1)
           call r3frac(epslat, v1, iv)
           ! Get corresponding non-reduced 3d index of -k grid
           iv = nint(v1*self%kset%ngridk-self%kset%vkloff)
           ! Get non-reduced 1d index form 3d index
           idxnr = self%kset%ikmapnr(iv(1), iv(2), iv(3))
           ! Write map 
           self%ik2ikm_nr(ik) = idxnr
        end do

        ! Build map ikmnr --> iknr
        allocate(self%ikm2ik_nr(self%kset%nkptnr))
        call sortidx(self%kset%nkptnr, dble(self%ik2ikm_nr),self%ikm2ik_nr)

        ! Build map ik --> ikm
        allocate(self%ik2ikm(kset%nkpt))
        do ik = 1, kset%nkpt
          self%ik2ikm(ik) = &
            & self%kset%ik2ikp( self%ik2ikm_nr( kset%ikp2ik(ik) ) )
        end do

        ! Build map ikq --> ik 
        allocate(self%ikm2ik(self%kset%nkpt))
        do ik = 1, self%kset%nkpt
          self%ikm2ik(ik) = &
            & kset%ik2ikp( self%ikm2ik_nr( self%kset%ikp2ik(ik) ) )
        end do

    end subroutine generate_km_vectors

!-------------------------------------------------------------------------------
    subroutine delete_km_vectors(self)
        type(km_set), intent(INOUT) :: self
        call delete_k_vectors(self%kset)
        if(allocated(self%ik2ikm)) deallocate(self%ik2ikm)
        if(allocated(self%ikm2ik)) deallocate(self%ikm2ik)
        if(allocated(self%ik2ikm_nr)) deallocate(self%ik2ikm_nr)
        if(allocated(self%ikm2ik_nr)) deallocate(self%ikm2ik_nr)
    end subroutine delete_km_vectors
!-------------------------------------------------------------------------------

    subroutine print_km_vectors(self, kset, funit)
        implicit none
        type(km_set), intent(IN) :: self
        type(k_set), intent(IN) :: kset
        integer, intent(IN) :: funit
        integer :: ik, i, j

        call boxmsg(funit,'-','k set')
        call print_k_vectors(kset,funit)

        call boxmsg(funit,'-','-k set')
        call print_k_vectors(self%kset,funit)

        call boxmsg(funit,'-','Maps non-reduced')
        write(funit,*) 'Mapping from k to -k grid: < ik2ikm_nr >'
        write(funit,*) 'Mapping from -k to k grid: < ikm2ik_nr >'
        write(funit,*) '< iknr    ik2ikm_nr    ikm2ik_nr >'
        do ik = 1, kset%nkptnr
          write(funit,101) ik, self%ik2ikm_nr(ik), self%ikm2ik_nr(ik)
        end do
        101 format(3i11)

        call boxmsg(funit,'-','Maps reduced')
        write(funit,*) 'Mapping from k to -k grid: < ik2ikm >'
        write(funit,*) '< ik    ik2ikm >'
        do ik = 1, kset%nkpt
          write(funit,'(2i11)') ik, self%ik2ikm(ik)
        end do
        write(funit,*) 'Mapping from -k to k grid: < ikm2ik >'
        write(funit,*) '< ikm    ikm2ik >'
        do ik = 1, self%kset%nkpt
          write(funit,'(2i11)') ik, self%ikm2ik(ik)
        end do

        return

    end subroutine print_km_vectors

!-------------------------------------------------------------------------------
    subroutine generate_q_vectors(self, kset, kpset, gset, reduceq)
        type(q_set), intent(OUT) :: self
        type(k_set), intent(IN) :: kset
        type(k_set), intent(IN) :: kpset
        type(g_set), intent(IN) :: gset
        logical, intent(IN) :: reduceq

        real(8) :: delta_vkloff(3), vql(3), vkpl(3)
        real(8), parameter :: epslat=1.d-6
        integer(4) :: ivg(3), iv(3), ik, ikp, iq, nkkp, ikkp

        ! Libzint not supported rigth now
        if(kpset%usedlibzint == .true. .or. kset%usedlibzint == .true.) then
          write(*,*) 'ERROR(mod_kpoints::generate_kkq_vectors):&
            & libzint not supported'
          stop
        end if
        
        ! Check if k'-grid and k-grid are compatible
        if(any(kpset%bvec /= kset%bvec)) then
          write(*,*) 'ERROR(mod_kpoints::generate_kkq_vectors):&
            & k grids have differing basis vectors'
          stop
        end if
        if(any(kpset%ngridk /= kset%ngridk)) then
          write(*,*) 'ERROR(mod_kpoints::generate_kkq_vectors):&
            & k grids have differing grid spacings'
          stop
        end if

        ! Reset self
        call delete_q_vectors(self)

        ! Get vkloff of q grid
        delta_vkloff = kpset%vkloff-kset%vkloff
        call r3frac(epslat, delta_vkloff, iv)

        ! Make q-grid
        call generate_k_vectors(self%qset, kset%bvec, kset%ngridk,&
          & delta_vkloff, reduceq, uselibzint=kset%usedlibzint)

        ! Make index map ik,ik' --> iq + ig for non-reduced points
        allocate(self%ikikp2iq_nr(kset%nkptnr, kpset%nkptnr))
        allocate(self%ikikp2ig_nr(kset%nkptnr, kpset%nkptnr))
        
        do ik = 1, kset%nkptnr
          do ikp = 1, kpset%nkptnr

            ! Get difference vector q form k' and k
            vql = kpset%vklnr(:,ikp) - kset%vklnr(:,ik)
            ! Reduce it to [0,1) and shift vector
            call r3frac(epslat, vql, ivg)

            ! Save index of G shift 
            self%ikikp2ig_nr(ik,ikp) = gset%ivgig(ivg(1),ivg(2),ivg(3))

            ! Get corresponding non-reduced 3d index of q-grid
            iv = nint(vql*self%qset%ngridk-self%qset%vkloff)

            ! Get non-reduced 1d index form 3d index
            iq = self%qset%ikmapnr(iv(1),iv(2),iv(3))

            ! Write map ik ikp --> iq
            self%ikikp2iq_nr(ik,ikp) = iq

          end do
        end do

        ! Make index map ik, iq --> ik' + ig for non-reduced points
        allocate(self%ikiq2ikp_nr(kset%nkptnr, self%qset%nkptnr))
        allocate(self%ikiq2ig_nr(kset%nkptnr, self%qset%nkptnr))

        do ik = 1, kset%nkptnr
          do iq = 1, self%qset%nkptnr

            ! Get k' form k and q
            vkpl = kset%vklnr(:,ik) + self%qset%vklnr(:,iq)
            ! Reduce it to [0,1) and shift vector
            call r3frac(epslat, vkpl, ivg)

            ! Save index of G shift 
            self%ikiq2ig_nr(ik,iq) = gset%ivgig(ivg(1),ivg(2),ivg(3))

            ! Get corresponding non-reduced 3d index of kp-grid
            iv = nint(vkpl*kpset%ngridk-kpset%vkloff)

            ! Get non-reduced 1d index form 3d index
            ikp = kpset%ikmapnr(iv(1),iv(2),iv(3))

            ! Write map ik,iq --> ikp
            self%ikiq2ikp_nr(ik,iq) = ikp

          end do
        end do

        ! Make index map for combination index ikk' --> iq, where
        ! the combination index only includes ik' >= ik
        nkkp = kset%nkptnr*(kset%nkptnr+1)/2
        allocate(self%ikkp2iq_nr(nkkp))
        allocate(self%ikkp2ig_nr(nkkp))
        do ikkp = 1, nkkp
          call kkpmap(ikkp, kset%nkptnr, ik, ikp)
          self%ikkp2iq_nr(ikkp) = self%ikikp2iq_nr(ik, ikp)
          self%ikkp2ig_nr(ikkp) = self%ikikp2ig_nr(ik, ikp)
        end do

        ! Make q-ordered ikkp list
        allocate(self%ikkp_qordered(nkkp))
        call sortidx(nkkp, dble(self%ikkp2iq_nr), self%ikkp_qordered)

    end subroutine generate_q_vectors

!-------------------------------------------------------------------------------
    subroutine delete_q_vectors(self)
        type(q_set), intent(INOUT) :: self

        call delete_k_vectors(self%qset)

        ! ik,ik' --> iq + ig
        if(allocated(self%ikikp2iq_nr)) deallocate(self%ikikp2iq_nr)
        if(allocated(self%ikikp2ig_nr)) deallocate(self%ikikp2ig_nr)

        ! ik,iq --> ik' + ig
        if(allocated(self%ikiq2ikp_nr)) deallocate(self%ikiq2ikp_nr)
        if(allocated(self%ikiq2ig_nr)) deallocate(self%ikiq2ig_nr)

        ! ikkp --> iq + ig
        if(allocated(self%ikkp2iq_nr)) deallocate(self%ikkp2iq_nr)
        if(allocated(self%ikkp2ig_nr)) deallocate(self%ikkp2ig_nr)

        ! ikkp --> ikkp(q)
        if(allocated(self%ikkp_qordered)) deallocate(self%ikkp_qordered)

    end subroutine delete_q_vectors

!-------------------------------------------------------------------------------
    subroutine print_q_vectors(self, kset, kpset, gset, funit)
        implicit none
        type(q_set), intent(in) :: self
        type(k_set), intent(in) :: kset, kpset
        type(g_set), intent(in) :: gset
        integer(4), intent(in) :: funit

        integer(4) :: ik, ikp, iq, ikkp, nkkp, i, ig, ivg(3)

        ! Sanity checks
        if(self%qset%nkptnr /= kset%nkptnr .or. self%qset%nkptnr /= kpset%nkptnr) then
          write(*,*) "Error(mod_kpointset::print_q_vectors): k,kp,q sets have&
           & differing sizes, returning"
          return
        end if
        if( any(self%qset%bvec /= kset%bvec)&
          & .or. any(self%qset%bvec /= kpset%bvec)) then
          write(*,*) "Error(mod_kpointset::print_q_vectors): k,kp,q sets have&
           & differing lattice basis vectors."
          return
        end if
        if( any(self%qset%ngridk /= kset%ngridk)&
          & .or. any(self%qset%ngridk /= kpset%ngridk)) then
          write(*,*) "Error(mod_kpointset::print_q_vectors): k,kp,q sets have&
           & differing ngridk."
          return
        end if
          
        call boxmsg(funit,'-','k point set')
        call print_k_vectors(kset, funit)

        call boxmsg(funit,'-','kp point set')
        call print_k_vectors(kpset, funit)

        call boxmsg(funit,'-','q set')
        call print_k_vectors(self%qset,funit)

        call boxmsg(funit,'-','Mapping form k, kp grid to q grid')
        write(funit,*) 'The q points are kp - k mapped back to the unit cell'
        write(funit,*) '< iknr   ikpnr    ikikp2iq_nr    ikikp2ig_nr    ivg(ig) >'
        do ik = 1, kset%nkptnr
          do ikp = 1, kpset%nkptnr
            write(funit,'(7i11)') ik, ikp, self%ikikp2iq_nr(ik,ikp),&
              & self%ikikp2ig_nr(ik,ikp),& 
              & gset%ivg(:,self%ikikp2ig_nr(ik,ikp))
          end do
        end do
        !write(funit,*)
        !write(funit,*) 'Same mapping but computed with fuction instead of map'
        !write(funit,*) '< iknr   ikpnr    iqnr    ig    ivg(ig) >'
        !do ik = 1, kset%nkptnr
        !  do ikp = 1, kpset%nkptnr
        !    call ikpik2iqivgnr(.false., ikp, ik,&
        !      & kpset%vkloff, kset%vkloff, kset%ngridk, iq, ivg)
        !    write(funit,'(7i11)') ik, ikp, iq, gset%ivgig(ivg(1),ivg(2),ivg(3)), ivg
        !  end do
        !end do

        call boxmsg(funit,'-','Mapping form k, q grid to kp grid')
        write(funit,*) 'The kp points are k + q mapped back to the unit cell'
        write(funit,*) '< iknr    iqnr    ikiq2ikp_nr   ikiq2ig_nr    ivg(ig) >'
        do ik = 1, kset%nkptnr
          do iq = 1, self%qset%nkptnr
            write(funit,'(7i11)') ik, iq, self%ikiq2ikp_nr(ik,iq),&
              & self%ikiq2ig_nr(ik,iq), gset%ivg(1:3,self%ikiq2ig_nr(ik,iq))
          end do
        end do

        !write(funit,*)
        !write(funit,*) 'Same mapping but computed with fuction instead of map'
        !write(funit,*) '< iknr    iqnr    ikqnr    ig    ivg(ig) >'
        !do ik = 1, kset%nkptnr
        !  do iq = 1, kpset%nkptnr
        !    call ikpik2iqivgnr(.true., ik, iq,&
        !      & kset%vkloff, self%qset%vkloff, kset%ngridk, ikp, ivg)
        !    write(funit,'(7i11)') ik, iq, ikp, gset%ivgig(ivg(1),ivg(2),ivg(3)), ivg
        !  end do
        !end do

        call boxmsg(funit,'-','Mapping form kp-k combination index to q grid')
        nkkp = kset%nkptnr*(kset%nkptnr+1)/2
        write(funit,*) 'The combination index labels the nkptnr*(nkptnr+1)/2&
          & combinations of ikp and ik with ikp >= ik'

        write(funit,*) '< ikkp    ik    ikp   ikkp2iq_nr    ikkp2ig_nr    ivg(ig) >'
        do ikkp = 1, nkkp
          call kkpmap(ikkp, kset%nkptnr, ik, ikp)
          write(funit,'(8i11)') ikkp, ik, ikp, self%ikkp2iq_nr(ikkp),&
              & self%ikkp2ig_nr(ikkp),& 
              & gset%ivg(:,self%ikkp2ig_nr(ikkp))
        end do
        write(funit,*)
        write(funit,*) 'ikkp index list ordered accordind to the correspondinf iq'
        write(funit,*) '< i  ikkp_qordered  iq >'
        do i = 1, nkkp
          write(funit,'(6i11)') i, self%ikkp_qordered(i),&
           &  self%ikkp2iq_nr(self%ikkp_qordered(i)),&
           &  gset%ivg(:,self%ikkp2ig_nr(self%ikkp_qordered(i)))
        end do

        return

    end subroutine print_q_vectors

!-------------------------------------------------------------------------------
    subroutine generate_p_vectors(self, kset, kpset, gset, reducep)
        type(p_set), intent(OUT) :: self
        type(k_set), intent(IN) :: kset
        type(k_set), intent(IN) :: kpset
        type(g_set), intent(IN) :: gset
        logical, intent(IN) :: reducep

        real(8) :: delta_vkloff(3), vpl(3), vkpl(3)
        real(8), parameter :: epslat=1.d-6
        integer(4) :: ivg(3), iv(3), ik, ikp, ip, nkkp, ikkp

        ! Libzint not supported right now
        if(kpset%usedlibzint == .true. .or. kset%usedlibzint == .true.) then
          write(*,*) 'ERROR(mod_kpoints::generate_kkq_vectors):&
            & libzint not supported'
          stop
        end if
        
        ! Check if k'-grid and k-grid are compatible
        if(any(kpset%bvec /= kset%bvec)) then
          write(*,*) 'ERROR(mod_kpoints::generate_kkq_vectors):&
            & k grids have differing basis vectors'
          stop
        end if
        if(any(kpset%ngridk /= kset%ngridk)) then
          write(*,*) 'ERROR(mod_kpoints::generate_kkq_vectors):&
            & k grids have differing grid spacings'
          stop
        end if

        ! Reset self
        call delete_p_vectors(self)

        ! Get vkloff of p grid
        delta_vkloff = -(kpset%vkloff+kset%vkloff)
        call r3frac(epslat, delta_vkloff, iv)

        ! Make p-grid
        call generate_k_vectors(self%pset, kset%bvec, kset%ngridk,&
          & delta_vkloff, reducep, uselibzint=kset%usedlibzint)

        ! Make index map ik,ik' --> ip + ig for non-reduced points
        allocate(self%ikikp2ip_nr(kset%nkptnr, kpset%nkptnr))
        allocate(self%ikikp2ig_nr(kset%nkptnr, kpset%nkptnr))
        
        do ik = 1, kset%nkptnr
          do ikp = 1, kpset%nkptnr

            ! Get negative sum vector p form k' and k
            vpl = -kpset%vklnr(:,ikp) - kset%vklnr(:,ik)
            ! Reduce it to [0,1) and shift vector
            call r3frac(epslat, vpl, ivg)

            ! Save index of G shift 
            self%ikikp2ig_nr(ik,ikp) = gset%ivgig(ivg(1),ivg(2),ivg(3))

            ! Get corresponding non-reduced 3d index of p-grid
            iv = nint(vpl*self%pset%ngridk-self%pset%vkloff)

            ! Get non-reduced 1d index form 3d index
            ip = self%pset%ikmapnr(iv(1),iv(2),iv(3))

            ! Write map ik,ikp --> ip
            self%ikikp2ip_nr(ik,ikp) = ip

          end do
        end do

        ! Make index map ik, ip --> ik' + ig for non-reduced points
        allocate(self%ikip2ikp_nr(kset%nkptnr, self%pset%nkptnr))
        allocate(self%ikip2ig_nr(kset%nkptnr, self%pset%nkptnr))

        do ik = 1, kset%nkptnr
          do ip = 1, self%pset%nkptnr

            ! Get k' form k and p
            vkpl = -kset%vklnr(:,ik) - self%pset%vklnr(:,ip)
            ! Reduce it to [0,1) and shift vector
            call r3frac(epslat, vkpl, ivg)

            ! Save index of G shift 
            self%ikip2ig_nr(ik,ip) = gset%ivgig(ivg(1),ivg(2),ivg(3))

            ! Get corresponding non-reduced 3d index of kp-grid
            iv = nint(vkpl*kpset%ngridk-kpset%vkloff)

            ! Get non-reduced 1d index form 3d index
            ikp = kpset%ikmapnr(iv(1),iv(2),iv(3))

            ! Write map ik,ip --> ikp
            self%ikip2ikp_nr(ik,ip) = ikp

          end do
        end do

    end subroutine generate_p_vectors

!-------------------------------------------------------------------------------
    subroutine delete_p_vectors(self)
        type(p_set), intent(INOUT) :: self

        call delete_k_vectors(self%pset)

        ! ik,ik' --> ip + ig
        if(allocated(self%ikikp2ip_nr)) deallocate(self%ikikp2ip_nr)
        if(allocated(self%ikikp2ig_nr)) deallocate(self%ikikp2ig_nr)

        ! ik,ip --> ik' + ig
        if(allocated(self%ikip2ikp_nr)) deallocate(self%ikip2ikp_nr)
        if(allocated(self%ikip2ig_nr)) deallocate(self%ikip2ig_nr)

    end subroutine delete_p_vectors

!-------------------------------------------------------------------------------
    subroutine print_p_vectors(self, kset, kpset, gset, funit)
        implicit none
        type(p_set), intent(in) :: self
        type(k_set), intent(in) :: kset, kpset
        type(g_set), intent(in) :: gset
        integer(4), intent(in) :: funit

        integer(4) :: ik, ikp, ip, i, ig, ivg(3)

        ! Sanity checks
        if(self%pset%nkptnr /= kset%nkptnr .or. self%pset%nkptnr /= kpset%nkptnr) then
          write(*,*) "Error(mod_kpointset::print_p_vectors): k,kp,p sets have&
           & differing sizes, returning"
          return
        end if
        if( any(self%pset%bvec /= kset%bvec)&
          & .or. any(self%pset%bvec /= kpset%bvec)) then
          write(*,*) "Error(mod_kpointset::print_p_vectors): k,kp,p sets have&
           & differing lattice basis vectors."
          return
        end if
        if( any(self%pset%ngridk /= kset%ngridk)&
          & .or. any(self%pset%ngridk /= kpset%ngridk)) then
          write(*,*) "Error(mod_kpointset::print_p_vectors): k,kp,p sets have&
           & differing ngridk."
          return
        end if
          
        call boxmsg(funit,'-','k point set')
        call print_k_vectors(kset, funit)

        call boxmsg(funit,'-','kp point set')
        call print_k_vectors(kpset, funit)

        call boxmsg(funit,'-','p set')
        call print_k_vectors(self%pset,funit)

        call boxmsg(funit,'-','Mapping form k, kp grid to p grid')
        write(funit,*) 'The p points are -(kp + k) mapped back to the unit cell'
        write(funit,*) '< iknr   ikpnr    ikikp2ip_nr    ikikp2ig_nr    ivg(ig) >'
        do ik = 1, kset%nkptnr
          do ikp = 1, kpset%nkptnr
            write(funit,'(7i11)') ik, ikp, self%ikikp2ip_nr(ik,ikp),&
              & self%ikikp2ig_nr(ik,ikp),& 
              & gset%ivg(:,self%ikikp2ig_nr(ik,ikp))
          end do
        end do

        call boxmsg(funit,'-','Mapping form k, p grid to kp grid')
        write(funit,*) 'The kp points are -k-p mapped back to the unit cell'
        write(funit,*) '< iknr    ipnr    ikip2ikp_nr   ikip2ig_nr    ivg(ig) >'
        do ik = 1, kset%nkptnr
          do ip = 1, self%pset%nkptnr
            write(funit,'(7i11)') ik, ip, self%ikip2ikp_nr(ik,ip),&
              & self%ikip2ig_nr(ik,ip), gset%ivg(1:3,self%ikip2ig_nr(ik,ip))
          end do
        end do

        return

    end subroutine print_p_vectors

!-------------------------------------------------------------------------------
    subroutine igkshift(gkset, gset, ik, igshift, igk2igkp)
      type(Gk_set), intent(in) :: gkset
      type(G_set), intent(in) :: gset
      integer(4), intent(in) :: ik, igshift
      integer(4), intent(out) :: igk2igkp(:)

      integer(4) :: ivgshift(3), igk, igkp, ig, igp, ivg(3), ivgp(3)
      integer(4) :: igmax, nkmax, igigkmax
      integer(4), parameter :: ispin = 1
      
      igmax = size(gset%ivg,2)
      igigkmax = size(gkset%igigk,1)
      nkmax = size(gkset%igigk,2)

      if(ik > nkmax) then
        write(*,*) "igkshift: ik too large for passed gk set"
        write(*,*) "ik=", ik, " ikmax=", nkmax
        stop
      end if
      if(igshift > igmax) then 
        write(*,*) "igkshift: igshift too large for passed g set"
        write(*,*) "igshift=", igshift, " igmax=", igmax
        stop
      end if

      write(*,*) "igshift=", igshift

      ivgshift = gset%ivg(1:3, igshift)

      write(*,*) "ivgshift=", ivgshift

      do igk = 1, gkset%ngk(ispin, ik)

        ! Get G vector from G+k index
        ig = gkset%igkig(igk, ispin, ik)

        if(ig <= igmax) then 
          ivg = gset%ivg(1:3,ig)
        else
          write(*,*) "igkshift: ig > igmax"
          write(*,*) "ig=", ig, "igmax=", igmax
          stop
        end if

        ! Add shift G vector
        ivgp = ivg + ivgshift

        if(ivgp(1) > gset%intgv(1,2) .or. ivgp(1) < gset%intgv(1,1) &
          & .or. ivgp(2) > gset%intgv(2,2) .or. ivgp(2) < gset%intgv(2,1) & 
          & .or. ivgp(3) > gset%intgv(3,2) .or. ivgp(3) < gset%intgv(3,1)) then
          write(*,*) "shifted G vector to large for passed g set"
          write(*,*) "ivg=",ivgp
          write(*,*) "intgv=", gset%intgv
          stop
        end if

        ! Get index of shifted G vector
        igp = gset%ivgig(ivgp(1),ivgp(2),ivgp(3))

        ! Get G+k index of shifted G vector
        if(igp > igigkmax) then 
          write(*,*) "Shifted G vector not in original G+k set"
          write(*,*) "igk=",igk
          write(*,*) "ig=", ig
          write(*,*) "ivg=", ivg
          write(*,*) "igp=",igp
          write(*,*) "ivgp=",ivgp
          igkp = -1
        else
          igkp = gkset%igigk(igp, ik, ispin)
        end if

        ! Write map
        igk2igkp(igk) = igkp

      end do

    end subroutine igkshift

!-------------------------------------------------------------------------------
    subroutine ikpik2iqivgnr(fplus, ikpnr, iknr, vkploff, vkloff, ngridk, iqnr, ivg)

      logical, intent(in) :: fplus
      integer(4), intent(in) :: ikpnr, iknr
      real(8), intent(in) :: vkploff(3), vkloff(3)
      integer(4), intent(in) :: ngridk(3)

      integer(4), intent(out) :: iqnr, ivg(3)

      integer(4) :: j
      integer(4) :: ivknr(3), ivkpnr(3), ivqnr(3)
      integer(4) :: ishift(3)
      real(8) :: vqloff(3)
      real(8), parameter :: epslat=1.d-6

      ! Note: vkloff and vkploff are assumed to be given in k-coordinates,
      !       positive and smaller than 1 elementwise.
      if(fplus) then
        vqloff = vkploff+vkloff
      else
        vqloff = vkploff-vkloff
      end if

      ! vkploff +/- vkloff = vqloff + ishift
      ! Get offset corresponding q grid offset in the first
      ! k-parallelepiped and a shift
      call r3frac(epslat, vqloff, ishift)

      ! Note: 3D index of k point is generated with i_j = 0,...,ngridk(j)-1
      !       and j=1 is the fastest variing dimension.
      ivknr = i3dnr(iknr, ngridk)
      ivkpnr = i3dnr(ikpnr, ngridk)

      ! Note: kp +/- k = q + G
      !       with kp,k and q in the [0,1) unit cell

      ! Get \tilde{q} in k-coordinates
      if(fplus) then
        ivqnr = ivkpnr + ivknr + ishift
      else
        ivqnr = ivkpnr - ivknr + ishift
      end if

      ! Shift \tilde{q} to [0,1)
      ivg = 0
      do j = 1,3
        if(ivqnr(j) < 0) then
          ivg(j) = -1
          ivqnr(j) = ivqnr(j) + ngridk(j)
        else if(ivqnr(j) > ngridk(j)-1) then
          ivg(j) = 1
          ivqnr(j) = ivqnr(j) - ngridk(j)
        end if
      end do

      ! Index of the q grid point.
      iqnr = i1dnr(ivqnr, ngridk)

      contains

        pure function i3dnr(i, n)
          integer(4) :: i3dnr(3)
          integer(4), intent(in) :: i
          integer(4), intent(in) :: n(3)
          integer(4) :: t
          t = i-1
          i3dnr(3) = t/(n(2)*n(1))
          t = t-i3dnr(3)*n(2)*n(1)
          i3dnr(2) = t/n(1)
          i3dnr(1) = t-i3dnr(2)*n(1)
        end function i3dnr

        pure function i1dnr(iv, n)
          integer(4) :: i1dnr
          integer(4), intent(in) :: iv(3)
          integer(4), intent(in) :: n(3)
          i1dnr = iv(1) + iv(2)*n(1) + iv(3)*n(2)*n(3) + 1
        end function i1dnr
      
    end subroutine ikpik2iqivgnr



END MODULE

