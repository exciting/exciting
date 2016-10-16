module modscl
  use mod_constants, only: zzero
  use modmpi
  use modxs, only: unitout

  implicit none

  ! BLACS contexts
  integer(4) :: ictxt2d, ictxt1d_r, ictxt1d_c

  ! Number of used processes in all contexts
  integer(4) :: nproc
  ! Process grid shapes
  integer(4) :: nprow, nprow1d_r, nprow1d_c
  integer(4) :: npcol, npcol1d_r, npcol1d_c
  ! Coordinates of current process
  integer(4) :: myprow, myprow1d_r, myprow1d_c
  integer(4) :: mypcol, mypcol1d_r, mypcol1d_c
  ! Coordinated of remote process
  integer(4) :: reprow, reprow1d_r, reprow1d_c
  integer(4) :: repcol, repcol1d_r, repcol1d_c

#define BLOCKSIZE 8
  ! 2D blocking
  ! The hermitian EVP solver of ScaLAPACK needs 
  ! mblck = nblck.
  integer(4), parameter :: mblck = BLOCKSIZE
  integer(4), parameter :: nblck = BLOCKSIZE
  ! 1D blocking rows
  integer(4), parameter :: mblck1d_c = BLOCKSIZE
  integer(4), parameter :: nblck1d_c = 1
  ! 1D blocking cols
  integer(4), parameter :: mblck1d_r = 1
  integer(4), parameter :: nblck1d_r = BLOCKSIZE

  ! Auxilliary: Senders local array dimensions
  integer(4) :: sender_nrl, sender_ncl
  ! Auxilliary: error flag
  integer(4) :: sclierr

  ! Distributed complex matrix type
  type dzmat
    ! Distributed or not?
    logical :: isdistributed
    ! Global dimension of the matrix
    integer(4) :: nrows, ncols
    ! Dimension of local part
    integer(4) :: nrows_loc, ncols_loc
    ! Local block cyclic part
    complex(8), allocatable :: za(:,:)
    ! BLACS descriptor
    integer(4), allocatable :: desc(:)
    ! Convinience vars (included in desc)
    ! BLACS context and blocking sizes
    integer(4) :: context
    integer(4) :: mblck, nblck
    ! Convinience l2g index maps
    integer(4), allocatable :: r2g(:)
    integer(4), allocatable :: c2g(:)
  end type dzmat

#ifdef SCAL
  ! BLACS/ScaLAPACK routines
  integer(4), external :: numroc, indxl2g, indxg2l, indxg2p
  integer(4), external :: blacs_pnum
#endif

  contains

    subroutine setupblacs

      integer(4) :: iam, sysproc
      integer(4) :: prow, pcol

      ! Setup 2D process grid
      ! Make rectangular process grid.
      ! Warn if there are dangling processes.
      npcol = int(sqrt(dble(procs)))
      nprow = procs/npcol
      nproc = npcol*nprow
      ! 1D grid cols
      nprow1d_r = 1
      npcol1d_r = nproc
      ! 1D grid rows
      nprow1d_c = nproc
      npcol1d_c = 1

      ictxt2d = -1
      ictxt1d_r = -1
      ictxt1d_c = -1

#ifdef SCAL
      ! Get info about mpi environment
      call blacs_pinfo(iam, sysproc)
      if(iam /= rank .or. procs /= sysproc) then
        write(*,*) "setupblacs (ERROR): Something is fishy:", iam, rank, sysproc, procs
        call terminate
      end if

      call blacs_get(-1, 0, ictxt2d)
      ! Make 2D process grid with row major ordering of the ranks
      call blacs_gridinit(ictxt2d, 'R', nprow, npcol)
      call blacs_gridinfo(ictxt2d, nprow, npcol, myprow, mypcol)
       
      call blacs_get(-1, 0, ictxt1d_r)
      ! Make 1D process grid (one row) with row major ordering of the ranks
      call blacs_gridinit(ictxt1d_r, 'R', nprow1d_r, npcol1d_r)
      call blacs_gridinfo(ictxt1d_r, nprow1d_r, npcol1d_r, myprow1d_r, mypcol1d_r)

      call blacs_get(-1,0,ictxt1d_c)
      ! Make 1D process grid (one col) with column major ordering of the ranks
      call blacs_gridinit(ictxt1d_c, 'C', nprow1d_c, npcol1d_c)
      call blacs_gridinfo(ictxt1d_c, nprow1d_c, npcol1d_c, myprow1d_c, mypcol1d_c)

      if(rank == 0) then
        write(unitout,'("Info(setup2dblacs): Using ",i4," x",i4,&
          &" process grid. Ctxt(",i1,")")')&
          & nprow, npcol, ictxt2d
        if(nproc /= procs) then
          write(unitout,'("Info(setup2dblacs): Warning - Processes do not fit 2d grid")')
          write(unitout,'("Info(setup2dblacs): Warning - ",i2," processes not used")')&
            & procs-nproc
        end if
        write(unitout,'("Info(setup2dblacs): Aux. ",i4," x",i4,&
          &" process grid. Ctxt(",i1,")")')&
          & nprow1d_r, npcol1d_r, ictxt1d_r
        write(unitout,'("Info(setup2dblacs): Aux. ",i4," x",i4,&
          &" process grid. Ctxt(",i1,")")')&
          & nprow1d_c, npcol1d_c, ictxt1d_c
      end if
#else
      myprow=0
      mypcol=0
      myprow1d_r=0
      mypcol1d_r=0
      myprow1d_c=0
      mypcol1d_c=0
#endif

    end subroutine setupblacs

    subroutine exit2dblacs
#ifdef SCAL
      call blacs_barrier(ictxt2d, 'A')
      call blacs_gridexit(ictxt2d)
      call blacs_gridexit(ictxt1d_r)
      call blacs_gridexit(ictxt1d_c)
#endif
    end subroutine exit2dblacs

    subroutine unpacktoglobal(zmat, buff, iblck, jblck, prow, pcol, npr, npc)
      complex(8), intent(inout) :: zmat(:,:)
      complex(8), intent(in) :: buff(:,:)
      integer(4), intent(in) :: iblck, jblck, prow, pcol, npr, npc

      integer(4) :: ig, jg, il, jl

      do jl = 1, size(buff, 2)
#ifdef SCAL
          jg = indxl2g(jl, jblck, pcol, 0, npc)
#endif
        do il = 1, size(buff, 1)
#ifdef SCAL
          ig = indxl2g(il, iblck, prow, 0, npr)
#else
          ig = il
          jg = jl
#endif
          zmat(ig, jg) = buff(il,jl)
        end do
      end do
    end subroutine unpacktoglobal

    subroutine dzmat_send2global_root(mat, dmat)

      complex(8), allocatable, intent(out) :: mat(:,:)
      type(dzmat), intent(in) :: dmat

      integer(4) :: context, ip, prow, pcol, npr, npc
      complex(8), allocatable :: buff(:,:)
      integer(4) :: i,j,ig,jg
      integer(4), allocatable :: irbuff(:), icbuff(:)

      context = dmat%context

      if(rank == 0) then

        if(allocated(mat)) deallocate(mat)
        allocate(mat(dmat%nrows,dmat%ncols))

#ifdef SCAL
        if(context == ictxt2d) then
          prow = myprow
          pcol = mypcol
          npr = nprow
          npc = npcol
        else if(context == ictxt1d_r) then
          prow = myprow1d_r
          pcol = mypcol1d_r
          npr = nprow1d_r
          npc = npcol1d_r
        else if(context == ictxt1d_c) then
          prow = myprow1d_c
          pcol = mypcol1d_c
          npr = nprow1d_c
          npc = npcol1d_c
        end if
        
        ! Rank 0's part
        do j = 1, dmat%ncols_loc
          jg = dmat%c2g(j)
          do i = 1, dmat%nrows_loc
            ig = dmat%r2g(i)
            mat(ig,jg) = dmat%za(i,j)
          end do
        end do

        ! Get data from other processes
        do ip = 1, nproc-1

          ! Get info about sender and sender's local matrix
          call blacs_pcoord(context, ip, reprow, repcol)
          sender_nrl = numroc(dmat%nrows, dmat%mblck, reprow, 0, npr) 
          sender_ncl = numroc(dmat%ncols, dmat%nblck, repcol, 0, npc) 

          ! Make receive buffer
          if(allocated(buff)) deallocate(buff)
          allocate(buff(sender_nrl, sender_ncl))
          if(allocated(irbuff)) deallocate(irbuff)
          allocate(irbuff(sender_nrl))
          if(allocated(icbuff)) deallocate(icbuff)
          allocate(icbuff(sender_ncl))

          ! Receive the senders local matrix
          call zgerv2d(context, sender_nrl, sender_ncl,&
            & buff, sender_nrl, reprow, repcol) 
          ! Receive indices
          call igerv2d(context, sender_nrl, 1,&
            & irbuff, sender_nrl, reprow, repcol) 
          call igerv2d(context, sender_ncl, 1,&
            & icbuff, sender_ncl, reprow, repcol) 

          ! Write content to global matrix
          do j = 1, sender_ncl
            jg = icbuff(j)
            do i = 1, sender_nrl
              ig = irbuff(i)
              mat(ig,jg) = buff(i,j)
            end do
          end do

        end do
#else
        mat = dmat%za
#endif

      else

#ifdef SCAL
        call zgesd2d(context, dmat%nrows_loc, dmat%ncols_loc,&
          & dmat%za, dmat%nrows_loc, 0, 0)
        call igesd2d(context, dmat%nrows_loc, 1,&
          & dmat%r2g, dmat%nrows_loc, 0, 0)
        call igesd2d(context, dmat%ncols_loc, 1,&
          & dmat%c2g, dmat%ncols_loc, 0, 0)
#endif

      end if

    end subroutine dzmat_send2global_root

    subroutine dzmat_send2global_all(mat, dmat)

      complex(8), allocatable, intent(inout) :: mat(:,:)
      type(dzmat), intent(in) :: dmat

      integer(4) :: context, ip, prow, pcol, npr, npc
      complex(8), allocatable :: buff(:,:)
      integer(4) :: i,j,ig,jg
      integer(4), allocatable :: irbuff(:), icbuff(:)


      if(allocated(mat)) deallocate(mat)
      allocate(mat(dmat%nrows,dmat%ncols))

#ifdef SCAL
      context = dmat%context

      if(context == ictxt2d) then
        prow = myprow
        pcol = mypcol
        npr = nprow
        npc = npcol
      else if(context == ictxt1d_r) then
        prow = myprow1d_r
        pcol = mypcol1d_r
        npr = nprow1d_r
        npc = npcol1d_r
      else if(context == ictxt1d_c) then
        prow = myprow1d_c
        pcol = mypcol1d_c
        npr = nprow1d_c
        npc = npcol1d_c
      end if

      do ip = 0, nproc-1

        ! Get info about sender and sender's local matrix
        call blacs_pcoord(context, ip, reprow, repcol)
        sender_nrl = numroc(dmat%nrows, dmat%mblck, reprow, 0, npr) 
        sender_ncl = numroc(dmat%ncols, dmat%nblck, repcol, 0, npc) 

        ! Make receive buffers
        if(allocated(buff)) deallocate(buff)
        allocate(buff(sender_nrl, sender_ncl))
        if(allocated(irbuff)) deallocate(irbuff)
        allocate(irbuff(sender_nrl))
        if(allocated(icbuff)) deallocate(icbuff)
        allocate(icbuff(sender_ncl))
        
        ! Broadcast send
        if(prow == reprow .and. pcol == repcol) then


          ! Matrix content
          call zgebs2d(context, 'ALL', ' ', dmat%nrows_loc, dmat%ncols_loc,&
            & dmat%za, dmat%nrows_loc)

          ! Index maps
          call igebs2d(context, 'ALL', ' ', dmat%nrows_loc, 1,&
            & dmat%r2g, dmat%nrows_loc)

          call igebs2d(context, 'ALL', ' ', dmat%ncols_loc, 1,&
            & dmat%c2g, dmat%ncols_loc)

          ! Write local to global
          do j = 1, dmat%ncols_loc
            jg = dmat%c2g(j)
            do i = 1, dmat%nrows_loc
              ig = dmat%r2g(i)
              mat(ig,jg) = dmat%za(i,j)
            end do
          end do

        ! Broadcast receive
        else

          ! Matrix content
          call zgebr2d(context, 'ALL', ' ', sender_nrl, sender_ncl, buff,&
            & sender_nrl, reprow, repcol)

          ! Index maps
          call igebr2d(context, 'ALL', ' ', sender_nrl, 1, irbuff,&
            & sender_nrl, reprow, repcol)
          call igebr2d(context, 'ALL', ' ', sender_ncl, 1, icbuff,&
            & sender_ncl, reprow, repcol)

          ! Write local to global
          do j = 1, sender_ncl
            jg = icbuff(j)
            do i = 1, sender_nrl
              ig = irbuff(i)
              mat(ig,jg) = buff(i,j)
            end do
          end do

        end if

      end do
#else
      mat = dmat%za
#endif

    end subroutine dzmat_send2global_all

    subroutine new_dzmat(self, nrows, ncols, context, rblck, cblck)

      type(dzmat), intent(inout) :: self
      integer(4), intent(in) :: nrows, ncols
      integer(4), intent(in), optional :: context, rblck, cblck

      integer(4) :: con, i, j

      self%nrows = nrows
      self%ncols = ncols
      self%nrows_loc = nrows
      self%ncols_loc = ncols
      self%isdistributed = .false.
      self%context = -1
      self%mblck = 1
      self%nblck = 1

#ifdef SCAL

      if(present(context)) then
        if(context == -1) then
          con = -1
          self%context = con
          self%isdistributed = .false.
        else
          con = context
          self%context = con
          self%isdistributed = .true.
        end if
      else 
        con = -1
        self%context = con
        self%isdistributed = .false.
      end if

      if(self%isdistributed) then
        if(con == ictxt2d) then
          self%mblck = min(mblck, self%nrows)
          self%nblck = min(nblck, self%ncols)
          if(present(rblck)) self%mblck = min(rblck, self%nrows)
          if(present(cblck)) self%nblck = min(cblck, self%ncols)
          self%nrows_loc = numroc(self%nrows, self%mblck, myprow, 0, nprow)
          self%ncols_loc = numroc(self%ncols, self%nblck, mypcol, 0, npcol)
          ! Write maps
          if(allocated(self%r2g)) deallocate(self%r2g)
          if(allocated(self%c2g)) deallocate(self%c2g)
          allocate(self%r2g(self%nrows_loc))
          allocate(self%c2g(self%ncols_loc))
          do i = 1, self%nrows_loc
            self%r2g(i) = indxl2g(i, self%mblck, myprow, 0, nprow)
          end do
          do j = 1, self%ncols_loc
            self%c2g(j) = indxl2g(j, self%nblck, mypcol, 0, npcol)
          end do
        else if(con == ictxt1d_r) then
          self%mblck = min(mblck1d_r, self%nrows)
          self%nblck = min(nblck1d_r, self%ncols)
          if(present(rblck)) self%mblck = min(rblck, self%nrows)
          if(present(cblck)) self%nblck = min(cblck, self%ncols)
          self%nrows_loc = numroc(self%nrows, self%mblck, myprow1d_r, 0, nprow1d_r)
          self%ncols_loc = numroc(self%ncols, self%nblck, mypcol1d_r, 0, npcol1d_r)
          ! Write maps
          if(allocated(self%r2g)) deallocate(self%r2g)
          if(allocated(self%c2g)) deallocate(self%c2g)
          allocate(self%r2g(self%nrows_loc))
          allocate(self%c2g(self%ncols_loc))
          do i = 1, self%nrows_loc
            self%r2g(i) = indxl2g(i, self%mblck, myprow1d_r, 0, nprow1d_r)
          end do
          do j = 1, self%ncols_loc
            self%c2g(j) = indxl2g(j, self%nblck, mypcol1d_r, 0, npcol1d_r)
          end do
        else if(con == ictxt1d_c) then
          self%mblck = min(mblck1d_c, self%nrows)
          self%nblck = min(nblck1d_c, self%ncols)
          if(present(rblck)) self%mblck = min(rblck, self%nrows)
          if(present(cblck)) self%nblck = min(cblck, self%ncols)
          self%nrows_loc = numroc(nrows, self%mblck, myprow1d_c, 0, nprow1d_c)
          self%ncols_loc = numroc(ncols, self%nblck, mypcol1d_c, 0, npcol1d_c)
          ! Write maps
          if(allocated(self%r2g)) deallocate(self%r2g)
          if(allocated(self%c2g)) deallocate(self%c2g)
          allocate(self%r2g(self%nrows_loc))
          allocate(self%c2g(self%ncols_loc))
          do i = 1, self%nrows_loc
            self%r2g(i) = indxl2g(i, self%mblck, myprow1d_c, 0, nprow1d_c)
          end do
          do j = 1, self%ncols_loc
            self%c2g(j) = indxl2g(j, self%nblck, mypcol1d_c, 0, npcol1d_c)
          end do
        else
          write(*,'("new_dzmat@rank",i3," (Error): invalid context ",i3)') rank, con
          call terminate
        end if

        ! Make descriptor
        if(allocated(self%desc)) deallocate(self%desc)
        allocate(self%desc(9))
        call descinit(self%desc, self%nrows, self%ncols, &
          & self%mblck, self%nblck, 0, 0, self%context, max(1,self%nrows_loc), sclierr)
        if(sclierr /= 0) then
          write(*,'("new_dzmat@rank",i3," (Error):&
            & descinit returned non zero error code")') rank, sclierr
          call terminate
        end if

      end if
#endif

      ! Allocate local array for global matrix
      if(allocated(self%za)) deallocate(self%za)
      allocate(self%za(self%nrows_loc, self%ncols_loc))

      ! Zero it for good measure.
      if(self%nrows_loc > 0) self%za = cmplx(0,0,8)

    end subroutine new_dzmat

    subroutine del_dzmat(self)
      type(dzmat), intent(inout) :: self
      if(allocated(self%za)) deallocate(self%za)
      if(allocated(self%desc)) deallocate(self%desc)
      if(allocated(self%r2g)) deallocate(self%r2g)
      if(allocated(self%c2g)) deallocate(self%c2g)
    end subroutine del_dzmat

end module modscl
