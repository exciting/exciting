module modsclbse
  use mod_constants, only: zzero
  use modmpi
  use modxs, only: unitout

  implicit none

#ifdef SCAL
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

  ! 2D blocking
  ! The hermitian EVP solver of ScaLAPACK needs 
  ! mblck = nblck.
  integer(4), parameter :: mblck = 8
  integer(4), parameter :: nblck = 8
  ! 1D blocking rows
  integer(4), parameter :: mblck1d_c = 64
  integer(4), parameter :: nblck1d_c = 1
  ! 1D blocking cols
  integer(4), parameter :: mblck1d_r = 1
  integer(4), parameter :: nblck1d_r = 64

  ! Auxilliary: Senders local array dimensions
  integer(4) :: sender_nrl, sender_ncl
  ! Auxilliary: error flag
  integer(4) :: sclierr
#endif

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
  end type dzmat

#ifdef SCAL
  ! BLACS/ScaLAPACK routines
  integer(4), external :: numroc, indxl2g, indxg2l
  integer(4), external :: blacs_pnum
#endif

  contains

    subroutine setupblacs

      integer(4) :: iam, sysproc
      integer(4) :: prow, pcol

#ifdef SCAL
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
        write(unitout,'("setup2dblacs: Using ",i4," x",i4,&
          &" process grid. Ctxt(",i1,")")')&
          & nprow, npcol, ictxt2d
        if(nproc /= procs) then
          write(unitout,'("setup2dblacs: Warning - Processes do not fit 2d grid")')
          write(unitout,'("setup2dblacs: Warning - ",i2," processes not used")')&
            & procs-nproc
        end if
        write(unitout,'("setup2dblacs: Aux. ",i4," x",i4,&
          &" process grid. Ctxt(",i1,")")')&
          & nprow1d_r, npcol1d_r, ictxt1d_r
        write(unitout,'("setup2dblacs: Aux. ",i4," x",i4,&
          &" process grid. Ctxt(",i1,")")')&
          & nprow1d_c, npcol1d_c, ictxt1d_c
      end if
#endif

    end subroutine setupblacs

    subroutine exit2dblacs
#ifdef SCAL
      call barrier
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

      do il = 1, size(buff, 1)
        do jl = 1, size(buff, 2)
#ifdef SCAL
          ig = indxl2g(il, iblck, prow, 0, npr)
          jg = indxl2g(jl, jblck, pcol, 0, npc)
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

      context = dmat%context

      if(rank == 0) then

        if(allocated(mat)) deallocate(mat)
        allocate(mat(dmat%nrows,dmat%ncols))

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
        
        call unpacktoglobal(mat, dmat%za, dmat%mblck, dmat%nblck,&
          & prow, pcol, npr, npc)

        do ip = 1, nproc-1

          ! Get info about sender and sender's local matrix
          call blacs_pcoord(context, ip, reprow, repcol)
          sender_nrl = numroc(dmat%nrows, dmat%mblck, reprow, 0, npr) 
          sender_ncl = numroc(dmat%ncols, dmat%nblck, repcol, 0, npc) 

          ! Make receive buffer
          if(allocated(buff)) deallocate(buff)
          allocate(buff(sender_nrl, sender_ncl))

          ! Receive the senders local matrix
          call zgerv2d(context, sender_nrl, sender_ncl,&
            & buff, sender_nrl, reprow, repcol) 
          ! Write content to global matrix
          call unpacktoglobal(mat, buff, dmat%mblck, dmat%nblck,&
            & reprow, repcol, npr, npc)

        end do

      else

        call zgesd2d(context, dmat%nrows_loc, dmat%ncols_loc,&
          & dmat%za, dmat%nrows_loc, 0, 0)

      end if

    end subroutine dzmat_send2global_root

    subroutine dzmat_send2global_all(mat, dmat)

      complex(8), allocatable, intent(inout) :: mat(:,:)
      type(dzmat), intent(in) :: dmat

      integer(4) :: context, ip, prow, pcol, npr, npc
      complex(8), allocatable :: buff(:,:)


      if(allocated(mat)) deallocate(mat)
      allocate(mat(dmat%nrows,dmat%ncols))


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
        if(allocated(buff)) deallocate(buff)
        allocate(buff(sender_nrl, sender_ncl))
        buff = zzero
        
        if(prow == reprow .and. pcol == repcol) then
          call zgebs2d(context, 'ALL', ' ', dmat%nrows_loc, dmat%ncols_loc,&
            & dmat%za, dmat%nrows_loc)
          buff = dmat%za
        else
          call zgebr2d(context, 'ALL', ' ', sender_nrl, sender_ncl, buff,&
            & sender_nrl, reprow, repcol)
        end if

        call unpacktoglobal(mat, buff, dmat%mblck, dmat%nblck,&
          & reprow, repcol, npr, npc)

      end do

    end subroutine dzmat_send2global_all

    subroutine new_dzmat(self, nrows, ncols, context, rblck, cblck)

      type(dzmat), intent(inout) :: self
      integer(4), intent(in) :: nrows, ncols
      integer(4), intent(in), optional :: context, rblck, cblck

      integer(4) :: con

      self%nrows = nrows
      self%ncols = ncols
      self%nrows_loc = nrows
      self%ncols_loc = ncols
      self%isdistributed = .false.

#ifdef SCAL
      if(allocated(self%desc)) deallocate(self%desc)

      if(present(context)) then
        con = context
        self%context = con
        self%isdistributed = .true.
        allocate(self%desc(9))
      else
        con = -1
        self%context = con
        self%isdistributed = .false.
      end if

      if(con == ictxt2d) then
        self%mblck = mblck
        self%nblck = nblck
        if(present(rblck)) self%mblck = rblck
        if(present(cblck)) self%nblck = cblck
        self%nrows_loc = numroc(nrows, self%mblck, myprow, 0, nprow)
        self%ncols_loc = numroc(ncols, self%nblck, mypcol, 0, npcol)
      else if(con == ictxt1d_r) then
        self%mblck = mblck1d_r
        self%nblck = nblck1d_r
        if(present(rblck)) self%mblck = rblck
        if(present(cblck)) self%nblck = cblck
        self%nrows_loc = numroc(nrows, self%mblck, myprow1d_r, 0, nprow1d_r)
        self%ncols_loc = numroc(ncols, self%nblck, mypcol1d_r, 0, npcol1d_r)
      else if(con == ictxt1d_c) then
        self%mblck = mblck1d_c
        self%nblck = nblck1d_c
        if(present(rblck)) self%mblck = rblck
        if(present(cblck)) self%nblck = cblck
        self%nrows_loc = numroc(nrows, self%mblck, myprow1d_c, 0, nprow1d_c)
        self%ncols_loc = numroc(ncols, self%nblck, mypcol1d_c, 0, npcol1d_c)
      else if(con /= -1) then
        write(*,'("new_dzmat@rank",i3," (Error): invalid context ",i3)') rank, con
        call terminate
      end if

      if(self%isdistributed) then 
        call descinit(self%desc, self%nrows, self%ncols, &
          & self%mblck, self%nblck, 0, 0, self%context, self%nrows_loc, sclierr)
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
      self%za = cmplx(0,0,8)

    end subroutine new_dzmat

    subroutine del_dzmat(self)
      type(dzmat), intent(inout) :: self
      if(allocated(self%za)) deallocate(self%za)
#ifdef SCAL
      if(allocated(self%desc)) deallocate(self%desc)
#endif
    end subroutine del_dzmat

end module modsclbse
