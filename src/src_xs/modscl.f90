module modscl
  use mod_constants, only: zzero
  use modmpi
  use modxs, only: unitout

  implicit none

#define BLOCKSIZE 8

  type blacsinfo
    ! Underlying MPI communicator
    type(mpiinfo) :: mpi
    ! BLACS context
    integer(4) :: context
    ! Number of processes in grid
    integer(4) :: nprocs
    ! Number of processes in each row/column
    integer(4) :: nprows, npcols
    ! Is process (0,0)
    logical :: isroot
    ! Is active, i.e. positive process coordinates
    logical :: isactive
    ! Current process grid coordinates
    integer(4) :: myprow, mypcol
    ! Default block sizes
    integer(4) :: mblck, nblck
    ! Error flag
    integer(4) :: ierr
  end type

  type(blacsinfo) :: bi2d, bi1d

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
    ! Convinience local to global index maps
    integer(4), allocatable :: r2g(:) ! rows
    integer(4), allocatable :: c2g(:) ! columns
  end type dzmat

#ifdef SCAL
  ! BLACS/ScaLAPACK routines
  integer(4), external :: numroc, indxl2g, indxg2l, indxg2p
  integer(4), external :: blacs_pnum
#endif

  contains

    !BOP
    ! !ROUTINE: setupblacs
    ! !INTERFACE:
    subroutine setupblacs(mpicom, typechars, blacscom, np)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    !   type(mpiinfo) :: mpicom    ! Info type for the MPI communicator to be 
    !                              ! used as the basis for BLACS
    !   character(*)  :: typechars ! String to indicate 2d or 1d BLACS grid
    !   integer(4), optional :: np ! Number of ranks to be used
    ! OUT:
    !   type(blacsinfo) :: blacscom  ! Info type for the resulting BLACS grid
    !
    ! !DESCRIPTION:
    !   This routine creates an BLACS process grid on top of an existing 
    !   {\tt MPI} communicator. Grid shape options are ``2d", ``1dr" and ``1dc". 
    !   The {\tt MPI} ranks are distributed in a row-major fashion.
    !
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !
    !EOP
    !BOC

      type(mpiinfo), intent(in) :: mpicom
      character(*), intent(in) :: typechars
      integer(4), intent(in), optional :: np
      type(blacsinfo), intent(out) :: blacscom

      ! Number of used processes to be used for blacs
      integer(4) :: nprocs, nprocs2d
      ! BLACS context
      integer(4) :: ictxt
      ! Process grid shape
      integer(4) :: nprows
      integer(4) :: npcols
      ! Coordinates of current process
      integer(4) :: myprow
      integer(4) :: mypcol
      ! Default block sizes
      integer(4) :: mblck, nblck

      ! Check input
      if(present(np)) then
        nprocs = np
      else
        nprocs = mpicom%procs
      end if
      if(nprocs < 1) then 
        write(*,'("Error(setupblacs): np < 1: ",i3)') nprocs
        call terminate
      else if(nprocs > mpicom%procs) then
        write(*,'("Error(setupblacs): np > mpiprocs: ",2i4)') nprocs, mpicom%procs
        call terminate
      end if

#ifdef SCAL
      
      select case(trim(adjustl(typechars)))

        ! Make 2D process grid with row major ordering of the ranks
        case('2D','2d','grid','Grid','GRID')
          ! Make square'ish porcess grid out of nporcs processes
          npcols = int(sqrt(dble(nprocs)))
          nprows = nprocs/npcols
          nprocs2d = npcols*nprows
          ! Use MPI communicator as basis
          ictxt = mpicom%comm
          ! Make the blacs grid
          call blacs_gridinit(ictxt, 'R', nprows, npcols)
          call blacs_gridinfo(ictxt, nprows, npcols, myprow, mypcol)
          ! Set default block sizes
          mblck = BLOCKSIZE
          nblck = BLOCKSIZE
          ! Write info about BLACS grid
          if(mpicom%rank == 0) then
            write(unitout,'("Info(setupblacs): Using ", i4, " x", i4,&
              & " process grid. Ctxt(", i3, ") MPIcom(", i12, ")")')&
              & nprows, npcols, ictxt, mpicom%comm
            if(nprocs2d /= nprocs) then
              write(unitout,'("Info(setupblacs):&
                & Warning - Processes do not fit 2d grid")')
              write(unitout,'("Info (setup2dblacs):&
                & Warning - ",i2," processes not used")') nprocs-nprocs2d
            end if
          end if
          nprocs = nprocs2d

        ! Make 1D process grid (one row) 
        case('1Dr','1DR','1dr','1dR','row','Row','ROW')
          npcols = nprocs
          nprows = 1
          nprocs = nprocs
          ! Use MPI communicator as basis
          ictxt = mpicom%comm
          ! Make the blacs grid
          call blacs_gridinit(ictxt, 'C', nprows, npcols)
          call blacs_gridinfo(ictxt, nprows, npcols, myprow, mypcol)
          ! Set default block sizes
          mblck = BLOCKSIZE*BLOCKSIZE
          nblck = 1
          ! Write info about BLACS grid
          if(mpicom%rank == 0) then
            write(unitout,'("Info(setupblacs): Aux. ",i4," x",i4,&
              &" process grid. Ctxt(",i3,") MPIcom(",i12,")")')&
              & nprows, npcols, ictxt, mpicom%comm
          end if

        ! Make 1D process grid (one column) 
        case('1Dc','1DC','1dc','1dC','column','Column','COLUMN')
          npcols = 1
          nprows = nprocs
          nprocs = nprocs
          ! Use MPI communicator as basis
          ictxt = mpicom%comm
          ! Make the blacs grid
          call blacs_gridinit(ictxt, 'R', nprows, npcols)
          call blacs_gridinfo(ictxt, nprows, npcols, myprow, mypcol)
          ! Set default block sizes
          mblck = 1
          nblck = BLOCKSIZE*BLOCKSIZE
          ! Write info about BLACS grid
          if(mpicom%rank == 0) then
            write(unitout,'("Info (setupblacs): Aux. ",i4," x",i4,&
              &" process grid. Ctxt(",i12,") MPIcom(",i12,")")')&
              & nprows, npcols, ictxt, mpicom%comm
          end if

        case default
          write(*,'("Error (setupblacs): Using unknown type: ",a)') typechars
          call terminate

      end select

#else
      ictxt = -1
      nprocs = 0
      nprows = -1
      npcols = -1 
      myprow = -1 
      mypcol = -1
      mblck = 1
      nblck = 1
#endif

      ! Make output type
      blacscom%mpi = mpicom
      blacscom%context = ictxt
      blacscom%nprocs = nprocs
      blacscom%nprows = nprows
      blacscom%npcols = npcols
      blacscom%myprow = myprow
      blacscom%mypcol = mypcol
      blacscom%mblck = mblck
      blacscom%nblck = nblck
      if(myprow == 0 .and. mypcol == 0) then
        blacscom%isroot = .true.
      else
        blacscom%isroot = .false.
      end if
      if(myprow >= 0 .and. mypcol >= 0) then
        blacscom%isactive = .true.
      else
        blacscom%isactive = .false.
      end if


    end subroutine setupblacs
    !EOC

    !BOP
    ! !ROUTINE: exitblacs
    subroutine exitblacs(blacscom)
    ! !INPUT/OUTPUT PARAMTERS:
    ! IN:
    !   type(blacsinfo) :: blacscom
    !
    ! !DESCRIPTION:
    !   This routine waits untill all porcesses of
    !   the passed {\tt BLACS} grid have reached it and
    !   then terminated the grid.
    !
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !EOP
    !BOC
      type(blacsinfo), intent(in) :: blacscom
#ifdef SCAL
      call blacs_barrier(blacscom%context, 'A')
      call blacs_gridexit(blacscom%context)
#endif
    end subroutine exitblacs
    !EOC

    !BOP
    ! !ROUTINE: new_dzmat
    ! !INTERFACE:
    subroutine new_dzmat(self, nrows, ncols, binfo, rblck, cblck)
    ! !INPUTP/OUTPUT PARAMETERS:
    ! IN:
    !   integer(4) :: nrows      ! Global number of rows of self
    !   integer(4) :: ncols      ! Global number of columns of self
    !   type(blacsinfo) :: binfo ! Info type for BLACS grid the matrix 
    !                            ! is to be distrubuted over
    !   integer(4), optional :: rblck, cblck ! Blocking size of rows and columns
    ! IN/OUT:
    !   type(dzmat) :: self ! The distributed complex(8) matrix
    !
    ! !DESCRIPTION:
    !   This routine initialized a block cyclic distributed global 
    !   matrix for a given {\tt BLACS} process grid. It sets up
    !   descriptor arrays, builds index maps and allocates the local matrix.
    !   In the case of compilation without {\tt ScaLAPACK} the global matrix
    !   is identical to the local matrix.
    !
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !EOP
    !BOC

      type(dzmat), intent(inout) :: self
      integer(4), intent(in) :: nrows, ncols
      type(blacsinfo), intent(in) :: binfo
      integer(4), intent(in), optional :: rblck, cblck

      integer(4) :: con, i, j

      ! Defaults
      self%isdistributed = .false.
      self%nrows = nrows
      self%ncols = ncols
      self%nrows_loc = nrows
      self%ncols_loc = ncols
      self%context = binfo%context
      self%mblck = 1
      self%nblck = 1

#ifdef SCAL

      if(self%context == -1) then
        self%isdistributed = .false.
      else
        self%isdistributed = .true.
      end if

      if(self%isdistributed) then

        ! Set block sizes
        self%mblck = min(binfo%mblck, self%nrows)
        self%nblck = min(binfo%nblck, self%ncols)
        if(present(rblck)) self%mblck = min(rblck, self%nrows)
        if(present(cblck)) self%nblck = min(cblck, self%ncols)

        ! Get number of locas rows and columns
        self%nrows_loc = numroc(self%nrows, self%mblck, binfo%myprow, 0, binfo%nprows)
        self%ncols_loc = numroc(self%ncols, self%nblck, binfo%mypcol, 0, binfo%npcols)

        ! Write index maps
        if(allocated(self%r2g)) deallocate(self%r2g)
        if(allocated(self%c2g)) deallocate(self%c2g)
        allocate(self%r2g(self%nrows_loc))
        allocate(self%c2g(self%ncols_loc))
        do i = 1, self%nrows_loc
          self%r2g(i) = indxl2g(i, self%mblck, binfo%myprow, 0, binfo%nprows)
        end do
        do j = 1, self%ncols_loc
          self%c2g(j) = indxl2g(j, self%nblck, binfo%mypcol, 0, binfo%npcols)
        end do

        ! Make descriptor
        if(allocated(self%desc)) deallocate(self%desc)
        allocate(self%desc(9))
        call descinit(self%desc, self%nrows, self%ncols, &
          & self%mblck, self%nblck, 0, 0, self%context, max(1,self%nrows_loc), binfo%ierr)
        if(binfo%ierr /= 0) then
          write(*,'("new_dzmat@rank(",i3,1x,i3,") (Error):&
            & descinit returned non zero error code")')&
            & binfo%mpi%rank, mpiglobal%rank, binfo%ierr
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
    !EOC

    !BOP
    ! !ROUTINE: del_zmat
    ! !INTERFACE:
    subroutine del_dzmat(self)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    !   type(dzmat) :: self
    !
    ! !DESCRIPTION:
    !   This routine deallocated a distributed complex(8) matrix.
    ! 
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !EOP
    !BOC
      type(dzmat), intent(inout) :: self
      if(allocated(self%za)) deallocate(self%za)
      if(allocated(self%desc)) deallocate(self%desc)
      if(allocated(self%r2g)) deallocate(self%r2g)
      if(allocated(self%c2g)) deallocate(self%c2g)
    end subroutine del_dzmat
    !EOC

    !BOP
    ! !ROUTINE: dzmat_send2global_root
    ! !INTERFACE:
    subroutine dzmat_send2global_root(mat, dmat, binfo)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    !   type(dzmat) :: dmat      ! Distributed complex(8) matrix
    !   type(blacsinfo) :: binfo ! Info type for BLACS grid
    ! OUT:
    !   complex(8), allocatable :: mat(:,:) ! Global matrix
    !
    ! !DESCRIPTION:
    !   This routine collects a distributed complex matrix 
    !   in a global array on porcess $(0,0)$.
    !
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !EOP
    !BOC

      complex(8), allocatable, intent(out) :: mat(:,:)
      type(dzmat), intent(in) :: dmat
      type(blacsinfo), intent(in) :: binfo

      logical :: iamroot
      integer(4) :: i,j,ig,jg
      integer(4) :: context, ip, prow, pcol, npr, npc
      integer(4) :: reprow, repcol
      integer(4) :: sender_nrl, sender_ncl
      integer(4), allocatable :: irbuff(:), icbuff(:)
      complex(8), allocatable :: buff(:,:)

      context = dmat%context
      if(context /= binfo%context) then
        write(*,'("Error (dzmat_send2global_root):&
          & Matrix definded on different context:",2i4)') context, binfo%context
        call terminate
      end if

      iamroot = binfo%isroot

      if(iamroot) then

        if(allocated(mat)) deallocate(mat)
        allocate(mat(dmat%nrows,dmat%ncols))

#ifdef SCAL
        prow = binfo%myprow
        pcol = binfo%mypcol
        npr = binfo%nprows
        npc = binfo%npcols
        
        ! Root's part
        do j = 1, dmat%ncols_loc
          jg = dmat%c2g(j)
          do i = 1, dmat%nrows_loc
            ig = dmat%r2g(i)
            mat(ig,jg) = dmat%za(i,j)
          end do
        end do

        ! Get data from other processes
        do ip = 1, binfo%nprocs-1

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
    !EOC

    !BOP
    ! !ROUTINE: dzmat_send2global_all
    ! !INTERFACE:
    subroutine dzmat_send2global_all(mat, dmat, binfo)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    !   type(dzmat) :: dmat      ! Distributed complex(8) matrix
    !   type(blacsinfo) :: binfo ! Info type for BLACS grid
    ! OUT:
    !   complex(8), allocatable :: mat(:,:) ! Global matrix
    !
    ! !DESCRIPTION:
    !   This routine collects a distributed complex matrix 
    !   in a global array on all porcesses of the grid.
    !
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !EOP
    !BOC

      complex(8), allocatable, intent(out) :: mat(:,:)
      type(dzmat), intent(in) :: dmat
      type(blacsinfo), intent(in) :: binfo

      integer(4) :: i,j,ig,jg
      integer(4) :: context, ip, prow, pcol, npr, npc
      integer(4) :: reprow, repcol
      integer(4) :: sender_nrl, sender_ncl
      integer(4), allocatable :: irbuff(:), icbuff(:)
      complex(8), allocatable :: buff(:,:)

      context = dmat%context
      if(context /= binfo%context) then
        write(*,'("Error (dzmat_send2global_root):&
          & Matrix definded on different context:",2i4)') context, binfo%context
        call terminate
      end if

      if(allocated(mat)) deallocate(mat)
      allocate(mat(dmat%nrows,dmat%ncols))

#ifdef SCAL

      prow = binfo%myprow
      pcol = binfo%mypcol
      npr = binfo%nprows
      npc = binfo%npcols

      do ip = 0, binfo%nprocs-1

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
    !EOC

!    subroutine unpacktoglobal(zmat, buff, iblck, jblck, prow, pcol, npr, npc)
!      complex(8), intent(inout) :: zmat(:,:)
!      complex(8), intent(in) :: buff(:,:)
!      integer(4), intent(in) :: iblck, jblck, prow, pcol, npr, npc
!
!      integer(4) :: ig, jg, il, jl
!
!      do jl = 1, size(buff, 2)
!#ifdef SCAL
!          jg = indxl2g(jl, jblck, pcol, 0, npc)
!#endif
!        do il = 1, size(buff, 1)
!#ifdef SCAL
!          ig = indxl2g(il, iblck, prow, 0, npr)
!#else
!          ig = il
!          jg = jl
!#endif
!          zmat(ig, jg) = buff(il,jl)
!        end do
!      end do
!    end subroutine unpacktoglobal



end module modscl
