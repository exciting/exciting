module modsclbse
  use modmpi
  use modxs, only: unitout

  implicit none

  integer(4) :: ictxt2d, ictxt1d
  integer(4) :: myprow, mypcol, myprow1d, mypcol1d
  integer(4) :: sender_prow, sender_pcol, sender_prow1d, sender_pcol1d
  integer(4) :: nprow, npcol, nproc, nprow1d, npcol1d
  integer(4) :: sender_nrl, sender_ncl
  ! The hermitian EVP solver of ScaLAPACK needs 
  ! mblck = nblck.
  integer(4), parameter :: mblck = 8
  integer(4), parameter :: nblck = 8
  integer(4) :: sclierr

  type dzmat
    integer(4) :: nrows, ncols
    integer(4) :: nrows_loc, ncols_loc
    complex(8), pointer :: za(:,:)
    integer(4), allocatable :: desc(:)
  end type dzmat

  integer(4), external :: numroc, indxl2g

  contains

    subroutine unpacktoglobal(ham, buff, prow, pcol, npr, npc)
      complex(8), intent(inout) :: ham(:,:)
      complex(8), intent(in) :: buff(:,:)
      integer(4), intent(in) :: prow, pcol, npr, npc

      integer(4) :: ig, jg, il, jl

      do il = 1, size(buff, 1)
        do jl = 1, size(buff, 2)
          ig = indxl2g(il, mblck, prow, 0, npr)
          jg = indxl2g(jl, nblck, pcol, 0, npc)
          ham(ig, jg) = buff(il,jl)
        end do
      end do
    end subroutine unpacktoglobal

    subroutine setupblacs

      integer(4) :: iam, sysproc
      integer(4), external :: blacs_pnum

      ! Setup 2D process grid
      ! Make rectangular process grid.
      ! Warn if there are dangling processes.
      npcol = int(sqrt(dble(procs)))
      nprow = procs/npcol
      nproc = npcol*nprow
      ! 1D grid
      nprow1d = 1
      npcol1d = nproc

#ifdef SCAL
      call blacs_pinfo(iam, sysproc)
      if(iam /= rank .or. procs /= sysproc) then
        write(*,*) "Something is fishy:", iam, rank, sysproc, procs
        call terminate
      end if
      call blacs_get(-1,0,ictxt2d)
      call blacs_gridinit(ictxt2d,'R',nprow, npcol)
      call blacs_gridinfo(ictxt2d, nprow, npcol, myprow, mypcol)

      if(myprow == 0 .and. mypcol ==0) then
        write(unitout,'("setup2dblacs: Using ",i4," x",i4," process grid.")')&
          & nprow, npcol
        if(nproc /= procs) then
          write(unitout,'("setup2dblacs: Warning - Processes do not fit 2d grid")')
          write(unitout,'("setup2dblacs: Warning - ",i2," processes not used")')&
            & procs-nproc
        end if
      end if

      ! Setup corresponding 1D grid (one row)
      call blacs_pinfo(iam, sysproc)
      if(iam /= rank .or. procs /= sysproc) then
        write(*,*) "Something is fishy:", iam, rank, sysproc, procs
        call terminate
      end if
      call blacs_get(-1,0,ictxt1d)
      call blacs_gridinit(ictxt1d,'R', 1, nprow*npcol)
      call blacs_gridinfo(ictxt1d, nprow1d, npcol1d, myprow1d, mypcol1d)
#endif

    end subroutine setupblacs

    subroutine exit2dblacs
#ifdef SCAL
      call barrier
      call blacs_gridexit(ictxt2d)
#endif
    end subroutine exit2dblacs

    subroutine new_dzmat(self, nrows, ncols)

      type(dzmat), intent(inout) :: self
      integer(4), intent(in) :: nrows, ncols

      integer(4), external :: numroc

      self%nrows = nrows
      self%ncols = ncols

      allocate(self%desc(9))

#ifdef SCAL
      self%nrows_loc = numroc(nrows, mblck, myprow, 0, nprow)
      self%ncols_loc = numroc(ncols, nblck, mypcol, 0, npcol)

      call descinit(self%desc, nrows, ncols, &
        & mblck, nblck, 0, 0, ictxt2d, self%nrows_loc, sclierr)
      if(sclierr /= 0) then
        write(*,'("new_dzmat@rank",i3," (Error):&
          & descinit returned non zero error code")') rank, sclierr
      end if
#endif
      ! Allocate local array for global matrix
      allocate(self%za(self%nrows_loc, self%ncols_loc))

      ! Zero it for good measure.
      self%za = cmplx(0,0,8)

    end subroutine new_dzmat

    subroutine del_dzmat(self)
      type(dzmat), intent(inout) :: self
      deallocate(self%za)
      deallocate(self%desc)
    end subroutine del_dzmat

end module modsclbse
