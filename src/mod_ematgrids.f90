module mod_ematgrids
  use modmpi, only: terminate
  use modinput
  use mod_kpointset
  use mod_Gvector, only: intgv 
  use mod_lattice, only: bvec

  implicit none
  private

  logical :: initialized = .false.

  ! Lattice types
  !   The k-grid, k+qmt-grid and maps between them
  type(kkqmt_set) :: kkqmtset
  !   The q-grid, i.e. the difference vectors (jk+qmt) - ik = iq
  type(q_set) :: qset
  !   G vectors
  type(g_set)  :: gset
  !   G+k, G+k+qmt, G+q vectors 
  type(gk_set) :: gkset, gkqmtset, gqset

  public :: ematgrids_write_grids
  public :: ematgrids_init, ematgrids_finalize, ematgrids_initialized
  public :: kkqmtset, qset, gset, gkset, gkqmtset, gqset

  contains

    logical function ematgrids_initialized()
      ematgrids_initialized = initialized
    end function ematgrids_initialized

    subroutine ematgrids_init(vqmtl, gkmax, gqmax_, reducek_, reduceq_)
      real(8), intent(in) :: vqmtl(3)
      real(8), intent(in) :: gkmax
      real(8), intent(in), optional :: gqmax_
      logical, intent(in), optional :: reducek_, reduceq_

      integer(4) :: ngridk(3)
      real(8) :: gmaxvr, vkloff(3), gqmax
      logical :: reducek, reduceq

      !! Setting k-space grid parameters
      ngridk = input%groundstate%ngridk
      reducek = input%groundstate%reducek
      reduceq = reducek
      if(associated(input%xs)) then
        reducek = input%xs%reducek
        reduceq = input%xs%reduceq
      end if
      vkloff = input%groundstate%vkloff
      ! Reduce k-grid?
      if(present(reducek_)) then 
        reducek = reducek_
      end if
      ! Reduce q-grid?
      if(present(reduceq_)) then 
        reduceq = reduceq_
      end if
      ! Default maximal G vecort length
      gmaxvr = input%groundstate%gmaxvr
      ! G+q cutoff
      gqmax = 0.0d0 
      if(present(gqmax_)) then
        if(gqmax_ > 0.0d0) then 
          gqmax = gqmax_ 
        end if
      else
        if(associated(input%xs)) then 
          gqmax = input%xs%gqmax
        end if
      end if

      ! Generate G set
      call generate_G_vectors(gset, bvec, intgv, gmaxvr)
      !! Setup k-space grids (not using libzint)
      ! Generate k and k'=k+qmt grids and maps between them
      call generate_kkqmt_vectors(kkqmtset, gset, bvec, ngridk, vkloff, reducek,&
        & vqmtl, uselibzint=.false.)
      ! Generate q-grid as differences of jk'-ik
      call generate_q_vectors(qset, kkqmtset%kset, kkqmtset%kqmtset, gset, reduceq)
      ! Generate G+k set
      call generate_Gk_vectors(gkset, kkqmtset%kset, gset, gkmax)
      ! G+k+qmt set 
      call generate_Gk_vectors(gkqmtset, kkqmtset%kqmtset, gset, gkmax)
      ! G+q set 
      call generate_Gk_vectors(gqset, qset%qset, gset, gqmax)

      ! Info
      write(*,*) "Info(ematgrids_init):"
      write(*,*) "  gkmax = ", gkmax
      write(*,*) "  gqmax = ", gqmax
      write(*,*) "  gmaxvr = ", gmaxvr
      write(*,*) "  reducek = ", reducek
      write(*,*) "  reduceq = ", reduceq
      write(*,*) "  ngkmax = ", gkset%ngkmax
      write(*,*) "  ngkqmtmax = ", gkqmtset%ngkmax
      write(*,*) "  ngqmax = ", gqset%ngkmax

      ! Set module flag. 
      initialized = .true.

    end subroutine ematgrids_init

    subroutine ematgrids_finalize()

      if( .not. initialized) then
        write(*,*) "Warning(mod_ematgrids::ematgrids_filanlize):&
          & Module was not initialize."
        return
      end if

      ! Free k-space grids
      call delete_kkqmt_vectors(kkqmtset)
      call delete_q_vectors(qset)
      call delete_g_vectors(gset)
      call delete_gk_vectors(gkset)
      call delete_gk_vectors(gkqmtset)
      call delete_gk_vectors(gqset)

      ! Set module flag 
      initialized = .false.
    end subroutine ematgrids_finalize

    subroutine ematgrids_write_grids()
      use m_getunit

      integer(4) :: un, i

      call system('test ! -e EMATGRIDS && mkdir EMATGRIDS')

      if( .not. initialized) then
        write(*,*) "Error(mod_ematgrids::ematgrids_write_grids):&
          & Module was not initialize."
        call terminate
      end if

      call getunit(un)
      open(unit=un, file='EMATGRIDS/emg_kkqmt.out', action='write', status='replace')
      call print_kkqmt_vectors(kkqmtset, gset, un) 
      close(un)

      call getunit(un)
      open(unit=un, file='EMATGRIDS/emg_gset.out', action='write', status='replace')
      call print_G_vectors(gset, un)
      close(un)

      call getunit(un)
      open(unit=un, file='EMATGRIDS/emg_q.out', action='write', status='replace')
      call print_q_vectors(qset, kkqmtset%kset, kkqmtset%kqmtset, gset, un) 
      close(un)

      call getunit(un)
      open(unit=un, file='EMATGRIDS/emg_gkset.out', action='write', status='replace')
      do i=1, kkqmtset%kset%nkpt
        call print_Gk_vectors(gkset, i, un)
      end do
      close(un)

      call getunit(un)
      open(unit=un, file='EMATGRIDS/emg_gkqmtset.out', action='write', status='replace')
      do i=1, kkqmtset%kqmtset%nkpt
        call print_Gk_vectors(gkqmtset, i, un)
      end do
      close(un)

      call getunit(un)
      open(unit=un, file='EMATGRIDS/emg_gqset.out', action='write', status='replace')
      do i=1, qset%qset%nkpt 
        call print_Gk_vectors(gqset, i, un)
      end do
      close(un)

      call getunit(un)
      open(unit=un, file='EMATGRIDS/emg_gqset_nr.out', action='write', status='replace')
      do i=1, qset%qset%nkptnr
        call print_Gknr_vectors(gqset, i, un)
      end do
      close(un)

      return
    end subroutine ematgrids_write_grids

end module mod_ematgrids
