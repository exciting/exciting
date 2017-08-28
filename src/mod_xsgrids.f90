! !MODULE: mod_xsgrids
! !DESCRIPTION:
!   Global variables and setup routines for
!   the shifted $\vec{k}$ grids needed for the
!   q-dependent BSE.
!
! !REVISION HISTORY:
!   Created 2016 (Aurich)
module mod_xsgrids
  use modmpi
  use modinput
  use mod_kpointset
  use mod_Gvector, only: intgv 
  use mod_lattice, only: bvec

  implicit none
  private

  logical :: initialized = .false.

  ! Lattice types

  !   The ik=k-grid, ikp=(k+qmt/2)-grid and maps between them
  type(kkqmt_set) :: k_kqmtp
  !   The ik=k-grid, ikm=(k-qmt/2)-grid and maps between them
  type(kkqmt_set) :: k_kqmtm
  ! Simple map ikm -> ikp
  integer(4), allocatable :: ikm2ikp(:)
  ! Simple map ikp -> ikm
  integer(4), allocatable :: ikp2ikm(:)

  ! The q-grid, i.e. the difference vectors 
  !   jk - ik = iq mapped back to the unit cell
  type(q_set), target :: q
  ! The p-grid, i.e. the (-k-kp)
  !   -jkp - ikm = ip mapped back to the unit cell
  type(p_set), target :: pqmt

  !   G vectors
  type(g_set)  :: g

  !   G+k, G+(k+qmt/2), G+(k-qmt/2) vectors 
  type(gk_set) :: g_k, g_kqmtp, g_kqmtm

  ! G+q
  type(gk_set) :: g_q
  ! G+p
  type(gk_set) :: g_pqmt
  
  public :: xsgrids_write_grids
  public :: xsgrids_init, xsgrids_finalize, xsgrids_initialized

  public :: k_kqmtp, k_kqmtm
  public :: q
  public :: pqmt
  public :: g
  public :: g_k, g_kqmtp, g_kqmtm
  public :: g_q, g_pqmt
  public :: ikm2ikp, ikp2ikm

  logical :: makegk, makegq

  contains

    logical function xsgrids_initialized()
      xsgrids_initialized = initialized
    end function xsgrids_initialized

    !BOP
    ! !ROUTINE: xsgrids_init
    ! !INTERFACE:
    subroutine xsgrids_init(vqmtl, gkmax, gqmax_, reducek_,&
      & reduceq_, makegk_, makegq_)
    ! !USES:
      use modxs, only: unitout
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! real(8) :: vqmtl(3) ! Momentum transfer vector in lattice coordinates
    ! real(8) :: gkmax    ! Cutoff for G+k
    ! real(8), optional :: gqmax_    ! Cutoff for G+q
    ! logical, optional :: reducek_  ! Flag for k-point symmetry reduction
    ! logical, optional :: reduceq_  ! Flag for q-point symmetry reduction
    ! logical, optional :: makegk_   ! Construct the G+k vectors or not
    ! logical, optional :: makegq_   ! Construct the G+q vectors or not
    !
    ! !DESCRIPTION:
    !   Given a momentum transfer vector $\vec{Q}_\text{mt} = \vec{G}_\text{mt} + \vec{q}_\text{mt}$,
    !   the routine sets up the associated k-grids $\{\vec{k}\}$,
    !   $\{\vec{k}_+ = \vec{k}+\vec{q}_\text{mt}/2\}$ and $\{\vec{k}_- = \vec{k}-\vec{q}_\text{mt}/2\}$.
    !   Additionally, the q-grids formed from $\{\vec{q} = \vec{k}'-\vec{k}\}$ and 
    !   $\{\vec{q} = -\vec{k}'_+ - \vec{k}_-\}$ are set up. 
    !   Various index maps are created relating the points of the different
    !   grids.
    !   If requested also the corresponding $\vec{G}+\vec{k}$ and $\vec{G}+\vec{q}$
    !   are created.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      real(8), intent(in) :: vqmtl(3)
      real(8), intent(in) :: gkmax
      real(8), intent(in), optional :: gqmax_
      logical, intent(in), optional :: reducek_, reduceq_
      logical, intent(in), optional :: makegk_, makegq_

      integer(4) :: ngridk(3), ikm, ikp
      real(8) :: gmaxvr, vkloff(3), gqmax, t0, t1
      logical :: reducek, reduceq

      character(*), parameter :: thisname = "xsgrids_init"

      ! Clear preexisting arrays and structures
      call xsgrids_finalize()

      call timesec(t0)
      
      !! Setting k-space grid parameters
      ngridk = input%xs%ngridk
      reducek = input%xs%reducek
      reduceq = input%xs%reduceq
      vkloff = input%xs%vkloff
      ! Reduce k-grid?
      if(present(reducek_)) then 
        reducek = reducek_
      end if
      ! Reduce q-grid?
      if(present(reduceq_)) then 
        reduceq = reduceq_
      end if
      ! Build G+k arrays
      if(present(makegk_)) then
        makegk = makegk_
      else
        makegk = .false.
      end if
      ! Build G+q arrays
      if(present(makegq_)) then 
        makegq = makegq_
      else
        makegq = .false.
      end if
      ! Default maximal G vector length
      gmaxvr = input%groundstate%gmaxvr
      ! G+q cutoff
      gqmax = input%xs%gqmax
      if(present(gqmax_)) then
        if(gqmax_ > 0.0d0) then 
          gqmax = gqmax_ 
        end if
      end if

      ! Generate G set
      call generate_G_vectors(g, bvec, intgv, gmaxvr)

      !! Setup k-space grids (not using libzint)
      !!  Given a k-grid, generate the symmetrically shifted k-grids
      !!  k+qmt/2 and k-qmt/2

      ! Generate k and k'=k+qmt/2 grids and maps between them
      call generate_kkqmt_vectors(k_kqmtp, g, bvec, ngridk, vkloff, reducek,&
        & vqmtl/2.0d0, uselibzint=.false.)
      ! Generate k and k'=k-qmt/2 grids and maps between them
      call generate_kkqmt_vectors(k_kqmtm, g, bvec, ngridk, vkloff, reducek,&
        & -vqmtl/2.0d0, uselibzint=.false.)

      ! Make mapping from ikm -> ikp
      allocate(ikm2ikp(k_kqmtm%kset%nkptnr))
      allocate(ikp2ikm(k_kqmtm%kset%nkptnr))
      do ikm=1, k_kqmtm%kset%nkptnr
        ikp = k_kqmtp%ik2ikqmt_nr(k_kqmtm%ikqmt2ik_nr(ikm))
        ikm2ikp(ikm) = ikp
        ikp2ikm(ikp) = ikm
      end do

      ! Generate q-grid as differences vectors jk-ik of the original k-grid 
      ! (vkloff cancels here and one is left with an unshifted k-grid)
      call generate_q_vectors(q, k_kqmtp%kset, k_kqmtp%kset, g, reduceq)

      ! Generate p-grid as sum vectors -(jkp+ikm)
      call generate_p_vectors(pqmt, k_kqmtm%kqmtset, k_kqmtp%kqmtset, g, reduceq)

      if(makegk) then 
        ! G+k set
        call generate_Gk_vectors(g_k, k_kqmtp%kset, g, gkmax)
        ! G+(k+qmt/2) set 
        call generate_Gk_vectors(g_kqmtp, k_kqmtp%kqmtset, g, gkmax)
        ! G+(k-qmt/2) set 
        call generate_Gk_vectors(g_kqmtm, k_kqmtm%kqmtset, g, gkmax)
      end if

      if(makegq) then 
        ! G+q sets
        call generate_Gk_vectors(g_q, q%qset, g, gqmax)
        ! G+p sets
        call generate_Gk_vectors(g_pqmt, pqmt%pset, g, gqmax)
      end if

      ! Info
      if(.false.) then 
        write(*,*) "Info(xsgrids_init):"
        write(*,'(a, 3i5)') "  ngridk = ", ngridk
        write(*,'(a, 3E10.3)') "  kvkloff = ", k_kqmtp%kset%vkloff
        write(*,'(a, 3E10.3)') "  vqmtl*ngridk = ", vqmtl*ngridk
        write(*,'(a, 3E10.3)') "  kqmtpvkloff = ", k_kqmtp%kqmtset%vkloff
        write(*,'(a, 3E10.3)') "  kqmtmvkloff = ", k_kqmtm%kqmtset%vkloff
        write(*,'(a, E10.3)') "  gkmax = ", gkmax
        write(*,'(a, E10.3)') "  gqmax = ", gqmax
        write(*,'(a, E10.3)') "  gmaxvr = ", gmaxvr
        write(*,'(a, L)') "  reducek = ", reducek
        write(*,'(a, L)') "  reduceq = ", reduceq
        write(*,'(a, i5)') "  ngkmax = ", g_k%ngkmax
        write(*,'(a, i5)') "  ngkqmtpmax = ", g_kqmtp%ngkmax
        write(*,'(a, i5)') "  ngkqmtmmax = ", g_kqmtm%ngkmax
      end if

      ! Set module flag. 
      initialized = .true.

      call timesec(t1)
      write(unitout,'("Info(",a,"): Time needed/s =", f12.6)') trim(thisname), t1-t0

    end subroutine xsgrids_init
    !EOC

    !BOP
    ! !ROUTINE: xsgrids_finalize
    ! !INTERFACE:
    subroutine xsgrids_finalize()
    !
    ! !DESCRIPTION:
    !   Clears module variables.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      ! Free k-space grids
      call delete_kkqmt_vectors(k_kqmtp)
      call delete_kkqmt_vectors(k_kqmtm)

      if(allocated(ikm2ikp)) deallocate(ikm2ikp)
      if(allocated(ikp2ikm)) deallocate(ikp2ikm)

      call delete_q_vectors(q)
      call delete_p_vectors(pqmt)

      call delete_G_vectors(g)

      if(makegk) then
        call delete_gk_vectors(g_k)
        call delete_gk_vectors(g_kqmtp)
        call delete_gk_vectors(g_kqmtm)
      end if

      if(makegq) then
        call delete_gk_vectors(g_q)
        call delete_gk_vectors(g_pqmt)
      end if

      ! Set module flag 
      initialized = .false.
      makegk = .false.
      makegq = .false.
    end subroutine xsgrids_finalize
    !EOC

    !BOP
    ! !ROUTINE: xsgrids_write_grids
    ! !INTERFACE:
    subroutine xsgrids_write_grids(iqmt)
    ! !USES:
      use m_getunit
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: iqmt ! Considered Q-point
    !
    ! !DESCRIPTION:
    !   Writes out detailed information about the generated
    !   k-grids. Creates the folder {\tt XSGRIDS}.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      integer(4) :: un, i, iqmt, ikm
      character(256) :: fext, fiqmt, fdir, fname

      write(fiqmt,*) iqmt
      fext = trim(adjustl(fiqmt))//'.out'
      fdir = 'XSGRIDS/'
      call system('test ! -e XSGRIDS && mkdir XSGRIDS')

      if( .not. initialized) then
        write(*,*) "Error(mod_xsgrids::xsgrids_write_grids):&
          & Module was not initialize."
        call terminate
      end if

      ! k grids

      call getunit(un)
      fname = trim(adjustl(fdir))//'k_kqmtp_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_kkqmt_vectors(k_kqmtp, g, un) 
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'k_kqmtm_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_kkqmt_vectors(k_kqmtm, g, un) 
      close(un)

      ! km kp map
      call getunit(un)
      fname = trim(adjustl(fdir))//'kqmtm_kqmtp_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      write(un,*) 'Mapping from k-qmt/2 to k+qmt/2 grid k index: < ikm2ikp_nr >'
      write(un,*) '< ikmnr    ikm2ikp_nr>'
      do ikm = 1, k_kqmtp%kset%nkptnr
        write(un,'(2i11)') ikm, ikm2ikp(ikm)
      end do
      close(un)

      ! g grid

      call getunit(un)
      fname = trim(adjustl(fdir))//'g.out'
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_G_vectors(g, un)
      close(un)

      ! q grid

      call getunit(un)
      fname = trim(adjustl(fdir))//'q.out'
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_q_vectors(q, k_kqmtp%kset, k_kqmtp%kset, g, un) 
      close(un)

      ! p grid

      call getunit(un)
      fname = trim(adjustl(fdir))//'p_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_p_vectors(pqmt, k_kqmtm%kqmtset, k_kqmtp%kqmtset, g, un) 
      close(un)

      ! G+k 
      if(makegk) then

        call getunit(un)
        fname = trim(adjustl(fdir))//'g_k.out'
        open(unit=un, file=trim(fname), action='write', status='replace')
        do i=1, k_kqmtp%kset%nkpt
          call print_Gk_vectors(g_k, i, un)
        end do
        close(un)

        call getunit(un)
        fname = trim(adjustl(fdir))//'g_kqmtp_qmt'//fext
        open(unit=un, file=trim(fname), action='write', status='replace')
        do i=1, k_kqmtp%kqmtset%nkpt
          call print_Gk_vectors(g_kqmtp, i, un)
        end do
        close(un)

        call getunit(un)
        fname = trim(adjustl(fdir))//'g_kqmtm_qmt'//fext
        open(unit=un, file=trim(fname), action='write', status='replace')
        do i=1, k_kqmtm%kqmtset%nkpt
          call print_Gk_vectors(g_kqmtm, i, un)
        end do
        close(un)

      end if

      ! G+q  and G+p
      if(makegq) then 

        ! G+q
        call getunit(un)
        fname = trim(adjustl(fdir))//'g_q.out'
        open(unit=un, file=trim(fname), action='write', status='replace')
        do i=1, q%qset%nkpt 
          call print_Gk_vectors(g_q, i, un)
        end do
        close(un)

        call getunit(un)
        fname = trim(adjustl(fdir))//'g_q_nr.out'
        open(unit=un, file=trim(fname), action='write', status='replace')
        do i=1, q%qset%nkptnr 
          call print_Gknr_vectors(g_q, i, un)
        end do
        close(un)

        ! G+p
        call getunit(un)
        fname = trim(adjustl(fdir))//'g_pqmt_qmt'//fext
        open(unit=un, file=trim(fname), action='write', status='replace')
        do i=1, pqmt%pset%nkpt
          call print_Gk_vectors(g_pqmt, i, un)
        end do
        close(un)

        call getunit(un)
        fname = trim(adjustl(fdir))//'g_pqmt_nr_qmt'//fext
        open(unit=un, file=trim(fname), action='write', status='replace')
        do i=1, pqmt%pset%nkptnr
          call print_Gknr_vectors(g_pqmt, i, un)
        end do
        close(un)

      end if

      return
    end subroutine xsgrids_write_grids
    !EOC

end module mod_xsgrids
!EOC
