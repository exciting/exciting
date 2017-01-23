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

  !   The k-grid, (k+qmt)-grid and maps between them
  type(kkqmt_set) :: k_kqmtp
  !   The k-grid, (k-qmt)-grid and maps between them
  type(kkqmt_set) :: k_kqmtm
  !   The -k-grid, -(k+qmt)-grid and maps between them
  type(kkqmt_set) :: mk_mkqmtp

  ! Grids & Maps for time inversion
  !   The (-k)-grid and maps between k and -k grids
  ! Note: -k and k grids are not the same, if an offset is used.
  type(km_set) :: mk
  !   The -(k+qmt)-grid and maps between k+qmt and -(k+qmt) grids
  type(km_set) :: mkqmtp

  ! The q-grid, i.e. the difference vectors 
  !   jk - ik = iq mapped back to the unit cell
  type(q_set) :: q_q
  !   (jk+qmt) - (ik+qmt) = iq mapped back to the unit cell
  type(q_set) :: qmtp_qmtp
  !   (jk+qmt) - ik = iq mapped back to the unit cell
  type(q_set) :: q_qmtp
  !    ik - (jk+qmt) = iq mapped back to the unit cell
  type(q_set) :: qmtp_q
  !   (jk-qmt) - ik = iq mapped back to the unit cell
  type(q_set) :: q_qmtm
  !   -(jk+qmt) - ik = iq mapped back to the unit cell
  type(q_set) :: q_mqmtp
  !   -jk - (ik+qmt)
  type(q_set) :: qmtp_mq

  !   G vectors
  type(g_set)  :: g

  !   G+k, G+(k+qmt), G+(k-qmt) vectors 
  type(gk_set) :: g_k, g_kqmtp, g_kqmtm
  !   G+(-k), G+(-k-qmt) vectors 
  type(gk_set) :: g_mk, g_mkqmtp

  ! G+q
  type(gk_set) :: g_qq, g_qmtpqmtp
  type(gk_set) :: g_qqmtp, g_qmtpq
  type(gk_set) :: g_qqmtm
  type(gk_set) :: g_qmqmtp
  type(gk_set) :: g_qmtpmq
  
  public :: xsgrids_write_grids
  public :: xsgrids_init, xsgrids_finalize, xsgrids_initialized

  public :: k_kqmtp, k_kqmtm
  public :: mk_mkqmtp
  public :: mk, mkqmtp

  public :: q_q, qmtp_qmtp
  public :: q_qmtp, qmtp_q
  public :: q_qmtm
  public :: q_mqmtp, qmtp_mq

  public :: g

  public :: g_k, g_kqmtp, g_kqmtm
  public :: g_mk, g_mkqmtp

  public :: g_qq, g_qmtpqmtp
  public :: g_qqmtp, g_qmtpq
  public :: g_qqmtm
  public :: g_qmqmtp, g_qmtpmq

  contains

    logical function xsgrids_initialized()
      xsgrids_initialized = initialized
    end function xsgrids_initialized

    subroutine xsgrids_init(vqmtl, gkmax, gqmax_, reducek_, reduceq_)
      real(8), intent(in) :: vqmtl(3)
      real(8), intent(in) :: gkmax
      real(8), intent(in), optional :: gqmax_
      logical, intent(in), optional :: reducek_, reduceq_

      integer(4) :: ngridk(3)
      real(8) :: gmaxvr, vkloff(3), gqmax
      logical :: reducek, reduceq

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
      ! Default maximal G vecort length
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
      ! Generate k and k'=k+qmt grids and maps between them
      call generate_kkqmt_vectors(k_kqmtp, g, bvec, ngridk, vkloff, reducek,&
        & vqmtl, uselibzint=.false.)
      ! Generate k and k'=k-qmt grids and maps between them
      call generate_kkqmt_vectors(k_kqmtm, g, bvec, ngridk, vkloff, reducek,&
        & -vqmtl, uselibzint=.false.)

      ! Generate k and k'=-k grids and maps between them
      call generate_km_vectors(mk, k_kqmtp%kset)
      ! Generate k=k+qmt and k'=-(k+qmt) set
      call generate_km_vectors(mkqmtp, k_kqmtp%kqmtset)

      ! Generate -k and -k-qmt grids and maps between them
      call generate_kkqmt_vectors(mk_mkqmtp, g, bvec, ngridk, mk%kset%vkloff, reducek,&
        & -vqmtl, uselibzint=.false.)

      ! Generate q-grid as differences of jk-ik
      call generate_q_vectors(q_q, k_kqmtp%kset, k_kqmtp%kset, g, reduceq)
      ! Generate q-grid as differences of (jk+qmt)-(ik+qmt)
      call generate_q_vectors(qmtp_qmtp, k_kqmtp%kqmtset, k_kqmtp%kqmtset, g, reduceq)
      ! Generate q-grid as differences of (jk+qmt)-ik
      call generate_q_vectors(q_qmtp, k_kqmtp%kset, k_kqmtp%kqmtset, g, reduceq)
      ! Generate q-grid as differences of ik - (jk+qmt)
      call generate_q_vectors(qmtp_q, k_kqmtp%kqmtset, k_kqmtp%kset, g, reduceq)
      ! Generate q-grid as differences of (jk-qmt)-ik
      call generate_q_vectors(q_qmtm, k_kqmtm%kset, k_kqmtm%kqmtset, g, reduceq)
      ! Generate q-grid as differences of -(jk+qmt)-ik
      call generate_q_vectors(q_mqmtp, k_kqmtp%kset, mkqmtp%kset, g, reduceq)
      ! Generate q-grid as differences of -jk - (ik+qmt)
      call generate_q_vectors(qmtp_mq, k_kqmtp%kqmtset, mk%kset, g, reduceq)

      ! G+k set
      call generate_Gk_vectors(g_k, k_kqmtp%kset, g, gkmax)
      ! G+(k+qmt) set 
      call generate_Gk_vectors(g_kqmtp, k_kqmtp%kqmtset, g, gkmax)
      ! G+(k-qmt) set 
      call generate_Gk_vectors(g_kqmtm, k_kqmtm%kqmtset, g, gkmax)
      ! G+(-k) set 
      call generate_Gk_vectors(g_mk, mk%kset, g, gkmax)
      ! G+(-k-qmt) set 
      call generate_Gk_vectors(g_mkqmtp, mkqmtp%kset, g, gkmax)

      ! G+q sets
      call generate_Gk_vectors(g_qq, q_q%qset, g, gqmax)
      call generate_Gk_vectors(g_qmtpqmtp, qmtp_qmtp%qset, g, gqmax)
      call generate_Gk_vectors(g_qqmtp, q_qmtp%qset, g, gqmax)
      call generate_Gk_vectors(g_qmtpq, qmtp_q%qset, g, gqmax)
      call generate_Gk_vectors(g_qqmtm, q_qmtm%qset, g, gqmax)
      call generate_Gk_vectors(g_qmqmtp, q_mqmtp%qset, g, gqmax)
      call generate_Gk_vectors(g_qmtpmq, qmtp_mq%qset, g, gqmax)

      ! Info
      if(.false.) then 
        write(*,*) "Info(xsgrids_init):"
        write(*,'(a, 3i5)') "  ngridk = ", ngridk
        write(*,'(a, 3E10.3)') "  kvkloff = ", k_kqmtp%kset%vkloff
        write(*,'(a, 3E10.3)') "  mkvkloff = ", mk%kset%vkloff
        write(*,'(a, 3E10.3)') "  vqmtl*ngridk = ", vqmtl*ngridk
        write(*,'(a, 3E10.3)') "  kqmtpvkloff = ", k_kqmtp%kqmtset%vkloff
        write(*,'(a, 3E10.3)') "  kqmtmvkloff = ", k_kqmtm%kqmtset%vkloff
        write(*,'(a, 3E10.3)') "  mkqmtpvkloff = ", mkqmtp%kset%vkloff
        write(*,'(a, E10.3)') "  gkmax = ", gkmax
        write(*,'(a, E10.3)') "  gqmax = ", gqmax
        write(*,'(a, E10.3)') "  gmaxvr = ", gmaxvr
        write(*,'(a, L)') "  reducek = ", reducek
        write(*,'(a, L)') "  reduceq = ", reduceq
        write(*,'(a, i5)') "  ngkmax = ", g_k%ngkmax
        write(*,'(a, i5)') "  ngmkmax = ", g_mk%ngkmax
        write(*,'(a, i5)') "  ngkqmtpmax = ", g_kqmtp%ngkmax
        write(*,'(a, i5)') "  ngkqmtmmax = ", g_kqmtm%ngkmax
        write(*,'(a, i5)') "  ngmkqmtpmax = ", g_mkqmtp%ngkmax
      end if

      ! Set module flag. 
      initialized = .true.

    end subroutine xsgrids_init

    subroutine xsgrids_finalize()

      if( .not. initialized) then
        write(*,*) "Warning(mod_xsgrids::xsgrids_filanlize):&
          & Module was not initialize."
        return
      end if

      ! Free k-space grids
      call delete_kkqmt_vectors(k_kqmtp)
      call delete_kkqmt_vectors(k_kqmtm)
      call delete_kkqmt_vectors(mk_mkqmtp)

      call delete_km_vectors(mk)
      call delete_km_vectors(mkqmtp)

      call delete_q_vectors(q_q)
      call delete_q_vectors(qmtp_qmtp)
      call delete_q_vectors(q_qmtp)
      call delete_q_vectors(qmtp_q)
      call delete_q_vectors(q_qmtm)
      call delete_q_vectors(q_mqmtp)
      call delete_q_vectors(qmtp_mq)

      call delete_G_vectors(g)

      call delete_gk_vectors(g_k)
      call delete_gk_vectors(g_kqmtp)
      call delete_gk_vectors(g_kqmtm)

      call delete_gk_vectors(g_mk)
      call delete_gk_vectors(g_mkqmtp)

      call delete_gk_vectors(g_qq)
      call delete_gk_vectors(g_qmtpqmtp)
      call delete_gk_vectors(g_qqmtp)
      call delete_gk_vectors(g_qmtpq)
      call delete_gk_vectors(g_qqmtm)
      call delete_gk_vectors(g_qmqmtp)
      call delete_gk_vectors(g_qmtpmq)

      ! Set module flag 
      initialized = .false.
    end subroutine xsgrids_finalize

    subroutine xsgrids_write_grids(iqmt)
      use m_getunit

      integer(4) :: un, i, iqmt
      character(256) :: fext, fiqmt, fdir, fname

      write(fiqmt,*) iqmt
      fext = trim(adjustl(fiqmt))//'.out'
      fdir = 'XSGRIDS/'
      call system('[[ ! -e XSGRIDS ]] && mkdir XSGRIDS')

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

      call getunit(un)
      fname = trim(adjustl(fdir))//'mk_mkqmtp_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_kkqmt_vectors(mk_mkqmtp, g, un) 
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'mk.out'
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_km_vectors(mk, k_kqmtp%kset, un) 
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'mkqmtp_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_km_vectors(mkqmtp, k_kqmtp%kqmtset, un) 
      close(un)

      ! g grid

      call getunit(un)
      fname = trim(adjustl(fdir))//'g.out'
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_G_vectors(g, un)
      close(un)

      ! q grids 

      call getunit(un)
      fname = trim(adjustl(fdir))//'q_q.out'
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_q_vectors(q_q, k_kqmtp%kset, k_kqmtp%kset, g, un) 
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'qmtp_qmtp_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_q_vectors(qmtp_qmtp, k_kqmtp%kqmtset, k_kqmtp%kqmtset, g, un) 
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'q_qmtp_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_q_vectors(q_qmtp, k_kqmtp%kset, k_kqmtp%kqmtset, g, un) 
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'qmtp_q_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_q_vectors(qmtp_q, k_kqmtp%kqmtset, k_kqmtp%kset, g, un) 
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'q_qmtm_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_q_vectors(q_qmtm, k_kqmtm%kset, k_kqmtm%kqmtset, g, un) 
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'q_mqmtp_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_q_vectors(q_mqmtp, k_kqmtp%kset, mkqmtp%kset, g, un) 
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'qmtp_mq_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      call print_q_vectors(qmtp_mq, k_kqmtp%kqmtset, mk%kset, g, un) 
      close(un)

      ! G+k 

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

      call getunit(un)
      fname = trim(adjustl(fdir))//'g_mk.out'
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, mk%kset%nkpt
        call print_Gk_vectors(g_mk, i, un)
      end do
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'g_mkqmtp_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, mkqmtp%kset%nkpt
        call print_Gk_vectors(g_mkqmtp, i, un)
      end do
      close(un)

      ! G+q 

      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qq.out'
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, q_q%qset%nkpt 
        call print_Gk_vectors(g_qq, i, un)
      end do
      close(un)
      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qq_nr.out'
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, q_q%qset%nkptnr 
        call print_Gknr_vectors(g_qq, i, un)
      end do
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qmtpqmtp_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, qmtp_qmtp%qset%nkpt 
        call print_Gk_vectors(g_qmtpqmtp, i, un)
      end do
      close(un)
      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qmtpqmtp_nr_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, qmtp_qmtp%qset%nkptnr 
        call print_Gknr_vectors(g_qmtpqmtp, i, un)
      end do
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qqmtp_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, q_qmtp%qset%nkpt
        call print_Gk_vectors(g_qqmtp, i, un)
      end do
      close(un)
      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qqmtp_nr_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, q_qmtp%qset%nkptnr
        call print_Gknr_vectors(g_qqmtp, i, un)
      end do
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qmtpq_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, qmtp_q%qset%nkpt
        call print_Gk_vectors(g_qmtpq, i, un)
      end do
      close(un)
      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qmtpq_nr_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, qmtp_q%qset%nkptnr
        call print_Gknr_vectors(g_qmtpq, i, un)
      end do
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qqmtm_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, q_qmtm%qset%nkpt 
        call print_Gk_vectors(g_qqmtm, i, un)
      end do
      close(un)
      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qqmtm_nr_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, q_qmtm%qset%nkptnr 
        call print_Gknr_vectors(g_qqmtm, i, un)
      end do
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qmqmtp_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, q_mqmtp%qset%nkpt
        call print_Gk_vectors(g_qmqmtp, i, un)
      end do
      close(un)
      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qmqmtp_nr_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, q_mqmtp%qset%nkptnr
        call print_Gknr_vectors(g_qmqmtp, i, un)
      end do
      close(un)

      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qmtpmq_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, qmtp_mq%qset%nkpt
        call print_Gk_vectors(g_qmtpmq, i, un)
      end do
      close(un)
      call getunit(un)
      fname = trim(adjustl(fdir))//'g_qmtpmq_nr_qmt'//fext
      open(unit=un, file=trim(fname), action='write', status='replace')
      do i=1, qmtp_mq%qset%nkptnr
        call print_Gknr_vectors(g_qmtpmq, i, un)
      end do
      close(un)

      return
    end subroutine xsgrids_write_grids

end module mod_xsgrids
