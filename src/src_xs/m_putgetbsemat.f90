module m_putgetbsemat
  use modmpi
  use m_getunit
  use modbse

  implicit none

  contains

    !BOP
    ! !ROUTINE: b_putbsemat
    ! !INTERFACE:
    subroutine b_putbsemat(fname, tag, ikkp, iqmt, zmat, maxsize)
    ! !USES:
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! character(*) :: fname       ! Output file name 
    ! integer(4) :: tag           ! MPI communication tag
    ! integer(4) :: ikkp          ! Index of ik jk combination
    ! integer(4) :: iqmt          ! Index of momentum transfer q
    ! complex(8) :: zmat(:,:,:,:) ! Complex 4d-array
    ! integer(4) :: maxsize       ! Maximal size of zmat over all considered ikkp
    !
    ! !DESCRIPTION:
    !   The routine writes complex 4d-array to a direct access file and
    !   is intended for the use in be BSE part of the code. It is used
    !   to write the screened coulomb interaction {\tt SCCLI.OUT} and 
    !   the exchange interaction {\tt EXCLI.OUT} to file.
    !
    ! !REVISION HISTORY:
    !   Forked from {\tt putbsemat}. (Aurich)
    !EOP
    !BOC

      implicit none

      ! Arguments
      character(*), intent(in) :: fname
      integer(4), intent(in) :: tag
      integer(4), intent(in) :: ikkp, iqmt, maxsize
      complex(8), intent(in) :: zmat(:,:,:,:)

      ! Local variables
      integer(4) :: reclen, zreclen, un
      integer(4) :: ik, jk, iknr, jknr
      integer(4) :: inu, ino, jnu, jno
      complex(8) :: dummy
      integer(4) :: buffer(7)

#ifdef MPI
      integer(4) :: iproc, stat(mpi_status_size)
#endif
      ! Get dimensions of array
      inu = size(zmat, 1)
      ino = size(zmat, 2)
      jnu = size(zmat, 3)
      jno = size(zmat, 4)

      ! Get individual ik jk index form compined ikkp index
      call kkpmap(ikkp, nk_bse, ik, jk)
      ! Get absolute non-reduced iknr jknr indices
      iknr = kmap_bse_rg(ik)
      jknr = kmap_bse_rg(jk)

      ! Send/Receive buffer
      buffer = [ikkp, iknr, jknr, inu, ino, jnu, jno]

      call getunit(un)

      ! Get large enough record length (size of zmat can depend on ikkp)
      inquire(iolength=zreclen) dummy
      inquire(iolength=reclen) iqmt, wl, wu, econv, nk_bse, buffer,&
       & koulims(:,iknr), koulims(:,jknr)
      reclen = reclen + maxsize*zreclen

#ifdef MPI

      ! Send result to rank 0
      if(mpiglobal%rank .ne. 0) then
        call mpi_send(buffer, size(buffer), mpi_integer, 0,&
          & tag*2, mpiglobal%comm, mpiglobal%ierr)
        call mpi_send(zmat, size(zmat), mpi_double_complex, 0,&
          & tag, mpiglobal%comm, mpiglobal%ierr)
      end if

      if(mpiglobal%rank .eq. 0) then

        ! Recive results form current process column
        ! (see documentaion of lastproc)
        do iproc = 0, lastproc(ikkp, nkkp_bse)

          if(iproc .ne. 0) then
            ! Receive data from slaves
            call mpi_recv(buffer, size(buffer), mpi_integer, iproc, tag*2,&
              & mpiglobal%comm, stat, mpiglobal%ierr)
            call mpi_recv(zmat, size(zmat), mpi_double_complex, iproc, tag,&
              & mpiglobal%comm, stat, mpiglobal%ierr)
          end if
#endif
          ! Only the master is performing file i/o
          open(unit=un, file=trim(fname), form='unformatted', action='write',&
            & access='direct', recl=reclen)
          ! Use ikkp as record index
          write(un, rec=buffer(1)) iqmt, wl, wu, econv, nk_bse, buffer,&
            & koulims(:,buffer(2)), koulims(:,buffer(3)), zmat
          close(un)
#ifdef MPI
        end do
      end if
#endif
    end subroutine b_putbsemat
    !EOC

    !BOP
    ! !ROUTINE: b_getbsemat
    ! !INTERFACE:
    subroutine b_getbsemat(fname, iqmt, ikkp, zmat, blindread)
    ! !DESCRIPTION:
    !   This routine is used for reading the screened coulomb interaction
    !   and exchange interaction from file. It works for ou,ou combinations.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC

      implicit none

      ! Arguments
      character(*), intent(in) :: fname
      integer(4), intent(in) :: iqmt
      integer(4), intent(in) :: ikkp
      logical, intent(in), optional :: blindread
      complex(8), intent(out) :: zmat(:,:)

      ! Local variables
      real(8) :: wl_, wu_, econv_
      integer(4) :: un, reclen
      integer(4) :: iqmt_, nk_bse_, ikkp_, iknr_, jknr_, ikoulims_(4), jkoulims_(4)
      integer(4) :: inu_, ino_, jnu_, jno_
      integer(4) :: inu, ino, jnu, jno, inou, jnou, iknr, jknr, ik, jk
      integer(4) :: inu_offset, ino_offset, jnu_offset, jno_offset

      complex(8), allocatable :: zm(:, :, :, :)
      logical :: ishere, trustme

      ! Check if file exsists
      inquire(file=trim(fname), exist=ishere)
      if(.not. ishere) then 
        write(*,*) "Error (b_getbsemat): file does not exsist fname=",trim(fname)
        call terminate
      end if

      if(present(blindread)) then
        trustme = blindread
      end if

      ! Check content or trust user (what could possibly go wrong)
      if(.not. trustme) then

        ! Get array dimensions
        inou = size(zmat,1)
        jnou = size(zmat,2)

        ! Get requested quantities for comparison to stored ones
        call kkpmap(ikkp, nk_bse, ik, jk)
        iknr = kmap_bse_rg(ik)
        inu = koulims(2,iknr)-koulims(1,iknr)+1
        ino = koulims(4,iknr)-koulims(3,iknr)+1
        jknr = kmap_bse_rg(jk)
        jnu = koulims(2,jknr)-koulims(1,jknr)+1
        jno = koulims(4,jknr)-koulims(3,jknr)+1

        ! Input consitency check
        if(inou /= ino*inu .or. jnou /= inu*ino) then 
          write(*,*) "Error (b_getbsemat): Inconstistency in input matrix size"
          write(*,*) "shape(zmat) =", shape(zmat)
          write(*,*) "inu, ino, jnu, jno =", inu, ino, jnu, jno
          call terminate
        end if

        ! Inspect content
        ! Get info about stored matrix
        inquire(iolength=reclen) iqmt_, wl_, wu_, econv_, nk_bse_, ikkp_,&
          & iknr_, jknr_, inu_, ino_, jnu_, jno_, ikoulims_(4), jkoulims_(4)
        call getunit(un)
        open(unit=un, file=trim(fname), form='unformatted',&
          & action='read', access='direct', recl=reclen)
        read(un, rec=1) iqmt_, wl_, wu_, econv_, ikkp_, iknr_, jknr_,&
          & inu_, ino_, jnu_, jno_, ikoulims_(4), jkoulims_(4)
        close(un)

        ! Check if stored data is compatible with requested data
        if(iqmt /= iqmt_) then 
          write(*,*)
          write(*, '("Error (b_getbsemat):&
            & Requested momentum transfer index differs")')
          write(*, '(" requested iqmt: ", i4)') iqmt 
          write(*, '(" stored iqmt_: ", i4)') iqmt_
          write(*,*)
          call terminate
        end if
        if(wl-econv < wl_-econv_ .or. wu+econv > wu_+econv_) then 
          write(*,*)
          write(*, '("Error (b_getbsemat): Requested energy range incompatible")')
          write(*, '(" requested (wl, wu, econv): ", 3E12.3)') wl, wu, econv
          write(*, '(" stored (wl_, wu_, econv_): ", 3E12.3)') wl_, wu_, econv_
          write(*,*)
          call terminate
        end if
        if(nk_bse /= nk_bse_) then 
          write(*,*)
          write(*, '("Error (b_getbsemat): Differring number of relevant k-points")')
          write(*, '(" requested nk_bse: ", i4)') nk_bse 
          write(*, '(" stored nk_bse_: ", i4)') nk_bse_
          write(*,*)
          call terminate
        end if
        if(ikkp /= ikkp_ .or. iknr /= iknr_ .or. jknr /= jknr_) then 
          write(*,*)
          write(*, '("Error (b_getbsemat): Requested k-points differ")')
          write(*, '(" requested (ikkp, iknr, jknr): ", 3i8)') ikkp, iknr, jknr
          write(*, '(" stored (ikkp_, iknr_, jknr_): ", 3i8)') ikkp_, iknr_, jknr_
          write(*,*)
          call terminate
        end if
        if(inu > inu_ .or. ino > ino_ .or. jnu > jnu_ .or. jno > jno_) then
          write(*,*)
          write(*, '("Error (b_getbsemat): Requested matrix size out of range")')
          write(*, '(" requested size (inu,ino,jnu,jno)  : ", 4i8)')&
            & inu, ino, jnu, jno
          write(*, '(" stored size (inu_,ino_,jnu_,jno_) : ", 4i8)')&
            & inu_, ino_, jnu_, jno_
          write(*,*)
          call terminate
        end if
        if(koulims(1,iknr) < ikoulims_(1) .or. koulims(2,iknr) > ikoulims_(2)&
          & .or. koulims(3,iknr) < ikoulims_(3) .or. koulims(4,iknr) < ikoulims_(4)) then
          write(*,*)
          write(*, '("Error (b_getbsemat): Requested bands out of range")')
          write(*, '(" requested size (iu_min,iu_max,io_min,io_max)  : ", 4i8)')&
            & koulims(:,iknr)
          write(*, '(" stored size (iu_min_,iu_max_,io_min_,io_max_) : ", 4i8)')&
            & ikoulims_(:)
          write(*,*)
          call terminate
        end if
        if(koulims(1,jknr) < jkoulims_(1) .or. koulims(2,jknr) > jkoulims_(2)&
          & .or. koulims(3,jknr) < jkoulims_(3) .or. koulims(4,jknr) < jkoulims_(4)) then
          write(*,*)
          write(*, '("Error (b_getbsemat): Requested bands out of range")')
          write(*, '(" requested size (ju_min,ju_max,jo_min,jo_max)  : ", 4i8)')&
            & koulims(:,jknr)
          write(*, '(" stored size (ju_min_,ju_max_,jo_min_,jo_max_) : ", 4i8)')&
            & jkoulims_(:)
          write(*,*)
          call terminate
        end if

        ! Requested matrix is just a part of the stored one.
        if(inu /= inu_ .or. ino /= ino_ .or. jnu /= jnu_ .or. jno /= jno_) then

          ! Get full record length
          inquire(file=trim(fname), recl=reclen)

          ! Read full matrix matrix
          allocate(zm(inu_, ino_, jnu_, jno_))

          call getunit(un)
          open(unit=un, file=trim(fname), form='unformatted',&
            & action='read', access='direct', recl=reclen)
          
          read(un, rec=ikkp) iqmt_, wl_, wu_, econv_,&
            & ikkp_, iknr_, jknr_, inu_, ino_, jnu_, jno_, ikoulims_(4), jkoulims_(4),&
            & zm

          close(un)

          ! Get offsets
          inu_offset = koulims(1, iknr) - ikoulims_(1)
          ino_offset = koulims(3, iknr) - ikoulims_(3)
          jnu_offset = koulims(1, jknr) - jkoulims_(1)
          jno_offset = koulims(3, jknr) - jkoulims_(3)

          ! Copy requested data into outgoing array
          zmat(1:inou, 1:jnou) =&
            & reshape(zm(1+inu_offset:inu_offset+inu, 1+ino_offset:ino_offset+ino,&
            &            1+jnu_offset:jnu_offset+jnu, 1+jno_offset:jno_offset+jno),&
            &         [inou, jnou])

          deallocate(zm)

          return

        end if

      end if

      ! Requested matrix is (or is assumed to be) identical to full stored one

      ! Get full record length
      inquire(file=trim(fname), recl=reclen)

      call getunit(un)
      open(unit=un, file=trim(fname), form='unformatted',&
        & action='read', access='direct', recl=reclen)
      
      read(un, rec=ikkp) iqmt_, wl_, wu_, econv_,&
        & ikkp_, iknr_, jknr_, inu_, ino_, jnu_, jno_, ikoulims_(4), jkoulims_(4),&
        & zmat

      close(un)

    end subroutine b_getbsemat
    !EOC

end module m_putgetbsemat
