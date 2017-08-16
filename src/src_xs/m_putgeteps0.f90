module m_putgeteps0
  use modmpi
  use modinput, only: input
  ! Where to find the q point information
  use mod_qpoint, only: nqpt, vql, ivq
  use modxs, only: ngq, nqptr, vqlr, ivqr, ngqr, nwdf, eps0dirname
  ! I/O modules
  use m_getunit
  use m_genfilname

  implicit none

  private :: zerolength, qwidx, qwidx_back

  contains

    !BOP
    ! !ROUTINE: puteps0
    ! !INTERFACE:
    subroutine puteps0(reduced, iq, iw, w, eps0, eps0wg, eps0hd, fname)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! logical :: reduced          ! q-point from reduced set
    ! integer(4) :: iq            ! q-point index
    ! integer(4) :: iw            ! w-point index
    ! real(8) :: w                ! iw'th frequency
    ! complex(8) :: eps0(:,:)     ! Body of RPA epsilon matrix 
    ! complex(8) :: eps0wg(:,:,:) ! Wings of RPA epsilon matrix 
    ! complex(8), optional :: eps0hd(:,:)   ! Body of RPA epsilon matrix
    !
    ! !DESCRIPTION:
    !   Writes the microscopic V-symmetrized Kohn-Sham dielectric function/tensor
    !   $\tilde{\epsilon}^1_{\bf{GG'}}({\bf q},\omega) = \delta_{\bf{GG'}}
    !   - \tilde{\chi}^0_{\bf{GG'}}({\bf q},\omega)$ for given frequency and q point
    !   to a direct access file.
    !   The record index is equal to frequency index,
    !   for each q point a new file is used.
    !   In each record the following is saved: iq, vql, ngq, iw, w, eps0(, eps0wg, chhd).
    !   Head and wings only present for $q=0$.
    !
    ! !REVISION HISTORY:
    !   Created. (Aurich)
    !EOP
    !BOC

      implicit none

      ! Arguments
      logical, intent(in) :: reduced
      integer(4), intent(in) :: iq, iw
      real(8), intent(in) :: w
      complex(8), intent(in) :: eps0(:, :)
      complex(8), intent(in), optional :: eps0wg(:,:,:), eps0hd(:,:)
      character(*), intent(in), optional :: fname

      ! Local variables
      character(*), parameter :: thisnam = 'puteps0'
      character(256) :: filename
      integer(4) :: un, stat, reclen
      logical :: tq0
      integer(4), pointer :: ivq_p(:,:), nqpt_p, ngq_p(:)
      real(8), pointer :: vql_p(:,:)

      ! Check if iq references the reduced or non reduced set
      if(reduced) then
        ivq_p => ivqr
        vql_p => vqlr
        nqpt_p => nqptr
        ngq_p => ngqr
      else
        ivq_p => ivq
        vql_p => vql
        nqpt_p => nqpt
        ngq_p => ngq
      end if

      ! Check if q=0
      tq0 = zerolength(vql_p(:,iq))

      ! Check q=0 but head or wings missing
      if(tq0 .and. (.not. present(eps0wg) .or.  .not. present(eps0wg))) then
        write(*,*) 'Error(' // trim(thisnam) // '): q=0 but head or wings missing'
        call terminate
      end if

      ! Generate q dependend file name
      if(present(fname)) then 
        filename=trim(adjustl(fname))
      else
        call genfilname(basename='EPS0', iq=iq, filnam=filename)
      end if
      call getunit(un)

      if(tq0) then

        ! Get record length (q dependend trough ngq)
        inquire(iolength=reclen) iq, vql_p(:,iq), ngq_p(iq), iw, w,&
          & eps0hd, eps0wg, eps0

        ! Open with recl=1, so that record length is determined 
        open(unit=un, file=trim(filename), form='unformatted',&
          & action='write', access='direct', recl=reclen, iostat=stat)
        if(stat /= 0) then
          write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
            & "Error opening file, iostat:", stat
          call terminate
        end if

        ! Debug output
        if(input%xs%dbglev>2) then
          write(*,*) "========PUTEPS0========"
          write(*,*) "Writing file:", trim(filename)
          write(*,*) "Writing record:", iw
          write(*,*) "Record of length:", reclen
          write(*,*) "iq =",iq
          write(*,*) "vql =", vql_p(:,iq)
          write(*,*) "ngq =", ngq_p(iq)
          write(*,*) "iw =", iw
          write(*,*) "w =", w
          write(*,*) "shape(eps0hd) =", shape(eps0hd)
          write(*,*) "shape(eps0wg) =", shape(eps0wg)
          write(*,*) "shape(eps0) =", shape(eps0)
          write(*,*) "======================="
        end if

        ! Write output for frequency iw
        write(un, rec=iw, iostat=stat) iq, vql_p(:,iq), ngq_p(iq), iw, w,&
          & eps0hd, eps0wg, eps0
        if(stat /= 0) then
          write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
            & "Error writing file, iostat:", stat
          call terminate
        end if

      else ! not q=0

        inquire(iolength=reclen) iq, vql_p(:,iq), ngq_p(iq), iw, w,&
          & eps0hd, eps0wg, eps0

        open(unit=un, file=trim(filename), form='unformatted',&
          & action='write', access='direct', recl=reclen, iostat=stat)
        if(stat /= 0) then
          write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
            & "Error opening file, iostat:", stat
          call terminate
        end if

        ! Debug output
        if(input%xs%dbglev>2) then
          write(*,*) "========PUTEPS0========"
          write(*,*) "Writing file:", trim(filename)
          write(*,*) "Writing record:", iw
          write(*,*) "Record of length:", reclen
          write(*,*) "iq =",iq
          write(*,*) "vql =", vql_p(:,iq)
          write(*,*) "ngq =", ngq_p(iq)
          write(*,*) "iw =", iw
          write(*,*) "w =", w
          write(*,*) "shape(eps0hd) =", shape(eps0hd)
          write(*,*) "shape(eps0wg) =", shape(eps0wg)
          write(*,*) "shape(eps0) =", shape(eps0)
          write(*,*) "======================="
        end if

        write(un, rec=iw, iostat=stat) iq, vql_p(:,iq), ngq_p(iq), iw, w, eps0
        if(stat /= 0) then
          write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
            & "Error writing file, iostat:", stat
          call terminate
        end if

      end if

      close(un)

    end subroutine puteps0
    !EOC

    !BOP
    ! !ROUTINE: geteps0
    ! !INTERFACE:
    subroutine geteps0(reduced, iq, iw, w, eps0, eps0wg, eps0hd, fname)
    ! !USES:
      use m_getunit
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! logical :: reduced          ! q-point from reduced set
    ! integer(4) :: iq            ! q-point index
    ! integer(4) :: iw            ! w-point index
    ! real(8) :: w                ! iw'th frequency
    ! Out:
    ! complex(8) :: eps0(:,:)               ! Body of RPA epsilon matrix 
    ! complex(8), optional :: eps0wg(:,:,:) ! Wings of RPA epsilon matrix 
    ! complex(8), optional :: eps0hd(:,:)   ! Body of RPA epsilon matrix
    !
    ! !DESCRIPTION:
    !   Read the microscopic V-symmetrized Kohn-Sham dielectric function/tensor
    !   $\tilde{\epsilon}^0_{\bf{GG'}}({\bf q},\omega) = \delta_{\bf{GG'}}
    !   - \tilde{\chi}^0_{\bf{GG'}}({\bf q},\omega)$ for given frequency and q point
    !   from a direct access file.
    !   The record index is equal to frequency index
    !   In each record the following is saved: iq, vql, ngq, iw, w, eps0(, eps0wg, chhd).
    !   Head and wings only present for $q=0$.
    !
    ! !REVISION HISTORY:
    !   Created. (Aurich)
    !EOP
    !BOC

      implicit none

      ! Arguments
      logical, intent(in) :: reduced
      integer(4), intent(in) :: iq, iw
      real(8), intent(in) :: w
      complex(8), intent(out) :: eps0(:,:)
      complex(8), intent(out), optional :: eps0wg(:,:,:), eps0hd(:,:)
      character(*), intent(in), optional :: fname

      ! Local variables
      character(*), parameter :: thisnam = 'geteps0'
      character(256):: filename
      integer(4) :: un, stat, reclen
      integer(4) :: ngq_, iq_, iw_
      real(8) :: vql_(3), w_
      real(8), parameter :: epslat=1.0d-8
      logical :: tq0, existent
      integer(4), pointer :: ivq_p(:,:), nqpt_p, ngq_p(:)
      real(8), pointer :: vql_p(:,:)

      ! Generate q dependend file name
      if(present(fname)) then 
        filename=trim(adjustl(fname))
      else
        call genfilname(basename='EPS0', iq=iq, filnam=filename)
      end if
      call getunit(un)

      ! Check if file exists
      inquire(file=trim(filename), exist=existent)
      if( .not. existent) then
        write(*, '(a)') 'Error(' // trim(thisnam) // '):&
          & file does not exist:' // trim(filename)
        call terminate
      end if

      ! Check if iq references the reduced or non reduced set
      if(reduced) then
        ivq_p => ivqr
        vql_p => vqlr
        nqpt_p => nqptr
        ngq_p => ngqr
      else
        ivq_p => ivq
        vql_p => vql
        nqpt_p => nqpt
        ngq_p => ngq
      end if

      ! Check if q=0
      tq0 = zerolength(vql_p(:,iq))

      ! Check if q=0 but head or wings missing
      if(tq0 .and. (.not. present(eps0wg) .or. .not. present(eps0wg)) ) then
        write(*,*) 'Error(' // trim(thisnam) // '): q=0 but head or wings missing'
        call terminate
      end if

      if(tq0) then

        ! Get q-dependend record length
        inquire(iolength=reclen) iq_, vql_, ngq_, iw_, w_, eps0hd, eps0wg, eps0

        open(unit=un, file=trim(filename), status='old',&
          & form='unformatted', action='read',&
          & access='direct', recl=reclen, iostat=stat)
        if(stat /= 0) then
          write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
            & "Error opening file, iostat:", stat
          call terminate
        end if

        ! Debug output
        if(input%xs%dbglev>2) then
          write(*,*) "========GETEPS0========"
          write(*,*) "Reading file:", trim(filename)
          write(*,*) "Reading record:", iw
          write(*,*) "Reading record of length:", reclen
          write(*,*) "iq =",iq
          write(*,*) "vql =", vql_p(:,iq)
          write(*,*) "ngq =", ngq_p(iq)
          write(*,*) "iw =", iw
          write(*,*) "w =", w
          write(*,*) "shape(eps0hd) =", shape(eps0hd)
          write(*,*) "shape(eps0wg) =", shape(eps0wg)
          write(*,*) "shape(eps0) =", shape(eps0)
          write(*,*) "======================="
        end if

        ! Read from recpos onwards
        read(un, rec=iw, iostat=stat) iq_, vql_, ngq_, iw_, w_, eps0hd, eps0wg, eps0
        if(stat /= 0) then
          write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
            & "Error reading file, iostat:", stat
          call terminate
        end if

      else

        inquire(iolength=reclen) iq_, vql_, ngq_, iw_, w_, eps0

        open(unit=un, file=trim(filename), status='old',&
          & form='unformatted', action='read',&
          & access='direct', recl=reclen, iostat=stat)
        if(stat /= 0) then
          write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
            & "Error opening file, iostat:", stat
          call terminate
        end if

        ! Debug output
        if(input%xs%dbglev>2) then
          write(*,*) "========GETEPS0========"
          write(*,*) "Reading file:", trim(filename)
          write(*,*) "Reading record:", iw
          write(*,*) "Reading record of length:", reclen
          write(*,*) "iq =",iq
          write(*,*) "vql =", vql_p(:,iq)
          write(*,*) "ngq =", ngq_p(iq)
          write(*,*) "iw =", iw
          write(*,*) "w =", w
          write(*,*) "shape(eps0) =", shape(eps0)
          write(*,*) "======================="
        end if

        read(un, rec=iw, iostat=stat) iq_, vql_, ngq_, iw_, w_, eps0
        if(stat /= 0) then
          write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
            & "Error reading file, iostat:", stat
          call terminate
        end if

      end if

      close(un)

      ! Check consistency of requested data with saved data 
      if(ngq_ .ne. ngq_p(iq)&
        & .or. any(abs(vql_-vql_p(:,iq)) > epslat)&
        & .or. w_ .ne. w) then

        write(*, '(a)') 'Error(' // trim(thisnam) // '):&
          & differring parameters for matrix elements (current/file): '
        write(*, '(a, 2i6)') 'ngq', ngq_p(iq), ngq_
        write(*, '(a, 3f12.6, a, 3f12.6)') 'vql', vql_p(:,iq), ', ', vql_
        write(*, '(a, 2i6)') 'for q-point :', iq, iq_
        write(*, '(a, 2i6)') 'for w-point :', iw, iw_
        write(*, '(a, 2E13.6)') 'w', w, w_
        write(*, '(a)') ' file: ', trim(filename) 
        call terminate

      end if

    end subroutine geteps0
    !EOC

    ! Auxilliary routines
    integer(4) function qwidx(iq, iw, nq)
      implicit none
      integer(4), intent(in) :: iq, iw, nq
      qwidx = iq + (iw-1) * nq
    end function qwidx

    subroutine qwidx_back(iqw, nq, iq, iw)
      implicit none
      integer(4), intent(in) :: iqw, nq
      integer(4), intent(out) :: iq, iw
      iw = (iqw-1)/nq + 1
      iq = iqw - (iw-1)*nq
    end subroutine qwidx_back

    logical function zerolength(vq)
      real(8), intent(in) :: vq(3)
      real(8), parameter :: epsg = 1.d-12
      zerolength = .false.
      if (sum(abs(vq)) .lt. epsg) zerolength = .true.
    end function zerolength

end module m_putgeteps0
