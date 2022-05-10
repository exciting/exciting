module putgeteps0
  use modmpi
  use m_getunit, only: getunit
  use m_genfilname
  use precision, only: sp, dp

  implicit none

  private
  
  public :: geteps0_finite_q, geteps0_zero_q,&
            puteps0_finite_q, puteps0_zero_q

  contains

    !> Writes the microscopic V-symmetrized Kohn-Sham dielectric function/tensor
    !>   \[
    !>     \tilde{\epsilon}^1_{\bf{GG'}}({\bf q},\omega) = 
    !>     \delta_{\bf{GG'}}- \tilde{\chi}^0_{\bf{GG'}}({\bf q},\omega)
    !>                                                                  /]
    !> for \( \mathbf{q} =  0 \)for given frequency to a direct access file.
    !> The record index is equal to frequency index,
    !> for each q point a new file is used.
    !> In each record the following is saved:
    !> iq, qvec, numgq, iw, w, eps0(, eps0wg, chhd).
    !> Head and wings are present due to vanishing \( \mathbf{q} \). 
    subroutine puteps0_zero_q(iq, qvec, iw, w, eps0, eps0wg, eps0hd, &
                            fname, debug)
      use modmpi, only: terminate_if_false

      !> q-point index from reduced q-point set
      integer(sp), intent(in) :: iq
      !> q-point in cartesian coordinates
      real(dp), intent(in) :: qvec(3) 
      !> Frequency index
      integer(sp), intent(in) ::  iw
      !> Frequency
      real(dp), intent(in) :: w
      !> Body of RPA epsilon matrix 
      complex(dp), intent(in) :: eps0(:, :)
      !> Wings of RPA epsilon matrix 
      complex(dp), intent(in) :: eps0wg(:,:,:)
      !> Head of dielectric matrix
      complex(dp), intent(in) :: eps0hd(:,:)
      !> Filename
      character(*), intent(in), optional :: fname
      !> Debug mode
      logical, intent(in), optional :: debug
      
      !> Name of subroutine
      character(*), parameter :: thisnam = 'puteps0_zero_q'
      !> Filename
      character(256) :: filename
      !> File unit
      integer(sp) :: un
      !> I/O status
      integer(sp) ::  stat
      !> Record length 
      integer(sp) ::  reclen
      !> Number of (G+q)-points for current q-point
      integer(sp) :: numgq
      !> Local debug mode
      logical :: debug_local

      numgq = size(eps0, dim=1)

      ! Generate q dependent file name
      if(present(fname)) then 
        filename=trim(adjustl(fname))
      else
        call genfilname(basename='EPS0', iq=iq, filnam=filename)
      end if

      ! Check for debug mode
      if(present(debug)) then
        debug_local = debug
      else 
        debug_local = .false.
      end if

      call getunit(un)
      
      ! Get record length (q dependent through numgq)
      inquire (iolength=reclen) iq, qvec, numgq, iw, w,&
            & eps0hd, eps0wg, eps0

      ! Open with recl=1, so that record length is determined 
      open(unit=un, file=trim(filename), form='unformatted',&
        & action='write', access='direct', recl=reclen, iostat=stat)
      
        if(stat /= 0) then
        write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
          & "Error opening file, iostat:", stat
        call terminate
      end if

      ! Write output for frequency iw
      write(un, rec=iw, iostat=stat) iq, qvec, numgq, iw, w,&
        & eps0hd, eps0wg, eps0

      if(stat /= 0) then
        write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
          & "Error writing file, iostat:", stat
        call terminate
      end if

      close(un)

      ! Debug output
      if(debug_local) then
        write(*,*) "========puteps0_zero_q========"
        write(*,*) "Writing file:", trim(filename)
        write(*,*) "Writing record:", iw
        write(*,*) "Record of length:", reclen
        write(*,*) "iq =",iq
        write(*,*) "vqc =", qvec
        write(*,*) "numgq =", numgq
        write(*,*) "iw =", iw
        write(*,*) "w =", w
        write(*,*) "shape(eps0hd) =", shape(eps0hd)
        write(*,*) "shape(eps0wg) =", shape(eps0wg)
        write(*,*) "shape(eps0) =", shape(eps0)
        write(*,*) "======================="
      end if

    end subroutine puteps0_zero_q


    !> Writes the microscopic V-symmetrized Kohn-Sham dielectric function/tensor
    !>   \[
    !>     \tilde{\epsilon}^1_{\bf{GG'}}({\bf q},\omega) = 
    !>     \delta_{\bf{GG'}}- \tilde{\chi}^0_{\bf{GG'}}({\bf q},\omega)
    !>                                                                  /]
    !> for \( \mathbf{q} \neq 0 \)for given frequency to a direct access file.
    !> The record index is equal to frequency index,
    !> for each q point a new file is used.
    !> In each record the following is saved:
    !> iq, qvec, numgq, iw, w, eps0(, eps0wg, chhd).
    !> Head and wings are not present due to finite \( \mathbf{q} \). 
    subroutine puteps0_finite_q(iq, qvec, iw, w, eps0, &
                                fname, debug)
      use modmpi, only: terminate_if_false

      !> q-point index from reduced q-point set
      integer(sp), intent(in) :: iq
      !> q-point in cartesian coordinates
      real(dp), intent(in) :: qvec(3) 
      !> Frequency index
      integer(sp), intent(in) ::  iw
      !> Frequency
      real(dp), intent(in) :: w
      !> Body of RPA epsilon matrix 
      complex(dp), intent(in) :: eps0(:, :)
      !> Filename
      character(*), intent(in), optional :: fname
      !> Debug mode
      logical, intent(in), optional :: debug
      
      !> Name of subroutine
      character(*), parameter :: thisnam = 'puteps0_finite_q'
      !> Filename
      character(256) :: filename
      !> File unit
      integer(sp) :: un
      !> I/O status
      integer(sp) ::  stat
      !> Record length 
      integer(sp) ::  reclen
      !> Number of (G+q)-points for current q-point
      integer(sp) :: numgq
      !> Local debug mode
      logical :: debug_local

      numgq = size(eps0, dim=1)

      ! Generate q dependent file name
      if(present(fname)) then 
        filename=trim(adjustl(fname))
      else
        call genfilname(basename='EPS0', iq=iq, filnam=filename)
      end if

      ! Check for debug mode
      if(present(debug)) then
        debug_local = debug
      else 
        debug_local = .false.
      end if

      call getunit(un)
  
      ! Get record length (q dependent through numgq)
      inquire (iolength=reclen) iq, qvec, numgq, iw, w, eps0

      ! Open with recl=1, so that record length is determined 
      open(unit=un, file=trim(filename), form='unformatted',&
        & action='write', access='direct', recl=reclen, iostat=stat)
      if(stat /= 0) then
        write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
          & "Error opening file, iostat:", stat
        call terminate
      end if

      ! Write output for frequency iw
      write(un, rec=iw, iostat=stat) iq, qvec, numgq, iw, w, eps0
      if(stat /= 0) then
        write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
          & "Error writing file, iostat:", stat
        call terminate
      end if

      close(un)

      ! Debug output
      if(debug_local) then
        write(*,*) "========puteps0_finite_q========"
        write(*,*) "Writing file:", trim(filename)
        write(*,*) "Writing record:", iw
        write(*,*) "Record of length:", reclen
        write(*,*) "iq =",iq
        write(*,*) "vqc =", qvec
        write(*,*) "numgq =", numgq
        write(*,*) "iw =", iw
        write(*,*) "w =", w
        write(*,*) "shape(eps0) =", shape(eps0)
        write(*,*) "======================="
      end if
    end subroutine puteps0_finite_q

    !>   Reads the microscopic V-symmetrized Kohn-Sham dielectric function/tensor
    !> \[
    !>   \tilde{\epsilon}^0_{\bf{GG'}}({\bf q},\omega) = 
    !>    delta_{\bf{GG'}} - \tilde{\chi}^0_{\bf{GG'}}({\bf q},\omega)
    !>                                                                 \]
    !> for \( \mathbf{q} = 0 \) for given frequency from a direct access file.
    !> The record index is equal to frequency index
    !> In each record the following is saved:
    !> iq, qvec, numgq, iw, w, eps0(, eps0wg, chhd).
    !> Head and wings are present due to vanishing \( \mathbf{q} \).
    subroutine geteps0_zero_q(iq, qvec, iw, w, eps0, eps0wg, eps0hd, &
                           fname, debug)
      use math_utils, only: all_close, all_zero
      
      !> q-point index from reduced q-vector set if xs%reduceq = "true",
      !> else from non-reduced set
      integer(sp), intent(in) :: iq
      !> q-point in cartesian coordinates
      real(dp), intent(in) :: qvec(3) 
      !> Frequency index
      integer(sp), intent(in) ::  iw
      !> Frequency
      real(dp), intent(in) :: w
      !> Body of RPA epsilon matrix 
      complex(dp), intent(out) :: eps0(:, :)
      !> Wings of RPA epsilon matrix 
      complex(dp), intent(out) :: eps0wg(:,:,:)
      !> Head of dielectric matrix
      complex(dp), intent(out) :: eps0hd(:,:)
      !> Filename
      character(*), intent(in), optional :: fname
      !> Debug mode
      logical, intent(in), optional :: debug
      !> Number of (G+q)-points for current q-point
      integer(sp) :: numgq
      !> Routine name
      character(*), parameter :: thisnam = 'geteps0_zero_q'
      !> Filename
      character(256):: filename
      !> File unit
      integer(sp) :: un
      !> I/O status
      integer(sp) ::  stat
      !> Record length 
      integer(sp) ::  reclen      !> Number of (G+q)-points for current q-point as read from file
      integer(sp) :: numgq_read
      !> q-point index as read from file
      integer(sp) :: iq_read
      !> Frequency index  as read from file
      integer(sp) :: iw_read
      !> q-point in cartesian coordinates as read from file
      real(dp) :: qvec_read(3)
      !> Frequency as read from file
      real(dp) ::  w_read
      !> Tolerance for  comparing of requested and read data
      real(dp), parameter :: tol=1.e-8_dp
      logical :: existent
      !> Local debug mode
      logical :: debug_local


      ! Generate q dependent file name
      if(present(fname)) then 
        filename=trim(adjustl(fname))
      else
        call genfilname(basename='EPS0', iq=iq, filnam=filename)
      end if
      call getunit(un)

      ! Check for debug mode
      if(present(debug)) then
        debug_local = debug
      else 
        debug_local = .false.
      end if

      ! Check if file exists
      inquire(file=trim(filename), exist=existent)
      if( .not. existent) then
        write(*, '(a)') 'Error(' // trim(thisnam) // '):&
          & file does not exist:' // trim(filename)
        call terminate
      end if

      numgq = size(eps0, dim=1)

      ! Get q-dependent record length
      inquire(iolength=reclen) iq_read, qvec_read, numgq_read, &
                      iw_read, w_read, eps0hd, eps0wg, eps0

      open(unit=un, file=trim(filename), status='old',&
        & form='unformatted', action='read',&
        & access='direct', recl=reclen, iostat=stat)
      if(stat /= 0) then
        write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
          & "Error opening file, iostat:", stat
        call terminate
      end if

      ! Read from recpos onwards
      read(un, rec=iw, iostat=stat) iq_read, qvec_read, numgq_read,&
                    iw_read, w_read, eps0hd, eps0wg, eps0
      if(stat /= 0) then
        write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
          & "Error reading file, iostat:", stat
        call terminate
      end if

      close(un)

       ! Debug output
      if(debug_local) then
        write(*,*) "========geteps0_zero_q========"
        write(*,*) "Reading file:", trim(filename)
        write(*,*) "Reading record:", iw
        write(*,*) "Reading record of length:", reclen
        write(*,*) "iq =",iq
        write(*,*) "vqc =", qvec
        write(*,*) "numgq =", numgq
        write(*,*) "iw =", iw
        write(*,*) "w =", w
        write(*,*) "shape(eps0hd) =", shape(eps0hd)
        write(*,*) "shape(eps0wg) =", shape(eps0wg)
        write(*,*) "shape(eps0) =", shape(eps0)
        write(*,*) "======================="
      end if

      ! Check consistency of requested data with saved data 
      if (numgq_read /= numgq &
        .or. .not. all_close(qvec_read, qvec, tol)&
        .or. .not. all_close(w_read, w, tol)&
        .or. iq_read /= iq &
        .or. iw_read /= iw) then

          write (*, '(a)') 'Error('//trim(thisnam)//'):&
            & differing parameters for matrix elements (current/file): '
          write (*, '(a, 2i6)') 'numgq', numgq, numgq_read
          write (*, '(a, 3f12.6, a, 3f12.6)') 'qvec', qvec, ', ', qvec_read
          write (*, '(a, 2i6)') 'for q-point :', iq, iq_read
          write (*, '(a, 2i6)') 'for frequency-point :', iw, iw_read
          write (*, '(a, 2E13.6)') 'for frequency', w, w_read
          write (*, '(a)') ' file: ', trim(filename)
          call terminate

      end if

    end subroutine geteps0_zero_q


    !>   Reads the microscopic V-symmetrized Kohn-Sham dielectric function/tensor
    !> \[
    !>   \tilde{\epsilon}^0_{\bf{GG'}}({\bf q},\omega) = 
    !>    delta_{\bf{GG'}} - \tilde{\chi}^0_{\bf{GG'}}({\bf q},\omega)
    !>                                                                 \]
    !> for \( \mathbf{q} = 0 \) for given frequency from a direct access file.
    !> The record index is equal to frequency index
    !> In each record the following is saved:
    !> iq, qvec, numgq, iw, w, eps0(, eps0wg, chhd).
    !> Head and wings are not present due to finite \( \mathbf{q} \).
    subroutine geteps0_finite_q(iq, qvec, iw, w, eps0, &
                                fname, debug)
   
      use math_utils, only: all_close, all_zero

      !> q-point index from reduced q-vector set if xs%reduceq = "true",
      !> else from non-reduced set
      integer(sp), intent(in) :: iq
      !> q-point in cartesian coordinates
      real(dp), intent(in) :: qvec(3) 
      !> Frequency index
      integer(sp), intent(in) ::  iw
      !> Frequency
      real(dp), intent(in) :: w
      !> Body of RPA epsilon matrix 
      complex(dp), intent(out) :: eps0(:, :)
      !> Filename
      character(*), intent(in), optional :: fname
      !> Debug mode
      logical, intent(in), optional :: debug

      !> Number of (G+q)-points for current q-point
      integer(sp) :: numgq
      !> Routine name
      character(*), parameter :: thisnam = 'geteps0_finite_q'
      !> Filename
      character(256):: filename
      !> File unit
      integer(sp) :: un
      !> I/O status
      integer(sp) ::  stat
      !> Record length 
      integer(sp) ::  reclen
      !> Number of (G+q)-points for current q-point as read from file
      integer(sp) :: numgq_read
      !> q-point index as read from file
      integer(sp) :: iq_read
      !> Frequency index  as read from file
      integer(sp) :: iw_read
      !> q-point in cartesian coordinates as read from file
      real(dp) :: qvec_read(3)
      !> Frequency as read from file
      real(dp) ::  w_read
      !> Tolerance for  comparing of requested and read data
      real(dp), parameter :: tol=1.e-8_dp
      logical :: existent
      !> Local debug mode
      logical :: debug_local

      ! Generate q dependent file name
      if(present(fname)) then 
        filename=trim(adjustl(fname))
      else
        call genfilname(basename='EPS0', iq=iq, filnam=filename)
      end if

      ! Check for debug mode
      if(present(debug)) then
        debug_local = debug
      else 
        debug_local = .false.
      end if


      call getunit(un)

      ! Check if file exists
      inquire(file=trim(filename), exist=existent)
      if( .not. existent) then
        write(*, '(a)') 'Error(' // trim(thisnam) // '):&
          & file does not exist:' // trim(filename)
        call terminate
      end if

      numgq = size(eps0, dim=1)

      ! Get q-dependent record length
      inquire(iolength=reclen) iq_read, qvec_read, numgq_read, &
                      iw_read, w_read, eps0

      open(unit=un, file=trim(filename), status='old',&
        & form='unformatted', action='read',&
        & access='direct', recl=reclen, iostat=stat)
      if(stat /= 0) then
        write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
          & "Error opening file, iostat:", stat
        call terminate
      end if

      ! Read from recpos onwards
      read(un, rec=iw, iostat=stat) iq_read, qvec_read, numgq_read,&
                    iw_read, w_read, eps0
      if(stat /= 0) then
        write(*,'("Error: (",a,"):",a, i3)') trim(thisnam),&
          & "Error reading file, iostat:", stat
        call terminate
      end if

      close(un)

      ! Debug output
      if(debug_local) then
        write(*,*) "========geteps0_finite_q========"
        write(*,*) "Reading file:", trim(filename)
        write(*,*) "Reading record:", iw
        write(*,*) "Reading record of length:", reclen
        write(*,*) "iq =",iq
        write(*,*) "vql =", qvec
        write(*,*) "numgq =", numgq
        write(*,*) "iw =", iw
        write(*,*) "w =", w
        write(*,*) "shape(eps0) =", shape(eps0)
        write(*,*) "======================="
      end if

      ! Check consistency of requested data with saved data 
      if (numgq_read /= numgq &
        .or. .not. all_close(qvec_read, qvec, tol)&
        .or. .not. all_close(w_read, w, tol)&
        .or. iq_read /= iq &
        .or. iw_read /= iw) then

          write (*, '(a)') 'Error('//trim(thisnam)//'):&
            & differing parameters for matrix elements (current/file): '
          write (*, '(a, 2i6)') 'numgq', numgq, numgq_read
          write (*, '(a, 3f12.6, a, 3f12.6)') 'qvec', qvec, ', ', qvec_read
          write (*, '(a, 2i6)') 'for q-point :', iq, iq_read
          write (*, '(a, 2i6)') 'for frequency-point :', iw, iw_read
          write (*, '(a, 2E13.6)') 'for frequency', w, w_read
          write (*, '(a)') ' file: ', trim(filename)
          call terminate

      end if

    end subroutine geteps0_finite_q

end module putgeteps0
