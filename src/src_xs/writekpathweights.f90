subroutine writekpathweights
  use modxs, only: unitout
  use modinput
  use m_read_bandstructure
  use m_putgetexcitons
  use mod_constants, only: h2ev
  use bspline_module
  use, intrinsic :: iso_fortran_env, only: wp => real64

  implicit none

  character(*), parameter :: thisname = 'writekpathweights'
  character(256) :: exckpathdir, syscommand

  logical :: fwritegridweights
  integer(4) :: iqmt, iqmti, iqmtf, nqmt, iq1, iq2
  real(8) :: en1, en2
  integer(4) :: i1, i2

  integer(4) :: nvmax, ncmax, icmin, icmax, ivmin, ivmax
  integer(4) :: iv1, iv2, ic1, ic2, ib1, ib2, nsteps
  integer(4) :: lambda
  integer(4) :: intorder

  ! Interpolation arrays
  real(wp), allocatable :: x(:), y(:), z(:)
  real(wp), allocatable :: x0(:), y0(:), z0(:)
  real(wp), allocatable :: inputdata(:,:,:)

  ! Data on gird
  real(8), allocatable :: abs2(:), rvwgrid(:,:), rcwgrid(:,:),&
    & arvwgrid(:,:), arcwgrid(:,:) 

  ! Interpolated data
  real(8), allocatable :: rvw(:,:), rcw(:,:), arvw(:,:), arcw(:,:)

  if(mpiglobal%rank == 0) then 

    ! Set defaults if writeexcitons is not specified
    if( .not. associated(input%xs%writekpathweights)) then
      write(unitout,'("Error(",a,"):&
        & Specify the writekpathweights element in the input file.")') trim(thisname)
      call terminate
    end if

    !====================================================!
    ! Read in band structure data to read_bandstrucute   !
    ! module                                             !
    !====================================================!
    call read_bandstructure('bandstructure.dat')
    !====================================================!

    ! Make output directory
    exckpathdir='KPATHEXC'
    syscommand = 'test ! -e '//trim(adjustl(exckpathdir))&
      & //' && mkdir -p '//trim(adjustl(exckpathdir))
    call system(trim(adjustl(syscommand)))

    ! Write out excitonic weights on grid?
    fwritegridweights = input%xs%writekpathweights%printgridweights

    ! Q-point list entries
    !   Use all
    iqmti = 1
    iqmtf = size(input%xs%qpointset%qpoint, 2)
    !   or use only one
    if(input%xs%bse%iqmtrange(1) /= -1) then 
      iqmti=input%xs%bse%iqmtrange(1)
      iqmtf=input%xs%bse%iqmtrange(2)
    end if
    nqmt = iqmtf-iqmti+1
    iq1 = 1
    iq2 = nqmt

    do iqmt = iqmti+iq1-1, iqmti+iq2-1

      !====================================================!
      ! Read in data to putgetexcitons module              !
      !====================================================!
      ! Requested excition index range
      if(input%xs%writekpathweights%selectenergy) then 
        en1=input%xs%writekpathweights%minenergyexcitons
        en2=input%xs%writekpathweights%maxenergyexcitons
        if(input%xs%storeexcitons%useev) then 
          en1=en1/h2ev
          en2=en2/h2ev
        end if
        call get_excitons(iqmt=iqmt, e1=en1, e2=en2)
      else
        i1 = input%xs%writekpathweights%minnumberexcitons
        i2 = input%xs%writekpathweights%maxnumberexcitons
        call get_excitons(iqmt=iqmt, a1=i1, a2=i2)
      end if
      !====================================================!

      ! Set maximal number of valence and conduction states over all k-points
      nvmax = maxval(koulims_(4,:)-koulims_(3,:)+1)
      ncmax = maxval(koulims_(2,:)-koulims_(1,:)+1)
      ! Set minimal and maximal valence and conduction state index over all k-points
      icmin = minval(koulims_(1,:))
      icmax = maxval(koulims_(2,:))
      ivmin = minval(koulims_(3,:))
      ivmax = maxval(koulims_(4,:))

      call printline(unitout, "+")
      write(unitout,'("Info(",a,"):&
        & Considering momentum transver vector:", i8)')&
        & trim(thisname), iqmt
      write(unitout,'("Info(",a,"):&
        & Weights can be interpolated in the range:", 2i8)')&
        & trim(thisname), ivmin, icmax

      ! Set index ranges that agree with read in bandstructure and exciton data
      ib1 = brange_(1)
      ib2 = brange_(2)
      nsteps = brange_(3)
      iv1 = max(ivmin, ib1)
      iv2 = min(ivmax, ib2)
      ic1 = max(icmin, ib1)
      ic2 = min(icmax, ib2)

      write(unitout,'("Info(",a,"):&
        & Calculating valence weights in the range:", 2i8)') trim(thisname), iv1, iv2
      write(unitout,'("Info(",a,"):&
        & Calculating conduction weights in the range:", 2i8)') trim(thisname), ic1, ic2
      if(iv2 >= ic1) then 
        write(unitout,'("Error(",a,"):&
          & iv2 >= ic1 Non-Insulators not yet implemented.")') trim(thisname)
        call terminate
      end if
      call printline(unitout, "-")

      ! Squared modulus of coefficients
      allocate(abs2(hamsize_))

      ! Weight data on grid
      allocate(rvwgrid(nk_max_, ivmin:ivmax))
      allocate(rcwgrid(nk_max_, icmin:icmax))
      if(fcoup_) then
        allocate(arvwgrid(nk_max_, ivmin:ivmax))
        allocate(arcwgrid(nk_max_, icmin:icmax))
      end if

      ! Interpolation arrays on -1:1 supercell 
      ! to mimic periodicity
      call setup_xyz(ngridk_, ikmap_, vkl_, x, y, z)
      call setup_xyz(ngridk_, ikmap_, vkl0_, x0, y0, z0)
      allocate(inputdata(3*ngridk_(1),3*ngridk_(2),3*ngridk_(3)))

      ! Interpolated weigts on k-path points
      intorder = input%xs%writekpathweights%intorder
      allocate(rvw(nsteps,iv1:iv2))
      allocate(rcw(nsteps,ic1:ic2))
      if(fcoup_) then 
        allocate(arvw(nsteps,iv1:iv2))
        allocate(arcw(nsteps,ic1:ic2))
      end if

      ! Looping over all selected excitons
      do lambda = iex1_, iex2_

        write(unitout,'("Info(",a,"):&
          & Calculating interpolated weights for excition index:", i8, " /", i8)')&
          & trim(thisname), lambda, iex2_

        !====================================================!
        ! Calculating weigths on the k-grid                  !
        !====================================================!

        ! Zeroing weights
        rvwgrid=0.0d0
        rcwgrid=0.0d0
        if(fcoup_) then 
          arvwgrid=0.0d0
          arcwgrid=0.0d0
        end if

        !! Make resonant weights
        abs2 = abs(rvec_(1:hamsize_,lambda))**2
        !   Valence 
        call genweights(ivmin, ivmax, hamsize_, nk_max_,&
          & smap_(2,:), smap_(3,:), abs2, rvwgrid, ik2ikqmtp_)
        !   Conduction
        call genweights(icmin, icmax, hamsize_, nk_max_,&
          & smap_(1,:), smap_(3,:), abs2, rcwgrid, ik2ikqmtm_)

        !! Make anti-resonant weights
        if(fcoup_) then 
          abs2 = abs(avec_(1:hamsize_, lambda))**2
          !   Valence
          call genweights(ivmin, ivmax, hamsize_, nk_max_,&
            & smap_(2,:), smap_(3,:), abs2, arvwgrid, ik2ikqmtp_)
          !   Conduction
          call genweights(icmin, icmax, hamsize_, nk_max_,&
            & smap_(1,:), smap_(3,:), abs2, arcwgrid, ik2ikqmtm_)
        end if

        if(fwritegridweights) then 
          ! Print weights
          call writeweights()
        end if

        !====================================================!
        ! Interpolate weigths onto the banstructute path.    !
        !====================================================!

        !! Resonant weights
        !   Valence
        call interpolate_kpathweights(iv1, iv2, x, y, z,&
          & rvwgrid(:,iv1:iv2), rvw(:,iv1:iv2))
        !   Conduction
        call interpolate_kpathweights(ic1, ic2, x0, y0, z0,&
          & rcwgrid(:,ic1:ic2), rcw(:,ic1:ic2))

        if(fcoup_) then 
          !! Anti-resonant weights
          !   Valence
          call interpolate_kpathweights(iv1, iv2, x, y, z,&
            & arvwgrid(:,iv1:iv2), arvw(:,iv1:iv2))
          !   Conduction
          call interpolate_kpathweights(ic1, ic2, x0, y0, z0,&
            & arcwgrid(:,ic1:ic2), arcw(:,ic1:ic2))
        end if

        ! Writeout
        call writekpathplot()

        !====================================================!

      ! Exciton loop
      end do

      call clear_excitons()

      deallocate(abs2)
      deallocate(rvwgrid, rcwgrid)
      if(allocated(arvwgrid)) deallocate(arvwgrid)
      if(allocated(arcwgrid)) deallocate(arcwgrid)
      deallocate(rvw, rcw)
      if(allocated(arvw)) deallocate(arvw)
      if(allocated(arcw)) deallocate(arcw)

      deallocate(x,y,z)
      deallocate(x0,y0,z0)
      deallocate(inputdata)

    ! iqmt
    end do

    call clear_bandstructure()

    call barrier(callername=trim(thisname))

  else

    write(unitout,'("Info(",a,"): Rank ",i3," is waiting..")')&
      & trim(thisname), mpiglobal%rank

    call barrier(callername=trim(thisname))

  end if

  contains 

    subroutine setup_xyz(ngridk, ikmap, vkl, x, y, z)
      integer(4), intent(in) :: ngridk(3)
      integer(4), intent(in) :: ikmap(0:,0:,0:)
      real(8), intent(in) :: vkl(:,:)
      real(wp), allocatable, intent(out) :: x(:), y(:), z(:)

      integer(4) :: ikx, iky, ikz, iknr
      integer(4) :: nx, ny, nz

      nx = ngridk(1)
      ny = ngridk(2)
      nz = ngridk(3)

      ! Use supercell -1:2 as input to capture periodicity 
      ! and interpolate points inside 0:1 exclusively
      allocate(x(3*nx))
      allocate(y(3*ny))
      allocate(z(3*nz))

      ! Set 0:1 values
      do ikz=1,nz
        do iky=1,ny
          do ikx=1,nx
            iknr = ikmap(ikx-1,iky-1,ikz-1)
            z(ikz+nz) = vkl(3,iknr)
            y(iky+ny) = vkl(2,iknr)
            x(ikx+nx) = vkl(1,iknr)
          end do
        end do
      end do
      x(1:nx) = x(nx+1:2*nx) - 1.0d0
      y(1:ny) = y(ny+1:2*ny) - 1.0d0
      z(1:nz) = z(nz+1:2*nz) - 1.0d0
      x(2*nx+1:3*nx) = x(nx+1:2*nx) + 1.0d0
      y(2*ny+1:3*ny) = y(ny+1:2*ny) + 1.0d0
      z(2*nz+1:3*nz) = z(nz+1:2*nz) + 1.0d0

    end subroutine setup_xyz

    subroutine genweights(i1, i2, n, nk, imap, ikmap, abs2, weight, ikkpmap)
      integer(4), intent(in) :: i1, i2, n, nk
      integer(4), intent(in) :: imap(1:n), ikmap(1:n)
      integer(4), intent(in) :: ikkpmap(1:nk)
      real(8), intent(in) :: abs2(1:n)
      real(8), intent(inout) :: weight(1:nk,i1:i2)

      integer(4) :: i, ip, a, ik, ikkp

      do i = i1, i2
        do a = 1, n
          ip = imap(a)
          if(i == ip) then
            ! Get reference k index of transition a
            ik = ikmap(a)
            ! Get k'=k_c=k-qmt/2 or k'=k_v=k+qmt/2 index of transition a
            ikkp = ikkpmap(ik)
            ! Add contributon to sum
            !   w_{ikc,c} = \Sum_v |A_{(ik,v),(ik',c)}|^2
            ! or
            !   w_{ikv,v} = \Sum_c |A_{(ik,v),(ik',c)}|^2
            weight(ikkp, i) = weight(ikkp, i) + abs2(a)
          end if
        end do
      end do
    end subroutine genweights

    subroutine writeweights()
      use m_getunit
      use m_genfilname

      integer(4) :: un, iv, ic, iknr
      character(256) :: fname
      character(256) :: tdastring, bsetypestring, scrtypestring

      ! Set output name modifiers
      if(fcoup_) then
        tdastring=''
      else
        if(input%xs%bse%chibarq) then 
          tdastring="-TDA-BAR"
        else
          tdastring="-TDA"
        end if
      end if
      bsetypestring = '-'//trim(input%xs%bse%bsetype)//trim(tdastring)
      scrtypestring = '-'//trim(input%xs%screening%screentype)

      ! Make filename
      call genfilname(dirname=trim(exckpathdir), basename="WEIGHTS",&
        & lambda=lambda, iqmt=iq_,&
        & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
        & nar= .not. input%xs%bse%aresbse, filnam=fname)

      call getunit(un)

      open(unit=un, file=trim(fname), form='formatted', action='write')

      write(un,'("#",1x,"BSE excitonic weights")')
      write(un,'("#")')
      write(un,'("#",1x,"Momentum transfer q_mt:",3f12.7)') vqlmt_
      write(un,'("#",1x,"Size of RR part of the Hamiltonian:",i8)') hamsize_
      write(un,'("#")')
      write(un,'("#",1x,"Considered (partially) occupied states:",i8," to",i8)') ivmin, ivmax
      write(un,'("#",1x,"Considered (partially) unoccupied states:",i8," to",i8)') icmin, icmax
      write(un,'("#")')
      write(un,'("# Eigenvector number:",i6," with energy/eV: ", f12.7)'), lambda, evals_(lambda)*h2ev
      write(un,'("#")')
      write(un,'("# Resonant weights")')

      write(un,'("#  Valence band weights")')

      write(un,'("#",a11,2a12,a8,a8,1x,a23)') "vkv1", "vkv2", "vkv2", "iv", "ik_v", "w_v"
      do iv = ivmin, ivmax
        do iknr = 1,nk_max_
          write(un,'(3f12.7, i8, i8,1x, E23.16)')&
            & vkl_(:,iknr), iv, iknr, rvwgrid(iknr, iv)
        end do
      end do

      write(un,'("#  Conduction band weights")')

      write(un,'("#",a11,2a12,a8,a8,1x,a23)') "vkc1", "vkc2", "vkc2", "ic", "ik_c", "w_c"
      do ic = icmin, icmax
        do iknr = 1, nk_max_
          write(un,'(3f12.7, i8, i8, E23.16)')&
            & vkl0_(:,iknr), ic, iknr, rcwgrid(iknr, ic)
        end do
      end do

      if(fcoup_) then 
        write(un,'("#")')
        write(un,'("# Anti-resonant weights")')
        write(un,'("#  Valence band weights")')
        write(un,'("#",a11,2a12,a8,a8,1x,a23)') "vkv1", "vkv2", "vkv2", "iv", "ik_v", "w_v"
        do iv = ivmin, ivmax
          do iknr = 1, nk_max_
            write(un,'(3f12.7, i8, i8, E23.16)') vkl_(:,iknr), iv, iknr, arvwgrid(iknr, iv)
          end do
        end do
        write(un,'("#  Conduction band weights")')
        write(un,'("#",a11,2a12,a8,a8,1x,a23)') "vkc1", "vkc2", "vkc2", "ic", "ik_c", "w_c"
        do ic = icmin, icmax
          do iknr = 1, nk_max_
            write(un,'(3f12.7, i8, i8, E23.16)') vkl0_(:,iknr), ic, iknr, arcwgrid(iknr, ic)
          end do
        end do
      end if

      close(un)
    end subroutine writeweights

    subroutine interpolate_kpathweights(i1, i2, x, y, z, wsource, winterp)
      integer(4), intent(in) :: i1, i2
      real(8), intent(in) :: wsource(:,:)
      real(wp), intent(in) :: x(:), y(:), z(:)
      real(8), intent(out) :: winterp(:,:)

      integer(4) :: i, iflag, istep, nsteps, ivg(3)
      real(8) :: vklpath_remapped(3)
      real(wp) :: xp,yp,zp,intervalue
      type(bspline_3d) :: spline
      real(8), parameter :: epslat=1.0d-8

      ! No extrapolation needed, since data was expanded to surrounding cells
      logical, parameter :: fexrapol=.false.
      ! Use linear interpolation (spline order 2)
      integer(4) :: ordx
      integer(4) :: ordy
      integer(4) :: ordz
      ! Inderested in the 0'th order derivate
      integer(4), parameter :: derivx=0
      integer(4), parameter :: derivy=0
      integer(4), parameter :: derivz=0

      ordx=intorder
      ordy=intorder
      ordz=intorder

      nsteps = brange_(3)

      do i = 1, i2-i1+1

        call map_inputdata(ngridk_, ikmap_, inputdata, wsource(:,i))

        ! Setup interpolator
        call spline%initialize(x, y, z, inputdata, ordx, ordy, ordz, iflag, fexrapol)
        if(iflag /= 0) then 
          write(unitout,'("Error(m_writebevec):&
            & Spline init returned non zero iflag:", i8)') iflag
          write(unitout,'(" ",a)') get_status_message(iflag)
          call terminate
        end if

        ! Interpolate for all point along the path
        do istep = 1, nsteps

          ! If path k-point is outside the unit cell, map it back
          vklpath_remapped(:) = vklpath_(:, istep, i+i1-1)
          call r3frac(epslat, vklpath_remapped, ivg)

          ! Set point to be interpolated
          xp = real(vklpath_remapped(1), wp)
          yp = real(vklpath_remapped(2), wp)
          zp = real(vklpath_remapped(3), wp)

          call spline%evaluate(xp, yp, zp, derivx, derivy, derivz, intervalue, iflag)
          if(iflag /= 0) then 
            write(unitout,'("Error(m_writebevec):&
              & spline eval returned non zero iflag:", i8)') iflag
            write(unitout,'(" ",a)') get_status_message(iflag)
            call terminate
          end if

          winterp(istep, i) = intervalue

        end do

        call spline%destroy()

      end do

    end subroutine interpolate_kpathweights

    subroutine map_inputdata(ngridk, ikmap, inputdata, weights)
      integer(4), intent(in) :: ngridk(3)
      integer(4), intent(in) :: ikmap(0:,0:,0:)
      real(8), intent(in) :: weights(:)
      real(wp), intent(out) :: inputdata(:,:,:)

      integer(4) :: ikx, iky, ikz, iknr
      integer(4) :: nx, ny, nz
      integer(4) :: ix, iy, iz

      nx = ngridk(1)
      ny = ngridk(2)
      nz = ngridk(3)

      ! Copy to 0:1 part of input data
      do ikz=1, nz
        do iky=1, ny
          do ikx=1, nx
            iknr = ikmap(ikx-1,iky-1,ikz-1)
            inputdata(ikx+nx, iky+ny, ikz+nz) = weights(iknr)
          end do
        end do
      end do

      ! Copy from 0:1 to surrounding cells
      do iz = 1, 3*nz, nz
        do iy = 1, 3*ny, ny
          do ix = 1, 3*nx, nx
            inputdata(ix:ix+nx-1, iy:iy+ny-1, iz:iz+nz-1) =&
              & inputdata(nx+1:2*nx, ny+1:2*ny, nz+1:2*nz)
          end do
        end do
      end do

    end subroutine map_inputdata

    subroutine writekpathplot()
      use m_getunit
      use m_genfilname
      use modxs, only:escale

      integer(4) :: un, ib, istep, nsteps, ib1, ib2
      character(256) :: fname
      character(256) :: tdastring, bsetypestring, scrtypestring

      ! Set output name modifiers
      if(fcoup_) then
        tdastring=''
      else
        if(input%xs%bse%chibarq) then
          tdastring="-TDA-BAR"
        else
          tdastring="-TDA"
        end if
      end if
      bsetypestring = '-'//trim(input%xs%bse%bsetype)//trim(tdastring)
      scrtypestring = '-'//trim(input%xs%screening%screentype)

      call genfilname(dirname=trim(exckpathdir), basename="KPATH",&
        & lambda=lambda, iqmt=iq_,&
        & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
        & nar= .not. input%xs%bse%aresbse, filnam=fname)

      ib1=brange_(1)
      ib2=brange_(2)
      nsteps=brange_(3)

      call getunit(un)

      open(unit=un, file=trim(fname), form='formatted', action='write')
      write(un, '("#",1x,"Bandstructure with excitonic weights")')
      write(un, '("#",1x,"escale:",f12.6)') escale
      write(un, '("#")')
      if(fcoup_) then 
        write(un, '("#",a23,a24,a24,a24)')&
          & "kpathlength","energy","res. exc. weight", "ares. exc. weight"
      else
        write(un, '("#",a23,a24,a24)') "kpathlength","energy","res. exc. weight"
      end if

      do ib = ib1, ib2
        do istep = 1, nsteps

          write(un, '(2E24.16)', advance="no")&
            & kpathlength_(istep, ib), energyval_(istep, ib)*escale

          if(iv1 <= ib .and. ib <= iv2) then
            if(fcoup_) then 
              write(un, '(2E24.16)') rvw(istep, ib), arvw(istep, ib) 
            else
              write(un, '(E24.16)') rvw(istep, ib)
            end if
          else if(ic1 <= ib .and. ib <= ic2) then 
            if(fcoup_) then
              write(un, '(2E24.16)') rcw(istep, ib), arcw(istep, ib) 
            else
              write(un, '(E24.16)') rcw(istep, ib)
            end if
          else
            if(fcoup_) then 
              write(un, '(2E24.16)') 0.0d0, 0.0d0
            else
              write(un, '(E24.16)') 0.0d0
            end if
          end if

        end do
        write(un, *)
      end do

      close(un)

    end subroutine writekpathplot

end subroutine writekpathweights
