module mod_wannier_filehandling
  use mod_wannier_variables

  use mod_APW_LO,                only: apwordmax, apword, nlorb, lorbl, nlotot, lolmax, lolmmax, nlomax
  use mod_eigensystem,           only: idxlo
  use mod_atoms,                 only: natmtot, nspecies, natoms, idxas, spsymb
  use constants,                 only: zzero, zone, twopi
  use mod_muffin_tin,            only: idxlm, lmmaxapw, nrmtmax
  use mod_Gvector,               only: igfft, ngrid
  use mod_Gkvector,              only: ngkmax_ptr
  use mod_eigensystem,           only: nmatmax_ptr
  use mod_spin,                  only: nspnfv
  use mod_eigenvalue_occupancy,  only: nstfv
  use mod_potential_and_density, only: xctype
  use mod_misc,                  only: filext
  use modinput
  use m_getunit
  use m_plotmat
  use modmpi

  implicit none

! methods
  contains

!*********************************!
!          FILE HANDLING          !
!*********************************!
    ! reads transformation matrices from file
    subroutine wffile_readtransform( success)
      logical, intent( out) :: success

      ! local variables
      integer :: ik, ix, iy, iz, un, igroup
      integer :: fst_, lst_, nst_, nwf_, nkpt_, ngroups_
      real(8) :: vkl_( 3), vkl_tmp( 3, wf_kset%nkpt)
      character(32) :: method

      call getunit( un)

      success = .true.
      inquire( file=trim( wf_filename)//"_TRANSFORM"//trim( filext), exist=success)
      if( .not. success) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write(*, '("Error (wffile_readtransform): File '//trim( wf_filename)//"_TRANSFORM"//trim( filext)//' does not exist.")')
        end if
        return
      end if
      open( un, file=trim( wf_filename)//"_TRANSFORM"//trim( filext), action='READ', form='UNFORMATTED', status='OLD')

      ! global parameters
      read( un) fst_, lst_, nst_, nwf_, nkpt_, ngroups_, wf_nprojtot
      if( (fst_ .ne. wf_fst) .or. (lst_ .ne. wf_lst)) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *, '("Warning (wffile_readtransform): different band-ranges in input (",I4,":",I4,") and file (",I4,":",I4,").")') wf_fst, wf_lst, fst_, lst_
          write( *, '(" Use data from file.")')
        end if
        wf_fst = fst_
        wf_lst = lst_
      end if
      wf_nst = wf_lst - wf_fst + 1
      if( nwf_ .ne. wf_nwf) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *, '("Warning (wffile_readtransform): different number of Wannier functions in input (",I4,") and file (",I4,").")') wf_nwf, nwf_
          write( *, '(" Use data from file.")')
        end if
        wf_nwf = nwf_
      end if
      if( nkpt_ .ne. wf_kset%nkpt) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *, '("Error (wffile_readtransform): different number of k-points in input (",I4,") and file (",I4,").")') wf_kset%nkpt, nkpt_
        end if
        call terminate
      end if
      if( ngroups_ .ne. wf_ngroups) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *, '("Warning (wffile_readtransform): different groups of bands in input (",I4,") and file (",I4,").")') wf_ngroups, ngroups_
          write( *, '(" Use data from file.")')
          if( allocated( wf_groups)) deallocate( wf_groups)
          allocate( wf_groups( ngroups_))
        end if
        wf_ngroups = ngroups_
      end if

      ! allocate global arrays
      if( allocated( wf_projst)) deallocate( wf_projst)
      allocate( wf_projst( 7, wf_nprojtot))
      if( allocated( wf_omega)) deallocate( wf_omega)
      allocate( wf_omega( wf_nwf))
      if( allocated( wf_omega_i)) deallocate( wf_omega_i)
      allocate( wf_omega_i( wf_nwf))
      if( allocated( wf_omega_d)) deallocate( wf_omega_d)
      allocate( wf_omega_d( wf_nwf))
      if( allocated( wf_omega_od)) deallocate( wf_omega_od)
      allocate( wf_omega_od( wf_nwf))
      if( allocated( wf_omega)) deallocate( wf_omega)
      allocate( wf_omega( wf_nwf))
      if( allocated( wf_centers)) deallocate( wf_centers)
      allocate( wf_centers( 3, wf_nwf))
      if( allocated( wf_transform)) deallocate( wf_transform)
      allocate( wf_transform( wf_fst:wf_lst, wf_nwf, wf_kset%nkpt))
      !if( allocated( wf_evecphase)) deallocate( wf_evecphase)
      !allocate( wf_evecphase( wf_fst:wf_lst, wf_kset%nkpt))

      !projectors
      do ix = 1, wf_nprojtot
        read( un) wf_projst( :, ix)
      end do

      ! groups
      do igroup = 1, wf_ngroups
        read( un) method
        if( (trim( method) .ne. trim( wf_groups( igroup)%method)) .and. (trim( wf_groups( igroup)%method) .ne. '')) then
          write(*,*)
          write(*,'("Warning (wffile_readtransform): Different method found in file and input for group ",i2,".")') igroup
          write(*,'(t10,"file:  ",a)') trim( method)
          write(*,'(t10,"input: ",a)') trim( wf_groups( igroup)%method)
          write(*,'("Input file setting is used.")')
        end if
        read( un) wf_groups( igroup)%fst, wf_groups( igroup)%lst, wf_groups( igroup)%nst 
        read( un) wf_groups( igroup)%fwf, wf_groups( igroup)%lwf, wf_groups( igroup)%nwf 
        read( un) wf_groups( igroup)%nproj
        read( un) wf_groups( igroup)%win_i
        read( un) wf_groups( igroup)%win_o

        if( allocated( wf_groups( igroup)%projused)) deallocate( wf_groups( igroup)%projused)
        allocate( wf_groups( igroup)%projused( wf_nprojtot))
        read( un) wf_groups( igroup)%projused

        if( allocated( wf_groups( igroup)%win_ii)) deallocate( wf_groups( igroup)%win_ii)
        if( allocated( wf_groups( igroup)%win_io)) deallocate( wf_groups( igroup)%win_io)
        if( allocated( wf_groups( igroup)%win_ni)) deallocate( wf_groups( igroup)%win_ni)
        if( allocated( wf_groups( igroup)%win_no)) deallocate( wf_groups( igroup)%win_no)
        if( wf_groups( igroup)%method .eq. 'disSMV' .or. wf_groups( igroup)%method .eq. 'disFull') then
          read( un) fst_, lst_
          allocate( wf_groups( igroup)%win_ii( fst_, lst_))
          allocate( wf_groups( igroup)%win_io( fst_, lst_))
          read( un) wf_groups( igroup)%win_ii
          read( un) wf_groups( igroup)%win_io

          allocate( wf_groups( igroup)%win_ni( wf_kset%nkpt))
          allocate( wf_groups( igroup)%win_no( wf_kset%nkpt))
          read( un) wf_groups( igroup)%win_ni
          read( un) wf_groups( igroup)%win_no
        else
          allocate( wf_groups( igroup)%win_ii( wf_groups( igroup)%nwf, wf_kset%nkpt), &
                    wf_groups( igroup)%win_io( wf_groups( igroup)%nwf, wf_kset%nkpt))
          allocate( wf_groups( igroup)%win_ni( wf_kset%nkpt), wf_groups( igroup)%win_no( wf_kset%nkpt))
          wf_groups( igroup)%win_ni = wf_groups( igroup)%nwf
          wf_groups( igroup)%win_no = 0
          wf_groups( igroup)%win_io = 0
          do ix = 1, wf_groups( igroup)%nst
            wf_groups( igroup)%win_ii(ix,:) = ix + wf_groups( igroup)%fst - 1
          end do
        end if
      end do

      ! transformation matrices
      do ik = 1, wf_kset%nkpt
        read( un) vkl_
        vkl_tmp( 1, :) = wf_kset%vkl( 1, :) - vkl_( 1)
        vkl_tmp( 2, :) = wf_kset%vkl( 2, :) - vkl_( 2)
        vkl_tmp( 3, :) = wf_kset%vkl( 3, :) - vkl_( 3)
        iz = minloc( norm2( vkl_tmp( :, :), 1), 1)
        if( norm2( vkl_tmp( :, iz)) .gt. input%structure%epslat) then
          if( mpiglobal%rank .eq. 0) then
            write(*,*)
            write( *, '("Error (wffile_readtransform): k-point in file not in k-point-set.")')
            write( *, '(x,3F23.6)') vkl_
          end if
          call terminate
        end if
        do iy = 1, wf_nwf
          do ix = wf_fst, wf_lst
            read( un) wf_transform( ix, iy, iz)
          end do
        end do
      end do

      ! spread
      do ix = 1, wf_nwf
        read( un) wf_omega( ix)
      end do
      do ix = 1, wf_nwf
        read( un) wf_omega_i( ix)
      end do
      do ix = 1, wf_nwf
        read( un) wf_omega_d( ix)
      end do
      do ix = 1, wf_nwf
        read( un) wf_omega_od( ix)
      end do

      ! centers
      do iy = 1, wf_nwf
        do ix = 1, 3
          read( un) wf_centers( ix, iy)
        end do
      end do
      close( un)
      return
    end subroutine wffile_readtransform

    ! writes transformation matrices to file
    subroutine wffile_writetransform
      ! local variables
      integer :: ik, ix, iy, un, igroup
      
      if( mpiglobal%rank .gt. 0) return

      call getunit( un)

      open( un, file=trim( wf_filename)//"_TRANSFORM"//trim( filext), action='WRITE', form='UNFORMATTED')
      ! global parameters
      write( un) wf_fst, wf_lst, wf_nst, wf_nwf, wf_kset%nkpt, wf_ngroups, wf_nprojtot

      !projectors
      do ix = 1, wf_nprojtot
        write( un) wf_projst( :, ix)
      end do

      ! groups
      do igroup = 1, wf_ngroups
        write( un) wf_groups( igroup)%method
        write( un) wf_groups( igroup)%fst, wf_groups( igroup)%lst, wf_groups( igroup)%nst
        write( un) wf_groups( igroup)%fwf, wf_groups( igroup)%lwf, wf_groups( igroup)%nwf
        write( un) wf_groups( igroup)%nproj
        write( un) wf_groups( igroup)%win_i
        write( un) wf_groups( igroup)%win_o

        write( un) wf_groups( igroup)%projused

        if( wf_groups( igroup)%method .eq. 'disSMV' .or. wf_groups( igroup)%method .eq. 'disFull') then
          write( un) shape( wf_groups( igroup)%win_ii)
          write( un) wf_groups( igroup)%win_ii
          write( un) wf_groups( igroup)%win_io

          write( un) wf_groups( igroup)%win_ni
          write( un) wf_groups( igroup)%win_no
        end if
      end do

      ! transformation matrices
      do ik = 1, wf_kset%nkpt
        write( un) wf_kset%vkl( :, ik)
        do iy = 1, wf_nwf
          do ix = wf_fst, wf_lst
            write( un) wf_transform( ix, iy, ik)
          end do
        end do
      end do

      ! spread
      do ix = 1, wf_nwf
        write( un) wf_omega( ix)
      end do
      do ix = 1, wf_nwf
        write( un) wf_omega_i( ix)
      end do
      do ix = 1, wf_nwf
        write( un) wf_omega_d( ix)
      end do
      do ix = 1, wf_nwf
        write( un) wf_omega_od( ix)
      end do

      ! centers
      do iy = 1, wf_nwf
        do ix = 1, 3
          write( un) wf_centers( ix, iy)
        end do
      end do
      !do ik = 1, wf_kset%nkpt
      !  do iy = wf_fst, wf_lst
      !    write( un) wf_evecphase( iy, ik)
      !  end do
      !end do
      close( un)
      return
    end subroutine wffile_writetransform

    subroutine wffile_reademat( success)
      logical, intent( out) :: success

      ! local variables
      integer :: i, j, idxn, ik, ix, iy, iz, un
      integer :: fst_, lst_, nst_, ntot_, nkpt_
      real(8) :: vln(3), vkl_( 3), vkl_tmp( 3, wf_kset%nkpt)
      complex(8) :: ztmp

      call getunit( un)

      success = .true.
      inquire( file=trim( wf_filename)//"_EMAT"//trim( filext), exist=success)
      if( .not. success) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write(*, '("Error (wffile_reademat): File '//trim( wf_filename)//"_EMAT"//trim( filext)//' does not exist.")')
        end if
        return
      end if
      open( un, file=trim( wf_filename)//"_EMAT"//trim( filext), action='READ', form='UNFORMATTED', status='OLD')
      read( un) fst_, lst_, nst_, ntot_, nkpt_
      if( (fst_ .gt. wf_fst) .or. (lst_ .lt. wf_lst)) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *, '("Error (wffile_reademat): bands in input (",I4,":",I4,") out of file band range (",I4,":",I4,").")') wf_fst, wf_lst, fst_, lst_
        end if
        success = .false.
        return
      end if
      if( ntot_ .ne. wf_n_ntot) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *, '("Error (wffile_reademat): different number of BZ-neighbors in input (",I4,") and file (",I4,").")') wf_n_ntot, ntot_
        end if
        success = .false.
        return
      end if
      if( nkpt_ .ne. wf_kset%nkpt) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *, '("Error (wffile_reademat): different number of k-points in input (",I4,") and file (",I4,").")') wf_kset%nkpt, nkpt_
        end if
        success = .false.
        return
      end if
      if( allocated( wf_m0)) deallocate( wf_m0)
      allocate( wf_m0( wf_fst:wf_lst, wf_fst:wf_lst, wf_kset%nkpt, wf_n_ntot))
      !if( allocated( wf_evecphase)) deallocate( wf_evecphase)
      !allocate( wf_evecphase( wf_fst:wf_lst, wf_kset%nkpt))
      do i = 1, wf_n_ntot
        read( un) vln
        ! find index of neighbor
        idxn = 0
        do j = 1, wf_n_ntot
          if( norm2( wf_n_vl( :, j) - vln) .lt. input%structure%epslat) then
              idxn = j
              exit
          end if
        end do
        if( idxn .gt. 0) then
          do ik = 1, wf_kset%nkpt
            read( un) vkl_
            vkl_tmp( 1, :) = wf_kset%vkl( 1, :) - vkl_( 1)
            vkl_tmp( 2, :) = wf_kset%vkl( 2, :) - vkl_( 2)
            vkl_tmp( 3, :) = wf_kset%vkl( 3, :) - vkl_( 3)
            iz = minloc( norm2( vkl_tmp( :, :), 1), 1)
            if( norm2( wf_kset%vkl( :, iz) - vkl_(:)) .gt. input%structure%epslat) then
              if( mpiglobal%rank .eq. 0) then
                write(*,*)
                write( *, '("Error (wffile_reademat): k-point in file not in k-point-set.")')
                write( *, '(x,3F23.6)') vkl_
              end if
              success = .false.
              return
            end if
            do iy = fst_, lst_
              do ix = fst_, lst_
                read( un) ztmp
                if( (ix .ge. wf_fst) .and. (ix .le. wf_lst) .and. (iy .ge. wf_fst) .and. (iy .le. wf_lst)) wf_m0( ix, iy, ik, idxn) = ztmp
              end do
            end do
          end do
        else
          if( mpiglobal%rank .eq. 0) then
            write(*,*)
            write( *, '("Error (wffile_reademat): neighboring vector in file not consistent with input.")')
            write( *, '(x,3F23.6)') vln
          end if
          success = .false.
          return
        end if
      end do
      !do ik = 1, wf_kset%nkpt
      !  do iy = wf_fst, wf_lst
      !    read( un) wf_evecphase( iy, ik)
      !  end do
      !end do
      close( un)
      return
    end subroutine wffile_reademat

    ! writes transformation matrices to file
    subroutine wffile_writeemat
      ! local variables
      integer :: ik, idxn, ix, iy, un

      if( mpiglobal%rank .gt. 0) return

      call getunit( un)

      open( un, file=trim( wf_filename)//"_EMAT"//trim( filext), action='WRITE', form='UNFORMATTED')
      write( un) wf_fst, wf_lst, wf_nst, wf_n_ntot, wf_kset%nkpt
      do idxn = 1, wf_n_ntot
        write( un) wf_n_vl( :, idxn)
        do ik = 1, wf_kset%nkpt
          write( un) wf_kset%vkl( :, ik)
          do iy = wf_fst, wf_lst
            do ix = wf_fst, wf_lst
              write( un) wf_m0( ix, iy, ik, idxn)
            end do
          end do
        end do
      end do
      !do ik = 1, wf_kset%nkpt
      !  do iy = wf_fst, wf_lst
      !    write( un) wf_evecphase( iy, ik)
      !  end do
      !end do
      close( un)
      !write( *, '(a,a)') ' Plane-wave matrix-elements written to file ', trim( wf_filename)//"_EMAT"//trim( filext)
      return
    end subroutine wffile_writeemat

    subroutine wffile_delemat
      integer :: un
      logical :: exist

      inquire( file=trim( wf_filename)//"_EMAT"//trim( filext), exist=exist)
      if( exist) then
        call getunit( un)
        open( un, file=trim( wf_filename)//"_EMAT"//trim( filext))
        close( un, status='delete')
      end if

      return
    end subroutine wffile_delemat

    subroutine wffile_writeinfo_lo
      integer :: i, j

      if( mpiglobal%rank .gt. 0) return

      call printbox( wf_info, '*', "Local-orbitals for projection")
      write( wf_info, *)

      write( wf_info, '(7x,"#",5x,"species",1x,"atom",9x,"cell",5x,"  n /  l /  m ",3x,"dord",6x,"groups")')
      write( wf_info, '(80("-"))')
      do i = 1, wf_nprojtot
        write( wf_info, '(4x,i4,10x,a2,2x,i3,4x,3(i3),5x,2(i3," /"),i3,4x,i4,4x)', advance='no') &
            i, &
            spsymb( wf_projst( 1, i)), &
            wf_projst( 2, i), &
            wf_projst( 8:10, i), &
            wf_projst( 4:7, i)
        do j = 1, wf_ngroups
          if( wf_groups( j)%projused( i) .eq. 1) then
            write( wf_info, '(i2,x)', advance='no') j
          else
            write( wf_info, '(3x)', advance='no')
          end if
        end do
        write( wf_info, *)
      end do

      write( wf_info, '(80("-"))')
      write( wf_info, '(36x,"local-orbitals used in total:",4x,I4)') wf_nprojtot
      write( wf_info, *)
      call flushifc( wf_info)
    end subroutine wffile_writeinfo_lo
    
    subroutine wffile_writeinfo_geometry
      integer :: i
      real(8) :: d, v(3,1), m(3,3)

      if( mpiglobal%rank .gt. 0) return

      call printbox( wf_info, '*', "Brillouin zone neighbors for k-gradient")
      write( wf_info, *)

      write( wf_info, '(" vectors to neighboring k-points (cartesian)")')
      write( wf_info, *)

      m = 0.d0
      d = 0.d0
      do i = 1, wf_n_ntot
        if( abs( wf_n_dist( i) - d) .gt. input%structure%epslat) then
          d = wf_n_dist( i)
          write( wf_info, '(6x,"length: ",f16.10,14x,"weight: ",f16.10)') wf_n_dist( i), wf_n_wgt( i)
        end if
        write( wf_info, '(8x,I2,4x,3(F22.10))') i, wf_n_vc( :, i)
        v(:,1) = wf_n_vc( :, i)
        m = m + 2.0d0*wf_n_wgt( i)*matmul( v, transpose( v))
      end do

      write( wf_info, *)
      write( wf_info, '(" vectors to neighboring k-points (lattice)")')
      write( wf_info, *)

      d = 0.d0
      do i = 1, wf_n_ntot
        if( abs( wf_n_dist( i) - d) .gt. input%structure%epslat) then
          d = wf_n_dist( i)
          write( wf_info, '(6x,"length: ",f16.10,14x,"weight: ",f16.10)') wf_n_dist( i), wf_n_wgt( i)
        end if
        write( wf_info, '(8x,I2,4x,3(F22.10))') i, wf_n_vl( :, i)
      end do

      write( wf_info, *)
      write( wf_info, '(" consistency check (fullfilled if unity, consider grid-dimension)")')
      write( wf_info, *)
      write( wf_info, '(14x,3(F22.10))') m(1,:)
      write( wf_info, '(14x,3(F22.10))') m(2,:)
      write( wf_info, '(14x,3(F22.10))') m(3,:)


      write( wf_info, *)
      call flushifc( wf_info)
    end subroutine wffile_writeinfo_geometry
    
    subroutine wffile_writeinfo_overall
      if( mpiglobal%rank .gt. 0) return

      call printbox( wf_info, '*', "Overall set-up")
      write( wf_info, *)
      write( wf_info, '(" lowest band:",T30,13x,I3)') wf_fst
      write( wf_info, '(" highest band:",T30,13x,I3)') wf_lst
      write( wf_info, '(" #bands involved:",T30,13x,I3)') wf_nst
      write( wf_info, '(" #Wannier functions:",T30,13x,I3)') wf_nwf
      write( wf_info, '(" #groups:",T30,13x,I3)') wf_ngroups
      write( wf_info, '(" #k-points:",T30,11x,I5)') wf_kset%nkpt
      call flushifc( wf_info)
    end subroutine wffile_writeinfo_overall
    
    subroutine wffile_writeinfo_task
      character(16) :: string
     
      if( mpiglobal%rank .gt. 0) return

      write( string, '("Group: ",I2)') wf_group
      call printbox( wf_info, '-', string)
      string = wf_groups( wf_group)%method(1:16)
      write( wf_info, *)
      write( wf_info, '(" lowest band:",T30,13x,I3)') wf_groups( wf_group)%fst
      write( wf_info, '(" highest band:",T30,13x,I3)') wf_groups( wf_group)%lst
      write( wf_info, '(" #bands involved:",T30,13x,I3)') wf_groups( wf_group)%nst
      write( wf_info, '(" #Wannier functions:",T30,13x,I3)') wf_groups( wf_group)%nwf
      if( wf_groups( wf_group)%neighcells) then
        write( wf_info, '(" #projection functions (nc):",T30,11x,I5)') wf_groups( wf_group)%nproj
      else
        write( wf_info, '(" #projection functions:",T30,11x,I5)') wf_groups( wf_group)%nproj
      end if
      write( wf_info, '(" method:",T30,A16)') adjustr( string)
      if( wf_groups( wf_group)%method .eq. 'disSMV' .or. wf_groups( wf_group)%method .eq. 'disFull') then
        write( wf_info, '(" outer window:",T30,f7.4,2x,f7.4)') wf_groups( wf_group)%win_o
        write( wf_info, '(" inner window:",T30,f7.4,2x,f7.4)') wf_groups( wf_group)%win_i
        write( wf_info, '(" min/max bands inner window: ",T30,7x,I4,"/",I4)') minval( wf_groups( wf_group)%win_ni), maxval( wf_groups( wf_group)%win_ni)
        write( wf_info, '(" min/max bands outer window: ",T30,7x,I4,"/",I4)') minval( wf_groups( wf_group)%win_ni+wf_groups( wf_group)%win_no), maxval( wf_groups( wf_group)%win_ni+wf_groups( wf_group)%win_no)
      end if
      write( wf_info, *)
      call flushifc( wf_info)
    end subroutine wffile_writeinfo_task

    subroutine wffile_writeinfo_finish
      real(8) :: t
      character(64) :: string
      
      if( mpiglobal%rank .gt. 0) return

      call timesec( t)
      write( string, '("total duration (seconds):",T30,F16.1)') t-wf_t0
      call printbox( wf_info, '-', string)
      call flushifc( wf_info)
      call wffile_writeinfo_results
    end subroutine wffile_writeinfo_finish

    subroutine wffile_writeinfo_results
      use mod_lattice, only: ainv
      integer :: i, igroup
      real(8) :: vl(3)

      if( mpiglobal%rank .gt. 0) return

      call printbox( wf_info, '*', "Wannier functions")
      write( wf_info, *)
      write( wf_info, '(3x,"#",4x,"localization center (lattice)",8x,"Omega",3x,"Omega_I",3x,"Omega_D",2x,"Omega_OD")')
      write( wf_info, '(80("="))')

      do igroup = 1, wf_ngroups
        do i = wf_groups( igroup)%fwf, wf_groups( igroup)%lwf
          call r3mv( ainv, wf_centers( :, i), vl)
          write( wf_info, '(I4,3x,3F10.4,3x,4F10.4)') i, vl, wf_omega( i), wf_omega_i( i), wf_omega_d( i), wf_omega_od( i)
        end do
        if( wf_ngroups .gt. 1) then
          if( igroup .le. wf_ngroups) write( wf_info, '(80("-"))')
          write( wf_info, '(34x,"total:",4F10.4)') sum( wf_omega( wf_groups( igroup)%fwf:wf_groups( igroup)%lwf)), &
                                                   sum( wf_omega_i( wf_groups( igroup)%fwf:wf_groups( igroup)%lwf)), &
                                                   sum( wf_omega_d( wf_groups( igroup)%fwf:wf_groups( igroup)%lwf)), &
                                                   sum( wf_omega_od( wf_groups( igroup)%fwf:wf_groups( igroup)%lwf))
        end if
        if( igroup .lt. wf_ngroups) write( wf_info, '(80("-"))')
      end do

      write( wf_info, '(80("="))')
      write( wf_info, '(34x,"total:",4F10.4)') sum( wf_omega), sum( wf_omega_i), sum( wf_omega_d), sum( wf_omega_od)
      write( wf_info, '(32x,"average:",4F10.4)') sum( wf_omega)/wf_nwf, sum( wf_omega_i)/wf_nwf, sum( wf_omega_d)/wf_nwf, sum( wf_omega_od)/wf_nwf

      write( wf_info, *)
      call flushifc( wf_info)
      close( wf_info)
    end subroutine wffile_writeinfo_results

end module mod_wannier_filehandling
