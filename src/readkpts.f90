Subroutine readkpts
      Use modmain
      use mod_kpointset
      Implicit None

      integer :: nk
      real(8) :: vk( 3, nkpt), v1(3), v2(3)
      logical :: xfastest
      type( k_set) :: ksettmp

      Integer :: ik, iq, ix, iy, iz, d2
      real(8) :: d1

      Open (50, File='KPOINTS'//trim(filext), Action='read', Form='formatted')
      read( 50, '(I6, " : nkpt; k-point, vkl, wkpt, nmat below")') nk
      if( nk .ne. nkpt) then
        write( *, '("ERROR (readkpts): Different number of kpoints.")')
        stop
      end if
      Do iq = 1, nk
         read (50, *) ik, vk(:, ik), d1, d2
         !write(*,'(I6,4G18.10,I8)') ik, vk( :, ik), d1, d2
      End Do
      Close (50)

      v1 = vk( :, 2) - input%groundstate%vkloff/dble( input%groundstate%ngridk)
      ik = 0
      lx: do ix = 0, input%groundstate%ngridk( 1) -1
        ly: do iy = 0, input%groundstate%ngridk( 1) -1
          lz: do iz = 0, input%groundstate%ngridk( 1) -1
            ik = ik + 1
            v2 = dble((/ix, iy, iz/))/dble( input%groundstate%ngridk)
            if( ik .eq. 2) then
              if( norm2( v1 - v2) .lt. input%structure%epslat) then
                xfastest = .false.
              else
                xfastest = .true.
              end if
              exit lx
            end if
          end do lz
        end do ly
      end do lx

      if( xfastest) then
        call init1
      else
        call generate_k_vectors( ksettmp, bvec, input%groundstate%ngridk, input%groundstate%vkloff, input%groundstate%reducek)
        vkl = ksettmp%vkl
        vkc = ksettmp%vkc
        nkpt = ksettmp%nkpt
        call generate_k_vectors( ksettmp, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .false.)
        vklnr = ksettmp%vkl
        vkcnr = ksettmp%vkc
        nkptnr = ksettmp%nkpt
      end if
      Return
End Subroutine
!EOC
