!
!
!
subroutine getpmatkgw(ik)

    use modinput, only: input
    !use modmain
    use modgw, only: kset, kqset
    use mod_bands, only: numin, nomax, nstdf
    use mod_core_states, only: ncg
    use mod_dielectric_function, only: pmatcv, pmatvv, fname_pmatcv, fname_pmatvv
    use mod_large_io, only: inquire_large, open_direct_unformatted_large
    use mod_symmetry, only: lsplsymc, nsymcrys, symlat, symlatc
    use precision, only: i32, long_int, dp

    implicit none
    
    ! input
    integer(i32), intent(in) :: ik
    ! local
    integer(i32) :: ikp, isym, lspl, iv(3)
    integer(i32) :: ie1, ie2, icg, ias
    real(dp) :: s(3,3), v1(3), v2(3), v3(3), t1
    logical :: lfound
    logical(i32) :: calculate_core
    integer(long_int) :: recl
    integer(i32) :: fid_pmatvv
    integer(i32) :: fid_pmatcv
    

    !---------
    ! val-val
    !---------
    if (allocated(pmatvv)) deallocate(pmatvv)
    allocate(pmatvv(nomax,numin:nstdf,3))
    call inquire_large( recl, pmatvv )
    call open_direct_unformatted_large( fid_pmatvv, fname_pmatvv, "read", recl, "old" )
    !----------
    ! core-val
    !----------
    calculate_core = (input%gw%coreflag=='all')
    if ( calculate_core ) then
      if (allocated(pmatcv)) deallocate(pmatcv)
      allocate(pmatcv(ncg,numin:nstdf,3))
      call inquire_large( recl, pmatcv )
      call open_direct_unformatted_large( fid_pmatcv, fname_pmatcv, "read", recl, "old" )
    end if
    
    ikp = kset%ik2ikp(ik)
    
    !=========================
    ! Read the data from file
    !=========================
    
    read(fid_pmatvv,rec=ikp) pmatvv
    if (input%gw%coreflag=='all') then
      read(fid_pmatcv,rec=ikp) pmatcv
    end if
    
    !-----------------------------------
    ! check if the k-point is reducible
    !-----------------------------------
    if (ik /= kset%ikp2ik(ikp)) then
    
      !=====================================================
      ! Find symmetry operation connecting ik and ikp points
      !=====================================================
      lfound = .false.
      do isym = 1, nsymcrys
        lspl = lsplsymc(isym)
        s(:,:) = dble(symlat(:,:,lspl))
        call r3mtv(s,kqset%vkl(:,ik),v1)
        call r3frac(input%structure%epslat,v1,iv)
        v2(:) = kset%vkl(:,ikp)
        call r3frac(input%structure%epslat,v2,iv)
        t1 = dabs(v1(1)-v2(1)) + &
        &    dabs(v1(2)-v2(2)) + &
        &    dabs(v1(3)-v2(3))
        if (t1 < input%structure%epslat) then
          lfound = .true.
          exit ! isym loop
        end if
      end do ! isym
      if (.not.lfound) then
        write(*,*)
        write(*,'("Error(getpmatkgw): No symmetry operation that connects &
        &kset and kqset points is found!")')
        write(*,*)
        stop
      end if
      
      !-------------------------------------------------------------------
      ! rotate the matrix element from the reduced to non-reduced k-point
      ! (note that the inverse operation is used)
      !-------------------------------------------------------------------
      ! val-val
      do ie2 = numin, nstdf
        do ie1 = 1, nomax
          v1(:) = dble(pmatvv(ie1,ie2,:))
          call r3mv(symlatc(:,:,lspl),v1,v2)
          v1(:) = aimag(pmatvv(ie1,ie2,:))
          call r3mv(symlatc(:,:,lspl),v1,v3)
          pmatvv(ie1,ie2,:) = cmplx(v2(:),v3(:),8)
        end do ! ie1
      end do ! ie2
      
      if (input%gw%coreflag=='all') then
        ! cor-val
        do ie2 = numin, nstdf
          do icg = 1, ncg
            v1(:) = dble(pmatcv(icg,ie2,:))
            call r3mv(symlatc(:,:,lspl),v1,v2)
            v1(:) = aimag(pmatcv(icg,ie2,:))
            call r3mv(symlatc(:,:,lspl),v1,v3)
            pmatcv(icg,ie2,:) = cmplx(v2(:),v3(:),8)
          end do ! icg
        end do ! ie2
      end if ! core
      
    end if

    close(fid_pmatvv)
    if ( calculate_core ) close(fid_pmatcv)
end subroutine
