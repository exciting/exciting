!
!
!
subroutine getpmatkgw(ik)

    use modinput
    use modmain
    use modgw
    use mod_dielectric_function
    use mod_hdf5
    implicit none
    
    ! input
    integer, intent(in) :: ik
    ! local
    integer :: ikp, isym, lspl, iv(3)
    integer :: ie1, ie2, icg, is, ia, ias, ic
    real(8) :: s(3,3), v1(3), v2(3), v3(3), t1
    logical :: lfound
    integer :: recl
    
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
    if (ik .ne. kset%ikp2ik(ikp)) then
    
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
        write(*,'("Error(getpmatkgw): No symmetry operation connecting &
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

    return
end subroutine
