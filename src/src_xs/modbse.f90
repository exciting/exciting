module modbse
  
  use modinput, only: input
  use modxs, only: bcbs

  implicit none

  ! Band combinations in relative and absolute indices
  Type(bcbs) :: bcou, bcouabs
  ! Reference states for occupied and unoccupied state indices
  integer(4) :: ioref, iuref 

  contains

    subroutine setranges_modxs(iq)
      use modxs, only: istocc0, istocc, istunocc0, istunocc,&
                     & isto0, isto, istu0, istu
      implicit none
      integer(4), intent(in) :: iq
      
      ! Find occupation limits for k and k+q and set variables in modxs
      call findocclims(iq, istocc0, istocc, istunocc0,&
        & istunocc, isto0, isto, istu0, istu)

      if(input%xs%dbglev > 2) then
        write(*,*) "modbse:setranges_modxs"
        write(*,*) "iq", iq
        write(*,*) "istocc0, istunocc0"
        write(*,*) istocc0, istunocc0
        write(*,*) "isto0"
        write(*,*) isto0
        write(*,*) "istu0"
        write(*,*) istu0
        write(*,*) "istocc, istunocc"
        write(*,*) istocc, istunocc
        write(*,*) "isto"
        write(*,*) isto
        write(*,*) "istu"
        write(*,*) istu
      end if

    end subroutine setranges_modxs

    subroutine setbcbs_bse
      use modxs, only: istunocc0, unitout
      use modinput, only: input
      implicit none

      !!<-- READ INPUT AND SET "CONVECTIONAL" INDICES
      ! Lowest occupied state (absolute index)
      bcou%il1 = input%xs%bse%nstlbse(1)
      ! Highest occupied state (absolute index)
      bcou%iu1 = input%xs%bse%nstlbse(2)
      ! Lowest unoccupied state (counted from first unoccupied state)
      bcou%il2 = input%xs%bse%nstlbse(3)
      ! Highest unoccupied state (counted from first unoccupied state)
      bcou%iu2 = input%xs%bse%nstlbse(4)
      ! Number of occupied states
      bcou%n1 = bcou%iu1-bcou%il1+1
      ! Number of unoccupied states
      bcou%n2 = bcou%iu2-bcou%il2+1
      !!-->

      !!<-- Set reference states
      ioref = 1
      iuref = istunocc0
      !!-->

      !!<-- SET ABSOLUTE INDICES
      bcouabs%il1 = bcou%il1
      bcouabs%iu1 = bcou%iu1
      bcouabs%il2 = bcou%il2 + iuref - 1 
      bcouabs%iu2 = bcou%iu2 + iuref - 1 
      bcouabs%n1 = bcou%n1
      bcouabs%n2 = bcou%n2
      !!-->

      if(input%xs%dbglev > 2) then
        write(*,*) "modbse:setbcbs_bse"
        write(*,*) "bcou", bcou%n1, bcou%il1, bcou%iu1, bcou%n2, bcou%il2, bcou%iu2
        write(*,*) "bcouabs", bcouabs%n1, bcouabs%il1,&
         & bcouabs%iu1, bcouabs%n2, bcouabs%il2, bcouabs%iu2
      end if

      ! Write out state ranges to INFOXS.OUT
      write(unitout,*)
      write(unitout, '("Info(bse): Information on number of states:")')
      write(unitout, '("  Ranges of states used in construction of BSE matrix:")')
      write(unitout, '("    Range of occupied states and number   :", 2i6, 3x, i6)')&
        & bcouabs%il1, bcouabs%iu1, bcouabs%n1
      write(unitout, '("    Range of unoccupied states and number :", 2i6, 3x, i6)')&
        & bcouabs%il2, bcouabs%iu2, bcouabs%n2

    end subroutine setbcbs_bse

    !!<-- COMBINED INDEX MAPPING
    ! hamidx calculates a combined index that labels the states as
    ! hamidx = {o1u1k1, o1u2k1, ..., o1uMk1,
    !      o2uMk1, ..., oLuMk1, o1u1k2, ..., oLuMkN} -> {1,...,L*M*N}
    ! (if i1=o n1=L i2=u n2=M)
    integer(4) function hamidx(i1, i2, ik, n1, n2)
      implicit none
      integer(4), intent(in) :: i1, i2, ik, n1, n2
      hamidx = i2 + n2 * (i1-1) + n1 * n2 * (ik-1)
    end function hamidx
    ! hamidx_back calculates the individual indices given a combined one.
    subroutine hamidx_back(s, i1, i2, ik, n1, n2)
      implicit none
      integer(4), intent(in) :: s, n1, n2
      integer(4), intent(out) :: i1, i2, ik
      integer(4) :: n12, tmp
      n12 = n1*n2
      ik = (s-1)/n12 + 1
      tmp = s - (ik-1)*n12
      i1 = (tmp-1)/n2 + 1
      i2 = tmp - (i1-1)*n2
    end subroutine hamidx_back
    !!-->

    subroutine checkoccupancies(iq, nk, bc, hamsize, amap, occfactor)
      use modinput, only: input
      use modxs, only: unitout, bcbs
      use m_getunit

      implicit none 

      integer(4), intent(in) :: iq, nk 
      type(bcbs), intent(in) :: bc
      integer(4), intent(out) :: hamsize
      integer(4), allocatable, intent(out) :: amap(:,:)
      real(8), allocatable, intent(out) :: occfactor(:)

      integer(4) :: ik, i, j, a1, a2, hamsize_max, un
      real(8) :: kdocc, cutoff
      real(8), allocatable :: occfactor_t(:), docc(:,:)
      logical, allocatable :: zeroflag(:), negativeflag(:)

      ! Maximal size of Hamiltonian
      ! Size of BSE-Hamiltonian = #o * #u * #k
      hamsize_max = bc%n1 * bc%n2 * nk

      cutoff = input%groundstate%epsocc
      
      if(.false.) then
        write(*,*) "modbse:checkoccupancies"
        write(*,*) "Maximal Hamilton size:", hamsize_max
        write(*,*) "Cutoff=", cutoff
      end if

      ! Initialize flags for zero occupation difference
      ! and negative occupation difference.
      allocate(zeroflag(hamsize_max))
      allocate(negativeflag(hamsize_max))
      zeroflag = .false.
      negativeflag = .false.
      ! Allocate temporary occupation factor array
      allocate(occfactor_t(hamsize_max))
      ! Allocate helper array
      allocate(docc(bc%n1,bc%n2))
      ! Check occupancies
      do ik = 1, nk
        ! Get occupation numbers 
        !   Calculates docc(l,m)= f_{o_l ki+q} - f_{u_m ki}
        !   Remark: Currently the iq argument does nothing at all
        call getdocc(iq, ik, ik, bc%il1, bc%iu1, bc%il2, bc%iu2, docc)
        


        ! Band indices
        do i = bc%il1, bc%iu1 ! Occupied
          do j = bc%il2, bc%iu2 ! Unoccupied
            ! Map to combined index
            a1 = hamidx(i-bc%il1+1, j-bc%il2+1, ik, bc%n1, bc%n2)

            ! Get occupation difference for current combination
            ! kdocc = f_{o_{a},k_{a}} - f_{u_{a},k_{a}}
!! This assumes KS calculations with max. occupancies of 2
            kdocc = docc(i-bc%il1+1, j-bc%il2+1)/2.0d0
            ! Set zero occupation difference flag
            if(abs(kdocc) .lt. cutoff) then
              zeroflag(a1) = .true.
            end if
            ! Set occupation difference sign flag
            if(kdocc .le. 0.0d0) then
              negativeflag(a1) = .true.
            end if
            ! Occupation factor
            occfactor_t(a1) = sqrt(abs(kdocc))
          end do
        end do
      end do
      deallocate(docc)

      ! Print all occupancy factors with flags
      call printoccupancies(zeroflag, negativeflag, occfactor_t)

      ! Adjust selection for only non zero occupancy differences
      ! and only positive occupancy differences.
      call getunit(un)
      open(un, file='BCBS_SKIPPED.OUT', action='write', status='replace')
      write(un,*) "# Skipped band combinations"
      write(un,*) "# s, io, iu, ik, ochkz, ockn"
      hamsize = 0
      do a1 = 1, hamsize_max
        if(.not. zeroflag(a1) .and. .not. negativeflag(a1)) then
          hamsize = hamsize + 1
        else
          call hamidx_back(a1, i, j, ik, bc%n1, bc%n2)
          write(un,'(I8,1x,I3,1x,I3,1x,I3,1x,L,1x,L)')&
            & a1, i, j, ik, zeroflag(a1), negativeflag(a1)
        end if
      end do
      close(un)
      if (hamsize /= hamsize_max) then
        write(unitout,*) "(bse.f90) [WARNING] Zero occupancy differences and/or&
          & negative fo-fu occured. Reduced size of BSE Hamiltonian from, ",&
          & hamsize_max, "to ", hamsize, " ."
      end if 

      ! Write out skipped band combinations
      ! and make map of combined indices to use.
      allocate(amap(hamsize, 3))
      allocate(occfactor(hamsize))
      a2=0
      do a1 = 1, hamsize_max
        if(.not. zeroflag(a1) .and. .not. negativeflag(a1)) then
          a2=a2+1
          occfactor(a2) = occfactor_t(a1)
          call hamidx_back(a1, i, j, ik, bc%n1, bc%n2)
          ! Corresponds to smap(i,:) = [io, iu, ik]
          amap(a2,:) = [i, j, ik]
        end if
      end do
      deallocate(occfactor_t, zeroflag, negativeflag)
      
      ! Write out index mapping
      call getunit(un)
      open(un, file='BSE_SINDEX.OUT', action='write', status='replace')
      write(un,*) "# s ik io iu"
      do a1 = 1, hamsize
        write(un,'(I3,1x,I3,1x,I3,1x,I3)')&
          & a1, amap(a1,3), amap(a1,1), amap(a1,2)
      end do
      close(un)

    end subroutine checkoccupancies

    subroutine printoccupancies(zeroflag, negativeflag, occfactor)
      use m_getunit
      implicit none 
      logical, intent(in) :: zeroflag(:), negativeflag(:)
      real(8), intent(in) :: occfactor(:)

      integer(4) :: i, un

      call getunit(un)
      open(un, file='BSE_OCCUPANCIES_ALL.OUT', action='write', status='replace')
      write(un,'(a)') "zeroflag negativeflag occfac"
      do i = 1, size(occfactor)
        write(un, '(L,1x,L,1x,E23.16)')&
          & zeroflag(i), negativeflag(i), occfactor(i)
      end do
      close(un)

    end subroutine printoccupancies

end module modbse
