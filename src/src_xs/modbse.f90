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
      write(unitout, '("Info(setbcbs): Information on number of states:")')
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
    ! Combinded o-u index at one k
    integer(4) function subhamidx(i1, i2, n2)
      implicit none
      integer(4), intent(in) :: i1, i2, n2
      subhamidx = i2 + n2 * (i1-1)
    end function subhamidx
    subroutine subhamidx_back(s, i1, i2, n2)
      implicit none
      integer(4), intent(in) :: s, n2
      integer(4), intent(out) :: i1, i2
      i1 = (s-1)/n2 + 1
      i2 = s - (i1-1)*n2
    end subroutine subhamidx_back
    !!-->

    subroutine checkoccupancies(iq, nk, bcabs, hamsize, amap, kouflag, kkpblocksize, occfactor)
      use modinput, only: input
      use modxs, only: unitout, bcbs
      use m_getunit

      implicit none 

      integer(4), intent(in) :: iq, nk 
      type(bcbs), intent(in) :: bcabs
      integer(4), intent(out) :: hamsize
      integer(4), allocatable, intent(out) :: amap(:,:), kkpblocksize(:)
      logical, allocatable, intent(out) :: kouflag(:,:)
      real(8), allocatable, intent(out) :: occfactor(:)

      integer(4) :: ik, i, j, a1, a2, hamsize_max, un, nou, nkkp, check
      real(8) :: kdocc, cutoff
      real(8), allocatable :: occfactor_t(:), docc(:,:)
      logical, allocatable :: zeroflag(:), negativeflag(:)
      character(256) :: frmt

      ! Maximal size of Hamiltonian
      ! Size of BSE-Hamiltonian = #o * #u * #k
      nkkp = nk*(nk+1)/2
      nou = bcabs%n1 * bcabs%n2
      hamsize_max =  nou * nk

      cutoff = input%groundstate%epsocc
      
      ! Initialize flags for zero occupation difference
      ! and negative occupation difference.
      allocate(zeroflag(hamsize_max))
      allocate(negativeflag(hamsize_max))
      zeroflag = .false.
      negativeflag = .false.
      ! Allocate temporary occupation factor array
      allocate(occfactor_t(hamsize_max))
      ! Allocate helper array
      allocate(docc(bcabs%n1,bcabs%n2))

      allocate(kouflag(nou, nk))
      kouflag = .true.

      ! (TEST) write out all possible s indices
      call getunit(un)
      open(un, file='BSE_SINDEX_ALL.OUT', action='write', status='replace')
      write(un,'("#",1x,a3,1x,a5,1x,a5,1x,a5)') "s", "ik", "io" ,"iua"

      ! Check occupancies
      do ik = 1, nk
        ! Get occupation numbers 
        !   Calculates docc(l,m)= f_{o_l ki+q} - f_{u_m ki}
        !   Remark: Currently the iq argument does nothing at all
        call getdocc(iq, ik, ik, bcabs%il1, bcabs%iu1, bcabs%il2, bcabs%iu2, docc)
        
        ! Band indices
        do i = bcabs%il1, bcabs%iu1 ! Occupied
          do j = bcabs%il2, bcabs%iu2 ! Unoccupied
            ! Map to combined index
            a1 = hamidx(i-bcabs%il1+1, j-bcabs%il2+1, ik, bcabs%n1, bcabs%n2)

            ! (TEST)
            write(un,'(I5,1x,I5,1x,I5,1x,I5)')&
                & a1, ik, i, j

            ! Get occupation difference for current combination
            ! kdocc = f_{o_{a},k_{a}} - f_{u_{a},k_{a}}
!! This assumes KS calculations with max. occupancies of 2
            kdocc = docc(i-bcabs%il1+1, j-bcabs%il2+1)/2.0d0
            ! Set zero occupation difference flag
            if(abs(kdocc) .lt. cutoff) then
              zeroflag(a1) = .true.
              kouflag(subhamidx(i-bcabs%il1+1,j-bcabs%il2+1, bcabs%n2), ik) = .false.
            end if
            ! Set occupation difference sign flag
            if(kdocc .le. 0.0d0) then
              negativeflag(a1) = .true.
              kouflag(subhamidx(i-bcabs%il1+1,j-bcabs%il2+1, bcabs%n2), ik) = .false.
            end if
            ! Occupation factor
            occfactor_t(a1) = sqrt(abs(kdocc))
          end do
        end do
      end do
      deallocate(docc)

      ! (TEST)
      close(un)

      allocate(kkpblocksize(nk))
      do ik = 1, nk
        kkpblocksize(ik) = count(kouflag(:,ik))
      end do

      ! Print kouflag and kkpblocksize
      call printkouflag(nou, nk, kouflag, kkpblocksize)

      ! Print all occupancy factors with flags
      call printoccupancies(zeroflag, negativeflag, occfactor_t)

      ! Adjust selection for only non zero occupancy differences
      ! and only positive occupancy differences.
      ! Write out skipped band combinations
      ! and make map of combined indices to use.
      call getunit(un)
      open(un, file='BSE_SKIPPED_BCBS.OUT', action='write', status='replace')
      write(un,'("#",1x,a)') "Skipped band combinations"
      frmt='("#",a7,1x,a4,1x,a4,1x,a4,1x,a4,1x,a4,1x,a1,1x,a1)'
      write(un,frmt) "s", "ik", "io", "iu", "ioa", "iua", "z", "n"

      hamsize = 0
      do a1 = 1, hamsize_max
        if(.not. zeroflag(a1) .and. .not. negativeflag(a1)) then
          hamsize = hamsize + 1
        else
          call hamidx_back(a1, i, j, ik, bcabs%n1, bcabs%n2)
          write(un,'(I8,1x,I4,1x,I4,1x,I4,1x,I4,1x,I4,1x,L,1x,L)')&
            & a1, ik, i, j, i+bcabs%il1-1, j+iuref-1, zeroflag(a1), negativeflag(a1)
        end if
      end do
      close(un)
      if (hamsize /= hamsize_max) then
        write(unitout,*) "(bse.f90) [WARNING] Zero occupancy differences and/or&
          & negative fo-fu occured. Reduced size of BSE Hamiltonian from, ",&
          & hamsize_max, "to ", hamsize, " ."
      end if 

      allocate(amap(hamsize, 3))
      allocate(occfactor(hamsize))
      a2=0
      do a1 = 1, hamsize_max
        if(.not. zeroflag(a1) .and. .not. negativeflag(a1)) then
          a2=a2+1
          occfactor(a2) = occfactor_t(a1)
          call hamidx_back(a1, i, j, ik, bcabs%n1, bcabs%n2)
          amap(a2,:) = [ik, i, j]
        end if
      end do

      ! Check consistency 
      check = sum(kkpblocksize)
      if (check /= hamsize) then
        write(*,*) "INCONSITENCY! check=", check, "hamsize=", hamsize
      end if

      deallocate(occfactor_t, zeroflag, negativeflag)

      ! Write out index mapping
      call getunit(un)
      open(un, file='BSE_SINDEX.OUT', action='write', status='replace')
      write(un,'("#",1x,a3,1x,a5,1x,a5,1x,a5)') "s", "ik", "io" ,"iu"
      do a1 = 1, hamsize
        write(un,'(I5,1x,I5,1x,I5,1x,I5)')&
          & a1, amap(a1,1), amap(a1,2), amap(a1,3)
      end do
      close(un)

      ! Write out index mapping for ou at one k
      call getunit(un)
      open(un, file='BSE_OUINDEX.OUT', action='write', status='replace')
      write(un,'("#",1x,a3,1x,a5,1x,a5)') "s", "io" ,"iu"
      do a1 = 1, nou
        call subhamidx_back(a1, i, j, bcabs%n2)
        write(un,'(I5,1x,I5,1x,I5)')&
          & a1, i, j
      end do
      close(un)

    end subroutine checkoccupancies

    subroutine printkouflag(nou, nk, kouflag, kkpblocksize)
      use m_getunit
      implicit none 
      integer, intent(in) :: nou, nk, kkpblocksize(:)
      logical, intent(in) :: kouflag(:,:)

      integer(4) :: i,j, ikk, ik, jk, un

      call getunit(un)
      open(un, file='BSE_KOUFLAG.OUT', action='write', status='replace')
      write(un,'("#",1x,a)') "iou/ik"
      do i = 1, nou
        write(un, '(L)', advance="no") kouflag(i,1)
        do j = 2, nk
          write(un, '(1x,L)', advance="no") kouflag(i,j)
        end do
        write(un,*)
      end do
      close(un)

      call getunit(un)
      open(un, file='BSE_KKPBLOCKSIZE.OUT', action='write', status='replace')
      write(un,'("#",1x,a3,1x,a5,1x,a5,1x,a5,1x,a5)')&
        & "ik", "jk", "ikkp", "rows", "cols"
      do ikk = 1, nk*(nk+1)/2
        call kkpmap(ikk, nk, ik, jk)
        write(un,'(I5,1x,I5,1x,I5,1x,I5,1x,I5)')&
          & ik, jk, ikk, kkpblocksize(ik), kkpblocksize(jk)
      end do
      close(un)
    end subroutine printkouflag

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
