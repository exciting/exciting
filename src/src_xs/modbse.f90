!BOP
! !MODULE: modbse
! !DESCRIPTION:
!   Supporting global variables and routines for the BSE scope.
!
! !REVISION HISTORY:
!   Created 2016 (Aurich)
!EOP   
!BOC
module modbse
  
  use modmpi
  use modinput, only: input
  use modxs, only: bcbs

  implicit none

  ! Band combinations in relative and absolute indices
  Type(bcbs) :: bcou, bcouabs
  ! Reference states for occupied and unoccupied state indices
  integer(4) :: ioref, iuref 
  integer(4) :: no, nu, noo, nuu, nou
  ! Number of k-k' combinations with ik'>=ik
  integer(4) :: nkkp, nk
  ! BSE gap
  real(8) :: egap
  ! Occupation factors for Hamiltonian construction
  real(8), allocatable :: ofac(:)
  ! Maps, maps for everyone !
  integer(4), allocatable, dimension(:,:) :: smap
  integer(4), allocatable, dimension(:) :: kousize
  logical, allocatable, dimension(:,:) :: kouflag

  contains

    !BOP
    ! !ROUTINE: setranges_modxs
    ! !INTERFACE:
    subroutine setranges_modxs(iq)
    ! !USES:
      use modxs, only: istocc0, istocc, istunocc0, istunocc,&
                     & isto0, isto, istu0, istu
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: iq ! Considered q-point index (must be on the not shifted k-grid)
    !
    ! !DESCRIPTION:
    !   A small wrapper for the routine {\tt findocclims}. Used to initialize
    !   occupation limits for $\vec{k}$ and $\vec{k}+\vec{q}$ in some BSE related
    !   routines.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      implicit none
      integer(4), intent(in) :: iq
      
      ! Find occupation limits for k and k+q and set variables in modxs
      call findocclims(iq, istocc0, istocc, istunocc0,&
        & istunocc, isto0, isto, istu0, istu)

    end subroutine setranges_modxs
    !EOC

    !BOP
    ! !ROUTINE: setbcbs_bse
    ! !INTERFACE:
    subroutine setbcbs_bse
    ! !USES:
      use modxs, only: istunocc0, unitout
      use modinput, only: input
    ! !DESCRIPTION:
    !   This routine needs to be called after {\tt setranges\_modxs}.
    !   It sets the band ranges used in BSE calculation is {\tt modbse}.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
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
      ! Number of occupied states
      no = bcou%n1
      nu = bcou%n2
      ! Number of combinations
      noo = no * no
      nuu = nu * nu
      nou = no * nu

      if(rank == 0) then 
        ! Write out state ranges to INFOXS.OUT
        write(unitout,*)
        write(unitout, '("Info(setbcbs): Information on number of states:")')
        write(unitout, '("  Ranges of states used in construction of BSE matrix:")')
        write(unitout, '("    Range of occupied states and number   :", 2i6, 3x, i6)')&
          & bcouabs%il1, bcouabs%iu1, bcouabs%n1
        write(unitout, '("    Range of unoccupied states and number :", 2i6, 3x, i6)')&
          & bcouabs%il2, bcouabs%iu2, bcouabs%n2
        write(unitout, '("  Ranges of states used in construction of BSE matrix&
          & (relative):")')
        write(unitout, '("    Range of occupied states and number   :", 2i6, 3x, i6)')&
          & bcou%il1, bcou%iu1, bcou%n1
        write(unitout, '("    Range of unoccupied states and number :", 2i6, 3x, i6)')&
          & bcou%il2, bcou%iu2, bcou%n2
      end if
    end subroutine setbcbs_bse
    !EOC

    !BOP
    ! !ROUTINE: hamidx
    ! !INTERFACE:
    integer(4) function hamidx(i1, i2, ik, n1, n2)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: i1, i2, ik ! Indices counting from 1 continuously 
    ! integer(4) :: n1, n2     ! Maximum values of i1, i2 respectively
    ! Out:
    ! integer(4) :: hamidx     ! Combined index
    !
    ! !DESCRIPTION:
    !   The function returns a combined index given two band indices
    !   and a $\vec{k}$ index. It is use in the construction of
    !   the BSE Hamiltonian.\\
    !   Map:\\
    !   $\text{hamdix} = i_2 + (i_1 - 1) n_2 + n_2 n_1 (i_k-1)$\\
    !   Notes:\\
    !     $i_2$ is the fastest varying index, followed in order by $i_1$ and $i_k$.\\
    !     All indices are assumed to be counted from 1 onwards continuously.
    ! 
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: i1, i2, ik, n1, n2
      hamidx = i2 + n2 * (i1-1) + n1 * n2 * (ik-1)
    end function hamidx
    !EOC

    !BOP
    ! !ROUTINE: hamidx_back
    ! !INTERFACE:
    subroutine hamidx_back(s, i1, i2, ik, n1, n2)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: s          ! Combined index created with {\tt hamidx}
    ! integer(4) :: n1, n2     ! Maximum values of i1, i2 respectively
    ! Out:
    ! integer(4) :: i1, i2, ik ! Individual indices 
    !
    ! !DESCRIPTION:
    !   The subroutine does the inverse operation of the function {\tt hamidx}.
    ! 
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
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
    !EOC

    !BOP
    ! !ROUTINE: subhamidx
    ! !INTERFACE:
    integer(4) function subhamidx(i1, i2, n2)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: i1, i2     ! Indices counting from 1 continuously 
    ! integer(4) :: n2         ! Maximum value of i2 
    ! Out:
    ! integer(4) :: subhamidx  ! Combined index
    !
    ! !DESCRIPTION:
    !   The function return a combined index given two indices.
    !   It is use in the construction of the BSE Hamiltonian.\\
    !   Map:\\
    !   $\text{hamdix} = i_2 + (i_1 - 1) n_2$\\
    !   Notes:\\
    !     $i_2$ is the fastest varying index, followed in by $i_1$.\\
    !     All indices are assumed to be counted from 1 onwards continuously.
    ! 
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: i1, i2, n2
      subhamidx = i2 + n2 * (i1-1)
    end function subhamidx
    !EOC

    !BOP
    ! !ROUTINE: subhamidx_back
    ! !INTERFACE:
    subroutine subhamidx_back(s, i1, i2, n2)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: s    ! Combined index
    ! integer(4) :: n2   ! Maximum value of i2 
    ! Out:
    ! integer(4) :: i1, i2   ! Individual indices
    !
    ! !DESCRIPTION:
    !   The routine performs the inverse operation to 
    !   {\tt subhamidx}.
    ! 
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: s, n2
      integer(4), intent(out) :: i1, i2
      i1 = (s-1)/n2 + 1
      i2 = s - (i1-1)*n2
    end subroutine subhamidx_back
    !EOC

    !BOP
    ! !ROUTINE: checkoccupancies
    ! !INTERFACE:
    subroutine checkoccupancies(iq, hamsize)
    ! !USES:
      use modinput, only: input
      use modxs, only: unitout, bcbs
      use m_getunit
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: iq  ! q-point index (on unshifted k-mesh)
    ! Out:
    ! integer(4) :: hamsize            ! Dimension of the BSE Hamiltonian matrix
    ! Module out:
    ! real(8)    :: ofac(hamsize)      ! Occupation factors need for 
    !                                  ! the construction of the BSE matrix
    ! integer(4) :: smap(hamsize, 5)   ! Map between BSE matrix index and k,o,u indices
    ! integer(4) :: kousize(nkkp)      ! How many o u combinations allowed for each k
    ! logical    :: kouflag(no*nu, nk) ! False for o u k combinations that are filtered out
    !
    ! !DESCRIPTION:
    !   This routine checks the requested ik,io,iu combinations whether they
    !   contain problematic occupancy differences and sorts them out if need be.
    !   The simple treatment of fractional occupancy does not allow for transitions
    !   between states of the same partial occupancy. Also cases of occupancy inversion
    !   where the occupancy difference is negative are filtered out, since those break
    !   any kind of hermiticity of the BSE Hamiltonian.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      implicit none 

      integer(4), intent(in) :: iq
      integer(4), intent(out) :: hamsize

      integer(4) :: un
      integer(4) :: ik, i, j, a1, a2, hamsize_max, check
      integer(4) :: iul, iuu, iol, iou
      real(8) :: kdocc, cutoff
      real(8), allocatable :: occfactor_t(:), docc(:,:)
      logical, allocatable :: zeroflag(:), negativeflag(:)
      character(256) :: frmt

      ! Band ranges
      iol = bcouabs%il1
      iou = bcouabs%iu1
      iul = bcouabs%il2
      iuu = bcouabs%iu2

      ! Maximal size of Hamiltonian
      ! Size of BSE-Hamiltonian = #o * #u * #k
      hamsize_max = nou * nk

      ! What is considered to be occupied
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
      allocate(docc(no,nu))
      allocate(kouflag(nou, nk))
      kouflag = .true.


      ! Check occupancies
      do ik = 1, nk
      
        ! Get occupation numbers 
        !   Calculates docc(l,m)= f_{o_l ki+q} - f_{u_m ki}
        !   Remark: Currently the iq argument does nothing at all
        call getdocc(iq, ik, ik, iol, iou, iul, iuu, docc)
        
        ! Band indices
        do i = iol, iou ! Occupied (absolute)
          do j = iul, iuu ! Unoccupied (absolute)

            ! Map to combined index
            a1 = hamidx(i-iol+1, j-iul+1, ik, no, nu)

            ! Get occupation difference for current combination
            ! kdocc = f_{o_{a},k_{a}} - f_{u_{a},k_{a}}
!! This assumes KS calculations with max. occupancies of 2
            kdocc = docc(i-iol+1, j-iul+1)/2.0d0
            ! Set zero occupation difference flag
            if(abs(kdocc) .lt. cutoff) then
              zeroflag(a1) = .true.
              kouflag(subhamidx(i-iol+1,j-iul+1, nu), ik) = .false.
            end if
            ! Set occupation difference sign flag
            if(kdocc .le. 0.0d0) then
              negativeflag(a1) = .true.
              kouflag(subhamidx(i-iol+1,j-iul+1, nu), ik) = .false.
            end if
            ! Occupation factor
            occfactor_t(a1) = sqrt(abs(kdocc))
          end do
        end do

      end do
      deallocate(docc)

      allocate(kousize(nk))
      do ik = 1, nk
        ! Number of io iu combinations valid for ik
        kousize(ik) = count(kouflag(:,ik))
      end do
      if(input%xs%dbglev > 2 .and. rank == 0) then 
        ! Print kouflag and kousize
        call printkouflag(nou, nk, kouflag, kousize)
        ! Print all occupancy factors with flags
        call printoccupancies(zeroflag, negativeflag, occfactor_t)
      end if

      ! Adjust selection for only non zero occupancy differences
      ! and only positive occupancy differences.
      ! Write out skipped band combinations
      ! and make map of combined indices to use.
      if(rank == 0) then
        call getunit(un)
        open(un, file='BSE_SKIPPED_BCBS.OUT', action='write', status='replace')
        write(un,'("#",1x,a)') "Skipped band combinations"
        frmt='("#",a7,1x,a4,1x,a4,1x,a4,1x,a4,1x,a4,1x,a1,1x,a1)'
        write(un,frmt) "s", "ik", "io", "iu", "ioa", "iua", "z", "n"
      end if

      hamsize = 0
      do a1 = 1, hamsize_max
        if(.not. zeroflag(a1) .and. .not. negativeflag(a1)) then
          hamsize = hamsize + 1
        else
          if(rank == 0) then
          call hamidx_back(a1, i, j, ik, bcouabs%n1, bcouabs%n2)
            write(un,'(I8,1x,I4,1x,I4,1x,I4,1x,I4,1x,I4,1x,L,1x,L)')&
              & a1, ik, i, j, i+iol-1, j+iul-1, zeroflag(a1), negativeflag(a1)
          end if
        end if
      end do

      if(rank == 0) then
        close(un)
        if (hamsize /= hamsize_max) then
          write(unitout,*) "(bse.f90) [WARNING] Zero occupancy differences and/or&
            & negative fo-fu occured. Reduced size of BSE Hamiltonian from, ",&
            & hamsize_max, "to ", hamsize, " ."
        end if 
      end if

      allocate(smap(hamsize, 5))
      allocate(ofac(hamsize))
      a2=0
      do a1 = 1, hamsize_max
        if(.not. zeroflag(a1) .and. .not. negativeflag(a1)) then
          a2=a2+1
          ofac(a2) = occfactor_t(a1)
          call hamidx_back(a1, i, j, ik, bcouabs%n1, bcouabs%n2)
          smap(a2,:) = [ik, i, j, i+iol-1, j+iul-1]
        end if
      end do

      ! Check consistency 
      if(rank == 0) then 
        check = sum(kousize)
        if (check /= hamsize) then
          write(*,*) "INCONSITENCY! check=", check, "hamsize=", hamsize
          call terminate
        end if
      end if

      deallocate(occfactor_t, zeroflag, negativeflag)

      ! Write out index mapping
      if(rank == 0) then 
        call getunit(un)
        open(un, file='BSE_SINDEX.OUT', action='write', status='replace')
        write(un,'("#",1x,a3,1x,a5,1x,a5,1x,a5)') "s", "ik", "io" ,"iu"
        do a1 = 1, hamsize
          write(un,'(I5,1x,I5,1x,I5,1x,I5)')&
            & a1, smap(a1,1), smap(a1,4), smap(a1,5)
        end do
        close(un)
      end if

    end subroutine checkoccupancies
    !EOC


    ! Auxiliary routines
    subroutine printkouflag(nou, nk, kouflag, kousize)
      use m_getunit
      implicit none 
      integer, intent(in) :: nou, nk, kousize(:)
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
          & ik, jk, ikk, kousize(ik), kousize(jk)
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
!EOC
