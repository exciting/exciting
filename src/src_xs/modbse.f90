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
  
  use m_filedel
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
  ! Energy range of spectrum
  real(8) :: wl, wu
  integer(8) :: nw
  ! BSE gap
  real(8) :: evalshift
  ! GW eigenvalue backup 
  real(8), allocatable, dimension(:,:) :: eval0
  ! Occupation factors for Hamiltonian construction
  real(8), allocatable :: ofac(:)
  ! Maps, maps for everyone !
  integer(4), allocatable, dimension(:,:) :: smap
  integer(4), allocatable, dimension(:) :: kousize
  logical, allocatable, dimension(:,:) :: kouflag
  ! Filenames
  character(256) :: scclifname = "SCCLI.OUT"
  character(256) :: exclifname = "EXCLI.OUT"
  character(256) :: scclicfname = "SCCLIC.OUT"
  character(256) :: exclicfname = "EXCLIC.OUT"

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
        write(unitout,*)
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
    !   $\text{hamdix} = i_1 + (i_2 - 1) n_1 + n_1 n_2 (i_k-1)$\\
    !   Notes:\\
    !     $i_1$ is the fastest varying index, followed in order by $i_2$ and $i_k$.\\
    !     All indices are assumed to be counted from 1 onwards continuously.
    ! 
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Changed fastes index to i1. (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: i1, i2, ik, n1, n2
      hamidx = i1 + n1 * (i2-1) + n1 * n2 * (ik-1)
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
      i2 = (tmp-1)/n1 + 1
      i1 = tmp - (i2-1)*n1
    end subroutine hamidx_back
    !EOC

    !BOP
    ! !ROUTINE: subhamidx
    ! !INTERFACE:
    integer(4) function subhamidx(i1, i2, n1)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: i1, i2     ! Indices counting from 1 continuously 
    ! integer(4) :: n1         ! Maximum value of i1 
    ! Out:
    ! integer(4) :: subhamidx  ! Combined index
    !
    ! !DESCRIPTION:
    !   The function return a combined index given two indices.
    !   It is use in the construction of the BSE Hamiltonian.\\
    !   Map:\\
    !   $\text{hamdix} = i_1 + (i_2 - 1) n_1$\\
    !   Notes:\\
    !     $i_1$ is the fastest varying index, followed by $i_2$.\\
    !     All indices are assumed to be counted from 1 onwards continuously.
    ! 
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !   Changed fastest index to i1. (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: i1, i2, n1
      subhamidx = i1 + n1 * (i2-1)
    end function subhamidx
    !EOC

    !BOP
    ! !ROUTINE: subhamidx_back
    ! !INTERFACE:
    subroutine subhamidx_back(s, i1, i2, n1)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: s    ! Combined index
    ! integer(4) :: n1   ! Maximum value of i1 
    ! Out:
    ! integer(4) :: i1, i2   ! Individual indices
    !
    ! !DESCRIPTION:
    !   The routine performs the inverse operation to 
    !   {\tt subhamidx}.
    ! 
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !   Changed fastest index to i1. (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: s, n1
      integer(4), intent(out) :: i1, i2
      i2 = (s-1)/n1 + 1
      i1 = s - (i2-1)*n1
    end subroutine subhamidx_back
    !EOC

    !BOP
    ! !ROUTINE: select_transitions
    ! !INTERFACE:
    subroutine select_transitions(iq, hamsize)
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
    ! integer(4) :: kousize(nk)        ! How many o u combinations allowed for each k
    !
    ! !DESCRIPTION:
    !   Given an selected energy range for the spectrum, this routine will 
    !   select relevant transitions for each k point. Appart from the KS transition
    !   energies the routine checks whether the transition contain problematic 
    !   occupancy differences and sorts them out if need be.
    !   The simple treatment of fractional occupancy does not allow for transitions
    !   between states of the same partial occupancy. Also cases of occupancy inversion
    !   where the occupancy difference is negative are filtered out, since those break
    !   any kind of hermiticity of the BSE Hamiltonian.\\
    !   The routine crates the compined index map
    !   $\alpha \leftrightarrow \{u, o, \vec{k}\}$, determins the size of 
    !   the resulting hamiltonian and auxilliary maps.
    !   In all cases u is the fastest index followed by o and k.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      implicit none 

      integer(4), intent(in) :: iq
      integer(4), intent(out) :: hamsize

      integer(4) :: ik, s
      integer(4) :: io, iu
      real(8) :: cutoff, econv
      real(8), parameter :: maxocc = 2.0d0
      real(8), allocatable :: docc(:,:)

      integer(4) :: nk_loc, hamsize_loc
      integer(4) :: no_max, nu_max, nou_max
      integer(4), allocatable :: smap_loc(:,:), kousize(:)
      real(8), allocatable :: ofac_loc(:)
      logical, allocatable :: sflag(:)
      integer(4) :: buflen, buflen_r
      integer(4) :: il, iu, il_r, iu_r

      ! Considered energy window plus extra convergence energy
      econv = abs(input%xs%bse%econv)/h2ev

      ! What is considered to be occupied
      cutoff = input%groundstate%epsocc

      ! Sizes local/maximal
      nk_loc = ceiling(real(nk,8)/real(mpiglobal%procs,8))
      no_max = istocc0
      nu_max = nstsv-istunocc0+1
      nou_max = no_max*nu_max 
      hamsize_loc = nk_loc*nou_max

      ! The index mapping we want to build 
      ! s(1) = iuabs, s(2) = ioabs, s(3) = iknr
      allocate(smap_loc(hamsize_loc, 3))

      ! Flag map whether to use k-o-u combination or not
      allocate(sflag(hamsize_loc))
      sflag = .false.

      ! Occupation factor
      allocate(ofac_loc(hamsize_loc))

      ! How many o-u combinations at ik 
      allocate(kousize(nk))
      kousize = 0
           
      ! Occupation differences 
      ! Maximum occupancy is set by maxocc parameter 
      allocate(docc(istocc0,nstsv-istunocc0+1))

      ! Loop over kpoints (non-reduced)
      !   MPI distributed loop over k points (continuous k on one thread)
      call genparidxran('k', nk_max)

      buflen = 0
      ik: do ik = kpari, kparf

        ! Get occupation numbers 
        !   Calculates docc(l,m)= f_{o_l ki} - f_{u_m ki+q}
        !   Remark: Currently the iq argument does nothing at all
        call getdocc(iq, ik, ik, 1, istocc0, istunocc0, nstsv, docc)

        kous = 0

        ! Loop over KS transition energies (q = 0) 
        !$OMP PARALLEL DO 
        !$OMP& COLLAPSE(2)
        !$OMP& DEFAULT(SHARED) PRIVATE(io,iu,s)
        !$OMP& REDUCTION(+:kous)
        do io = istocc0, 1, -1 
          do iu = istunocc0, nstsv 
            
            ! Only consider transitions which are in the energy window
            if( evalsv(iu, ik) - evalsv(io, ik) <= wu+econv&
              & .and. evalsv(iu, ik) - evalsv(io, ik) >= max(wl-econv,0.0d0) ) then

              ! Only consider transitions which have a positve non-zero 
              ! occupancy difference
              if( docc(io, iu-istunocc0+1) > cutoff) then 

                ! Shift global index to local index
                s = hamidx(iu-istunocc0+1, io, ik, nu_max, no_max) - (kpari-1)*nou_max

                sflag(s) = .true.

                smap_max(s,1) = iu
                smap_max(s,2) = io
                smap_max(s,3) = ik

                ofac_max(s) = sqrt(docc(io, iu-istunocc0+1)/maxocc)

                kous = kous + 1

              end if
              
            end if

          end do
        end do
        !$OMP END PARALLEL DO

        kousize(ik) = kous
        
        buflen = buflen + kous

      end do ik

      deallocate(docc)

#ifdef MPI
      ! Broadcast kousize entrys to every one
      do iproc=0, mpiglobal%procs-1
        if(mpiglobal%rank == iproc) then
          if(kparf > 0) then 
            call mpi_bcast(kousize(kpari:kparf), kparf-kpari+1, mpi_integer, iproc,&
              & mpiglobal%comm, mpiglobal%ierr)
          end if
        else
          kpari_r = firstofset(iproc, nk_max)
          kparf_r = lastofset(iproc, nk_max)
          if(kparf_r > 0) then 
            call mpi_bcast(kousize(kpari_r:kparf_r), kparf_r-kpari_r+1,&
              & mpi_integer, iproc, mpiglobal%comm, mpiglobal%ierr)
          end if
        end if
      end do
#endif
      
      ! Collect results
      hamsize = sum(kousize)
      allocate(smap(hamsize, 3))
      allocate(ofac(hamsize))

      if(kparf > 0) then 
        ! Apply selection flags / prepare send buffers
        il = sum(kousize(1:kpari-1))+1
        iu = sum(kousize(kpari:kparf))+il-1
        smap(il:iu,1) = pack(smap_loc(:,1),sflag)
        smap(il:iu,2) = pack(smap_loc(:,2),sflag)
        smap(il:iu,3) = pack(smap_loc(:,3),sflag)
        ofac(il:iu) = pack(ofac_loc,sflag)
      end if

      deallocate(smap_loc)
      deallocate(ofac_loc)
      deallocate(sflag)

#ifdef MPI
      ! Broadcast results
      do iproc=0, mpiglobal%procs-1
        ! Send local data
        if(mpiglobal%rank == iproc) then
          ! Only if local procuded something
          if(kparf > 0) then 
            call mpi_bcast(smap(il,1), buflen, mpi_integer, iproc,&
              & mpiglobal%comm, mpiglobal%ierr)
            call mpi_bcast(smap(il,2), buflen, mpi_integer, iproc,&
              & mpiglobal%comm, mpiglobal%ierr)
            call mpi_bcast(smap(il,3), buflen, mpi_integer, iproc,&
              & mpiglobal%comm, mpiglobal%ierr)
            call mpi_bcast(ofac(il), buflen, mpi_double_precision, iproc,&
              & mpiglobal%comm, mpiglobal%ierr)
          end if
        ! Receive remote data
        else
          kpari_r = firstofset(iproc, nk_max)
          kparf_r = lastofset(iproc, nk_max)
          ! Only if remote procued something
          if(kparf_r > 0) then 
            il_r = sum(kousize(1:kpari_r-1))+1
            iu_r = sum(kousize(kpari_r:kparf_r))+il_r-1
            buflen_r = iu_r-il_r+1
            call mpi_bcast(smap(il_r,1), buflen_r, mpi_integer, iproc,&
              & mpiglobal%comm, mpiglobal%ierr)
            call mpi_bcast(smap(il_r,2), buflen_r, mpi_integer, iproc,&
              & mpiglobal%comm, mpiglobal%ierr)
            call mpi_bcast(smap(il_r,3), buflen_r, mpi_integer, iproc,&
              & mpiglobal%comm, mpiglobal%ierr)
            call mpi_bcast(ofac(il_r), buflen_r, mpi_double_precision, iproc,&
              & mpiglobal%comm, mpiglobal%ierr)
          end if
        end if
      end do
#endif

      if(rank == 0) then 
        call printso
      end if

    end subroutine select_transitions
    !EOC

    subroutine printso
      use m_getunit
      implicit none 

      integer(4) :: i, un

      call getunit(un)
      open(un, file='BSE_SINDEX.OUT', action='write', status='replace')
      write(un,'("# Combined BSE index @ q =", 3(E10.3,1x))')  0.0d0, 0.0d0, 0.0d0
      write(un,'("# s ik io iu occ")')
      do i = 1, size(ofac)
        write(un, '(4(I8,1x),1x,E23.16)')&
          & i, smap(i,3), smap(i,2), smap(i,1), ofac(i)
      end do
      close(un)

    end subroutine printso

end module modbse
!EOC
