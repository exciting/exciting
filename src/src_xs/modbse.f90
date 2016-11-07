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
  use modxs, only: evalsv0, occsv0, ikmapikq
  use mod_eigenvalue_occupancy, only: evalsv, occsv
  use modmpi
  use modinput, only: input

  implicit none

  ! Reference states for occupied and unoccupied state indices
  integer(4) :: ioref, iuref 
  ! Size of the resonant-resonant block of the Hamiltonian
  integer(4) :: hamsize
  ! Number of all k-points and k-points/k-k'-combinations
  ! to be considered in BSE
  integer(4) :: nk_max, nk_bse, nkkp_bse
  ! Number of occupied and unoccupied bands present
  integer(4) :: no_max, nu_max, nou_max
  ! Maximum of needed number of o/u bands over all needed k
  integer(4) :: no_bse_max, nu_bse_max, nou_bse_max
  ! Energy range of spectrum
  real(8) :: wl, wu
  integer(4) :: nw
  ! Cutoff for occupation
  real(8) :: cutoffocc
  ! Convergence energy
  real(8) :: econv, ewidth
  ! BSE gap
  real(8) :: egap, evalshift
  ! Occupation factors for Hamiltonian construction
  real(8), allocatable :: ofac(:)
  ! Combined BSE index map
  integer(4), allocatable :: smap(:,:)
  integer(4), allocatable :: smap_rel(:,:)
  ! Contributing k points
  ! relative index ik -> global index iknr
  integer(4), allocatable :: kmap_bse_rg(:)
  ! global index iknr -> relative index ik
  integer(4), allocatable :: kmap_bse_gr(:)
  ! Number of allowed o-u transitions at each k-point
  integer(4), allocatable :: kousize(:)
  ! Lower and upper bound for occupied and unoccupied band indices for each k
  integer(4), allocatable :: koulims(:,:)

  ! Whether or not to read eigenvalues etc. form file
  logical :: read_eval_occ_qmt = .true.
  logical :: seltrans = .true.

  ! Filebasenames
  character(256) :: scclifbasename = "SCCLI"
  character(256) :: exclifbasename = "EXCLI"
  character(256) :: scclicfbasename = "SCCLIC"
  character(256) :: exclicfbasename = "EXCLIC"
  character(256) :: scclifname, scclicfname
  character(256) :: exclifname, exclicfname

  ! Legacy
  ! GW eigenvalue backup 
  real(8), allocatable, dimension(:,:) :: eval0

  contains

    !+++++++++++++++++++++++++++++++++++++++++++!
    ! Routines to setup the combinded index of  ! 
    ! the BSE hamiltonian.                      !
    !+++++++++++++++++++++++++++++++++++++++++++!

    !BOP
    ! !ROUTINE: setranges_modxs
    ! !INTERFACE:
    subroutine setranges_modxs(iq)
    ! !USES:
      use modxs, only: istocc0, istocc, istunocc0, istunocc,&
                     & isto0, isto, istu0, istu, ksgapval
      use mod_kpoint, only: nkptnr
      use mod_eigenvalue_occupancy, only: nstsv
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: iq ! Considered q-point index (must be on the not shifted k-grid)
    !
    ! !DESCRIPTION:
    !   A small wrapper for the routine {\tt findocclims}. Used to initialize
    !   modxs module variables for occupation limits for $\vec{k}$ and $\vec{k}+\vec{q}$ 
    !   in some BSE related routines. 
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      implicit none
      integer(4), intent(in) :: iq
      
      ! Find occupation limits for k and k+q and set variables in modxs
      ! This also reads in 
      ! (QMTXXX)
      ! mod_eigevalue_occupancy:evalsv, mod_eigevalue_occupancy:occsv 
      ! (QMT000)
      ! modxs:evalsv0, modxs:occsv0
      ! and determines the ks gap modxs:ksgap.
      call findocclims(iq, istocc0, istocc, istunocc0,&
        & istunocc, isto0, isto, istu0, istu)

      ! Set read flag
      read_eval_occ_qmt = .false.

      ! Set additional modbse variables
      iuref = istunocc0 
      ioref = 1
      nk_max = nkptnr
      no_max = istocc0
      nu_max = nstsv-istunocc0+1
      nou_max = no_max*nu_max
      egap = ksgapval

    end subroutine setranges_modxs
    !EOC

    !BOP
    ! !ROUTINE: select_transitions
    ! !INTERFACE:
    subroutine select_transitions(iqmt)
    ! !USES:
      use modinput, only: input
      use modxs, only: unitout, istocc0, istunocc0
      use mod_eigenvalue_occupancy, only: nstsv
      use m_getunit
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: iqmt  ! q-point index (on unshifted k-mesh)
    ! Module out:
    ! integer(4) :: hamsize            ! Dimension of the BSE Hamiltonian matrix
    ! real(8)    :: ofac(hamsize)      ! Occupation factors need for 
    !                                  ! the construction of the BSE matrix
    ! integer(4) :: smap(hamsize, 3)   ! Map between BSE matrix index and u,o,k indices
    ! integer(4) :: kousize(nk)        ! How many u o combinations allowed for each k
    !
    ! !DESCRIPTION:
    !   Given a selected energy range for the spectrum, this routine will 
    !   select relevant transitions for each k point. Appart from the KS transition
    !   energies the routine checks whether the transition contain problematic 
    !   occupancy differences and sorts them out if need be.
    !   The simple treatment of fractional occupancy does not allow for transitions
    !   between states of the same partial occupancy. Also cases of occupancy inversion
    !   where the occupancy difference is negative are filtered out, since those break
    !   any kind of hermiticity of the BSE Hamiltonian.\\
    !   The routine crates the compined index map
    !   $\alpha \leftrightarrow \{u_\alpha, o_\alpha, \vec{k}_\alpha\}$, auxilliary maps
    !   and determins the size of the resulting hamiltonian.
    !   In all cases is u the fastest index followed by o and k.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      implicit none 

      integer(4), intent(in) :: iqmt

      integer(4) :: ik, ikq, s, iknr
      integer(4) :: io, iu, kous
      integer(4) :: iomax, iumax
      integer(4) :: iomin, iumin
      real(8), parameter :: maxocc = 2.0d0

      integer(4) :: nk_loc, hamsize_loc
      integer(4), allocatable :: smap_loc(:,:)
      real(8), allocatable :: ofac_loc(:)
      logical, allocatable :: sflag(:)
      integer(4) :: buflen, buflen_r
      integer(4) :: k1, k2, k1_r, k2_r
      integer(4) :: iproc, i1, i2, il_r, iu_r

      ! What energy range is of interest?
      nw = input%xs%energywindow%points
      wl = input%xs%energywindow%intv(1)
      wu = input%xs%energywindow%intv(2)

      ! Exciton energies determine spectrum.
      ! Respect broadening of peaks and include additional 
      ! excition solutions within of a lorentian boradening decay of 1/100 
      ewidth = input%xs%broad*sqrt(99.0d0)

      ! Excitons are build from KS states and the KS transition energies 
      ! dominatly determine exciton energies. 
      ! Select KS transitions withing the energy energy window for the 
      ! spectrum plus extra convergence energy.
      econv = abs(input%xs%bse%econv)
      econv = econv + ewidth

      ! What is considered to be occupied
      cutoffocc = input%groundstate%epsocc

      ! Sizes local/maximal
      nk_loc = ceiling(real(nk_max,8)/real(mpiglobal%procs,8))
      hamsize_loc = nk_loc*nou_max

      ! The index mapping we want to build 
      ! s(1) = iuabs, s(2) = ioabs, s(3) = iknr
      allocate(smap_loc(3, hamsize_loc))

      ! Flag map whether to use k-o-u combination or not
      allocate(sflag(hamsize_loc))
      sflag = .false.

      ! Occupation factor
      allocate(ofac_loc(hamsize_loc))

      ! How many o-u combinations at ik 
      allocate(kousize(nk_max))
      ! and limits of band ranges
      allocate(koulims(4,nk_max))
           
      ! Loop over kpoints (non-reduced)
      !   Distribute k-loop over global MPI communicator.
      !   Each participating rank gets a continuous ik interval.
      !   Not participating ranks have k1=0 and k2=-1.
      k1 = firstofset(mpiglobal%rank, nk_max, mpiglobal%procs)
      k2 = lastofset(mpiglobal%rank, nk_max, mpiglobal%procs)

      ikloop: do ik = k1, k2

        ! Get k+q index form k index
        ikq = ik
        if(iqmt .ne. 0) ikq = ikmapikq(ik, iqmt)

        kous = 0
        iomax = 0
        iumax = 0
        iomin = istocc0+1
        iumin = nstsv+1

        ! Loop over KS transition energies 
        !$OMP PARALLEL DO &
        !$OMP& COLLAPSE(2),&
        !$OMP& DEFAULT(SHARED), PRIVATE(io,iu,s),&
        !$OMP& REDUCTION(+:kous),&
        !$OMP& REDUCTION(max:iomax),&
        !$OMP& REDUCTION(min:iomin),&
        !$OMP& REDUCTION(max:iumax),&
        !$OMP& REDUCTION(min:iumin)
        do io = 1, istocc0
          do iu = istunocc0, nstsv 
            
            ! Only consider transitions which are in the energy window
            ! \Delta E = \epsilon_{u ki+q} - \epsilon_{o ki}
            if( evalsv(iu, ikq) - evalsv0(io, ik) <= wu+econv&
              & .and. evalsv(iu, ikq) - evalsv0(io, ik) >= max(wl-econv,0.0d0) ) then

              ! Only consider transitions which have a positve non-zero 
              ! occupancy difference f_{o ki} - f_{u ki+q}
              if( occsv0(io, ik) - occsv(iu, ikq) > cutoffocc) then 

                ! Combine u, o and k index
                ! u is counted from lumo=1 upwards
                ! o is counted from homo=1 downwards
                ! ik index is shifted, due to MPI parallelization
                s = hamidx(iu-istunocc0+1, istocc0-io+1, ik-k1+1, nu_max, no_max)

                ! Use that u-o-k combination
                sflag(s) = .true.
                
                ! Write to combinded index map
                smap_loc(:,s) = [iu, io, ik]

                ! Save occupation factor
                ofac_loc(s) = sqrt((occsv0(io, ik) - occsv(iu, ikq))/maxocc)

                ! Keep track of how many valid transitions
                ! are considered at current k point.
                kous = kous + 1

                ! Keep track of the minimal/maximal io/iu
                iomax = max(iomax, io)
                iomin = min(iomin, io)
                iumax = max(iumax, iu)
                iumin = min(iumin, iu)

              end if
              
            end if

          end do
        end do
        !$OMP END PARALLEL DO

        koulims(:,ik) = [ iumin, iumax, iomin, iomax ] 

        kousize(ik) = kous
        
      end do ikloop

#ifdef MPI
      ! Collect kousize on all processes 
      call mpi_allgatherv_ifc(set=nk_max, rlen=1, ibuf=kousize,&
        & inplace=.true., comm=mpiglobal)
      ! Collect koulims on all processes 
      call mpi_allgatherv_ifc(set=nk_max, rlen=4, ibuf=koulims,&
        & inplace=.true., comm=mpiglobal)
#endif
      ! Calculate maximal no(k) and nu(k)
      no_bse_max=0
      nu_bse_max=0
      do ik = 1, nk_max
        no_bse_max = max(koulims(4, ik)-koulims(3, ik)+1, no_bse_max)
        nu_bse_max = max(koulims(2, ik)-koulims(1, ik)+1, nu_bse_max)
      end do
      nou_bse_max=no_bse_max*nu_bse_max
      
      ! Global results
      ! Number of contributing k points
      nk_bse = count(kousize /= 0)
      nkkp_bse = nk_bse*(nk_bse+1)/2
      ! List of participating k-points
      allocate(kmap_bse_rg(nk_bse))
      allocate(kmap_bse_gr(nk_max))
      ik = 0
      do iknr = 1, nk_max
        if(kousize(iknr) /= 0) then
          ik = ik + 1
          kmap_bse_rg(ik) = iknr
          kmap_bse_gr(iknr) = ik
        else
          kmap_bse_gr(iknr) = 0
        end if
      end do

      ! Size of the resulting resonant BSE Hamiltonian
      hamsize = sum(kousize)
      ! Combined index map
      allocate(smap(3,hamsize))
      ! Occupation factors
      allocate(ofac(hamsize))

      ! If rank was participating in k-loop.
      if(k2 > 0) then 
        ! Apply selection flag to local arrays
        ! and store result in global counterparts.
        i1 = sum(kousize(1:k1-1))+1
        i2 = sum(kousize(k1:k2))+i1-1
        smap(1,i1:i2) = pack(smap_loc(1,:),sflag)
        smap(2,i1:i2) = pack(smap_loc(2,:),sflag)
        smap(3,i1:i2) = pack(smap_loc(3,:),sflag)
        ofac(i1:i2) = pack(ofac_loc,sflag)
      end if

      ! Local auxiliary local arrays not needed anymore
      deallocate(smap_loc)
      deallocate(ofac_loc)
      deallocate(sflag)

#ifdef MPI
      ! Collect ofac on all processes 
      call mpi_allgatherv_ifc(set=nk_max, rlenv=kousize, rbuf=ofac,&
        & inplace=.true., comm=mpiglobal)
      ! Collect smap on all processes 
      call mpi_allgatherv_ifc(set=nk_max, rlenv=kousize*3, ibuf=smap,&
        & inplace=.true., comm=mpiglobal)
#endif
      
      ! Make relative combinded index map
      allocate(smap_rel(3,hamsize))
      do i1 = 1, hamsize
        iknr = smap(3, i1)
        smap_rel(1,i1) = smap(1, i1) - koulims(1, iknr) + 1
        smap_rel(2,i1) = smap(2, i1) - koulims(3, iknr) + 1
        smap_rel(3,i1) = kmap_bse_gr(iknr)
      end do

      if(rank == 0) then 
        call printso
      end if

      ! Set flag 
      seltrans = .false.

    end subroutine select_transitions
    !EOC

    subroutine printso
      use m_getunit
      implicit none 

      integer(4) :: i, un

      call getunit(un)
      open(un, file='BSE_SINDEX.OUT', action='write', status='replace')
      write(un,'("# Combined BSE index @ q =", 3(E10.3,1x))')  0.0d0, 0.0d0, 0.0d0
      write(un,'("# s ik io iu ik_rel io_rel iu_rel occ")')
      do i = 1, size(ofac)
        write(un, '(7(I8,1x),1x,E23.16)')&
          & i, smap(i,:), smap_rel(i,:), ofac(i)
      end do
      close(un)

      call getunit(un)
      open(un, file='KOU.OUT', action='write', status='replace')
      write(un,'("# k-o-u ranges used in combined BSE index @ q =", 3(E10.3,1x))')&
        &  0.0d0, 0.0d0, 0.0d0
      write(un,'("# ik iu1 iu2 io1 io2 nou")')
      do i = 1, nk_max
        write(un, '(6(I8,1x))')&
          & i, koulims(i,1), koulims(i,2), koulims(i,3), koulims(i,4), kousize(i)
      end do
      close(un)

    end subroutine printso

    !+++++++++++++++++++++++++++++++++++++++++++!
    ! Simple continuous combined index mappers  !
    !+++++++++++++++++++++++++++++++++++++++++++!

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

end module modbse
!EOC
