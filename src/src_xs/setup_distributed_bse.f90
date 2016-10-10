!BOP
! !ROUTINE: setup_distributed_bse
! !INTERFACE:
subroutine setup_distributed_bse(ham)
! !USES:
  use modmpi
  use modbse
  use modsclbse
  use modinput
  use m_getunit
! !INPUT/OUTPUT PARAMETERS:
! In/Out:
! type(dzmat) :: ham  ! 2D block cyclic distributed BSE-Hamiltonian matrix
! 
! !DESCRIPTION:
!   The routine sets up the content of the local array of the 
!   2d block cyclic distributed BSE-Hamiltonian matrix.
!
! !REVISION HISTORY:
!   Created. (Aurich)
!EOP
!BOC

  !! I/O
  type(dzmat), intent(inout) :: ham

  !! Local variables
  integer(4) :: ikkp, ikkp_loc
  integer(4) :: ik1, ik2, ik1_loc, ik2_loc
  integer(4) :: mytarget, mysource
  integer(4) :: id, iprc, ishift, a1, a2, i, j
  integer(4) :: ndj, nall, rest
  integer(4), allocatable :: myreadlist(:), myreceivelist(:,:), mysendlist(:,:)
  complex(8), allocatable, dimension(:,:,:,:) :: excli, sccli

  ! External function
!#ifdef SCAL
!  integer(4), external :: indxl2g
!#endif

  integer(4) :: un
  logical :: verbose
  character(256) :: fnam

  verbose = .true.

  ! Allocate work arrays
  allocate(excli(no,nu,no,nu)) ! RR part of V
  allocate(sccli(no,nu,no,nu)) ! RR part of W

  ! Read cycles
  !   The nkkp registers of W and V are read by npoc processes
  !   in N or N+1 read cycles. Where nkkp/nproc = N rest R, so
  !   that in N cycles every process reads a register, while in
  !   the N+1'th and last cycle only the first R processes read 
  !   from file.
  !   After the read is done, the data gets shifted nporc-1 times
  !   from neighbour to next neighbour in one direction (along 
  !   a process ring, in mathematically negative direction),
  !   since each process needs data form every 
  !   file register to build the local Hamiltonian matrix.
  ndj  = ceiling(dble(nkkp)/dble(nproc)) 
  nall = nkkp/nproc
  rest = mod(nkkp, nproc)

  ! Determine from which rank to receive data (always the process previous in ring)
  mysource = mypcol1d - 1 
  if(mysource < 0) mysource = nproc -1 

  ! Determine to which rank to send data (always next process in ring)
  mytarget = mypcol1d + 1
  if(mytarget > nproc-1) mytarget = 0 

  if(verbose) then
    if(rank == 0) then
      write(*,'(" Number of data injections: ", i3, 1x, i3, 1x, i3)'), ndj, nall, rest
    end if
    write(*,'(" Rank", i3," Rank1d", i3 " Form rank:", i3, " To rank:", i3)') rank, mypcol1d, mysource, mytarget
  end if

  ! The registers to read by each process are distributed
  ! around the process ring one by one.
  allocate(myreadlist(ndj))
  do id = 1, ndj
    myreadlist(id) = mypcol1d + (id-1)*nproc + 1
    ! An entry of 0 means no read from file
    ! in cycle N+1.
    if(myreadlist(id) > nkkp) then
      myreadlist(id) = 0
    end if
  end do

  !!*************************************************!!
  !! Make lists of which ikkps to receive each cycle !!
  !! and which to send.                              !!
  !!*************************************************!!
  allocate(myreceivelist(ndj, nproc))
  allocate(mysendlist(ndj, nproc))
  ! If there are more ikkp then processes:
  ! All processes read new data and send it around
  do id = 1, nall

    myreceivelist(id, 1) = 0
    mysendlist(id, 1) = myreadlist(id)

    do ishift = 2, nproc
      myreceivelist(id, ishift) = mysendlist(id,ishift-1) - 1
      if(mod(myreceivelist(id, ishift), nproc) == 0) then
        myreceivelist(id, ishift) = myreceivelist(id, ishift) + nproc
      end if
      mysendlist(id, ishift) = myreceivelist(id, ishift)
    end do
    mysendlist(id, nproc) = -1

  end do
  ! If there are less ikkp then processes or the processes
  ! do not divide nkkp without rest.
  ! Adjust for potential last (or only) partial data injection.
  if(rest > 0) then
    myreceivelist(ndj,:) = -1
    do i = 1, rest
      myreceivelist(ndj, i) = nall*nproc + (rest-i+1)
    end do
    myreceivelist(ndj, :) = cshift(myreceivelist(ndj,:), mypcol1d+1-rest)
    mysendlist(ndj,:) = myreceivelist(ndj, :)
    if(myreadlist(ndj) > 0) then
      myreceivelist(ndj, 1) = 0
    end if
    if(mypcol1d+1 < rest) then
      mysendlist(ndj, nproc) = -1
    end if
    if(mypcol1d == nproc-1) then 
      mysendlist(ndj, nproc-1) = -1
    end if
  end if

  if(verbose) then 
    !! Test output
    call getunit(un)
    fnam = ''
    write(fnam,*) mypcol1d
    fnam = "rese_"//trim(adjustl(fnam))
    open(unit=un, file=trim(fnam), action='write', status='replace')
    write(un,*) "Receive:"
    do i = 1, ndj
      do j = 1, nproc
        write(un, '(I4,1x)', advance='no') myreceivelist(i,j)
      end do
      write(un, *)
    end do
    write(un,*) "Send:"
    do i = 1, ndj
      do j = 1, nproc
        write(un, '(I4,1x)', advance='no') mysendlist(i,j)
      end do
      write(un, *)
    end do
    close(un)

    write(*,*) "Rank",rank,"build lists."
  end if

  !!*************************!!
  !! DATA READ/TRANSFER AND  !!
  !! CONSTRUCTION OF LOCAL H !!
  !!*************************!!
  dataget: do id = 1, ndj
    datashift: do ishift = 1, nproc

      ! Get data this time?
      igetdata: if(myreceivelist(id, ishift) /= -1) then


        !!*************!!
        !! GET DATA    !!
        !!*************!!
        
        ! Receive from file
        if(myreceivelist(id, ishift) == 0) then

          if(verbose) then
            write(*,'(" Rank:", i3, " id:", i4, " ishift:", i3, " RECEIVES DATA FROM FILE")') rank, id, ishift
          end if

          !!***********!!
          !! READ DATA !!
          !!***********!!
          select case(trim(input%xs%bse%bsetype))
            case('singlet', 'triplet')
              ! Read RR part of screened coulomb interaction W_{ouki,o'u'kj}
              call getbsemat('SCCLI.OUT', mysendlist(id, ishift), no, nu, sccli)
          end select
          select case(trim(input%xs%bse%bsetype))
            case('RPA', 'singlet')
              ! Read RR part of exchange interaction v_{ouki,o'u'kj}
              call getbsemat('EXCLI.OUT', mysendlist(id, ishift), no, nu, excli)
          end select
          ! Set ikkp do own list element
          ikkp = mysendlist(id, ishift)

        ! Receive from previous process in ring
        else

          if(verbose) then
            write(*,'(" Rank:",i3, " id:", i4, " ishift:", i3, " RECEIVES DATA FROM RANK", i4)') rank, id, ishift, mysource
          end if

#ifdef SCAL
          !!**************!!
          !! RECEIVE DATA !!
          !!**************!!
          call zgerv2d(ictxt1d, nou, nou, excli(1,1,1,1), nou, 0, mysource)
          call zgerv2d(ictxt1d, nou, nou, sccli(1,1,1,1), nou, 0, mysource)
          ! Set ikkp do sources list element
          ikkp = myreceivelist(id, ishift)
#endif

        end if

        !!**************!!
        !! PROCESS DATA !!
        !!**************!!

        if(verbose) then
          write(*,'(" Rank:",i3, " id:", i4, " ishift:", i3, " PROCESSES DATA FOR IKKP", i4)') rank, id, ishift, ikkp
        end if
        ! Loop over local indices
        do i=1, ham%nrows_loc
          do j=1, ham%ncols_loc
            
#ifdef SCAL
            ! Get corresponding global indices
            a1 = indxl2g(i, mblck, myprow, 0, nprow)
            a2 = indxl2g(j, nblck, mypcol, 0, npcol)
#else
            a1 = i
            a2 = j
#endif
            ! Don't bother with elements corresponding to 
            ! the strictly lower triangular part of the global matrix.
            ! ( We are going to use a hermitian solver )
            if(a2 < a1) then
              cycle
            end if

            ! Current W and V data is for ikkp. Only write
            ! the local entries that correspond to global indices
            ! with that k kp combination.
            ik1_loc = smap(a1,1)
            ik2_loc = smap(a2,1)
            call kkpmap_back(ikkp_loc, nk, ik1_loc, ik2_loc)
            if(ikkp_loc /= ikkp) then 
              cycle
            end if

            !! Construct Hamiltonian elements
            ! For blocks on the diagonal, add the KS transition
            ! energies to the diagonal of the block.
            if(a1 .eq. a2) then
              ham%za(i,j) = kstrans(a1, .true., input%xs%scissor)
            end if
            ! Add exchange term
            ! + 2* v_{o_{s1} u_{s1} k_i, o_{s2} u_{s2} k_j}
            ! * sqrt(abs(f_{o_{s1} k_{s1}} - f_{u_{s1} k_{s1}}))
            ! * sqrt(abs(f_{o_{s2} k_{s2}} - f_{u_{s2} k_{s2}}))
            select case(trim(input%xs%bse%bsetype))
              case('RPA', 'singlet')
                ham%za(i,j) = ham%za(i,j) + exint(a1, a2)
            end select
            ! Add direct term
            ! - W_{o_{s1} u_{s1} k_i, o_{s2} u_{s2} k_j}
            ! * sqrt(abs(f_{o_{s1} k_{s1}} - f_{u_{s1} k_{s1}}))
            ! * sqrt(abs(f_{o_{s2} k_{s2}} - f_{u_{s2} k_{s2}}))
            select case(trim(input%xs%bse%bsetype))
              case('singlet', 'triplet')
                ham%za(i,j) = ham%za(i,j) + scint(a1, a2)
            end select

          end do
        end do

      end if igetdata

#ifdef SCAL
      ! Need to send data further ?
      if(mysendlist(id, ishift) /= -1 .and. nproc /= 1) then
        if(verbose) then
          write(*,'(" Rank:",i3, " id:", i4, " ishift:", i3, " SENDS DATA TO RANK", i4)') rank, id, ishift, mytarget
        end if
        !!***********!!
        !! SEND DATA !!
        !!***********!!
        call zgesd2d(ictxt1d, nou, nou, excli(1,1,1,1), nou, 0, mytarget)
        call zgesd2d(ictxt1d, nou, nou, sccli(1,1,1,1), nou, 0, mytarget)
      end if
#endif
      
    end do datashift
  end do dataget

  ! Work arrays
  deallocate(myreadlist)
  deallocate(myreceivelist)
  deallocate(mysendlist)
  deallocate(excli)
  deallocate(sccli)

  contains

    !BOP
    ! !ROUTINE: kstrans
    ! !INTERFACE:
    complex(8) function kstrans(a1, inbe, scissor)
    ! !USES:
      use modxs, only: bsed
      use mod_eigenvalue_occupancy, only: evalsv
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: a1   ! Combindex index of Hamiltonian ia={ik,iu,io}
    ! real(8) :: scissor ! Scissor correction of band gap
    ! logical :: inbe    ! Subtract gap energy or not
    ! Module in:
    ! integer(4) :: smap(hamsize, 5)  ! Mapping between combined and individual indices
    ! real(8) :: evalsv(nstsv, nkpnr) ! The eigenvalues
    ! real(8) :: egap                 ! The relevant gap of the system 
    ! Out:
    ! complex(8) :: kstrans ! Transition energy
    !
    ! !DESCRIPTION:
    !   Given an combined index $\alpha = \{ \vec{k}_\alpha, o_\alpha, u_\alpha \}$
    !   this routine returns the associated independent particle transition energy
    !   $\Delta\epsilon_\alpha = \epsilon_{\vec{k}_\alpha, u_\alpha} -
    !   \epsilon_{\vec{k}_\alpha, o_\alpha}$. Needs the eigenvalues present in 
    !   {\tt mod\_eigenvalue\_occupancy}.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
      implicit none

      integer(4), intent(in) :: a1
      logical, intent(in) :: inbe
      real(8), intent(in) :: scissor

      integer(4) :: ik, ioabs, iuabs
      real(8) :: deval

      ! Get absolute state band indices form 
      ! combinded index
      ik = smap(a1, 1)
      ioabs = smap(a1, 4)
      iuabs = smap(a1, 5)
      ! Calculate ks energy difference for q=0
      deval = evalsv(iuabs,ik)-evalsv(ioabs,ik)+scissor
      if(inbe) then
        deval = deval - egap - bsed
      end if
      kstrans = cmplx(deval, 0.0d0, 8)
    end function kstrans
    !EOC

    !BOP
    ! !ROUTINE: exint
    ! !INTERFACE:
    complex(8) function exint(a1, a2)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: a1, a2 ! Combined indices of Hamiltonian ia={ik,iu,io}
    ! Module in:
    ! integer(4) :: smap(hamsize, 5)      ! Mapping between combined and individual indices
    ! complex(8) :: excli(no, nu, no, nu) ! Exchange interaction for current k kp combination
    ! real(8) :: ofac(hamsize)            ! Occupation factors
    ! Out:
    ! complex(8) :: exint  ! Exchange part of $H_{\alpha_1,\alpha_2}$
    !
    ! !DESCRIPTION:
    !   Given an two combined indices $\alpha = \{ \vec{k}_\alpha, o_\alpha, u_\alpha \}$
    !   this routine returns the associated exchange interaction part of the BSE-Hamiltonian\\
    !   $\text{exint}_{\alpha_1,\alpha_2} = 2 \cdot \sqrt{\left| f_{\vec{k}_{\alpha_1} o_{\alpha_1}} - f_{\vec{k}_{\alpha_1} u_{\alpha_1}} \right|} 
    !   v_{\vec{k}_{\alpha_1} o_{\alpha_1} u_{\alpha_1}, \vec{k}_{\alpha_2} o_{\alpha_2} u_{\alpha_2}} 
    !   \sqrt{\left| f_{\vec{k}_{\alpha_2} o_{\alpha_2}} - f_{\vec{k}_{\alpha_2} u_{\alpha_2}} \right|}$.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
      integer(4), intent(in) :: a1, a2

      integer(4) :: io1, iu1
      integer(4) :: io2, iu2

      complex(8), parameter :: ztwo=(2.0d0,0.0d0) 

      ! Get individual indices
      io1 = smap(a1, 2) ! Relative
      iu1 = smap(a1, 3) !
      io2 = smap(a2, 2) ! Relative
      iu2 = smap(a2, 3) ! 

      exint = ztwo * excli(io1, iu1, io2, iu2) * ofac(a1) * ofac(a2)

    end function exint
    !EOC

    !BOP
    ! !ROUTINE: scint
    ! !INTERFACE:
    complex(8) function scint(a1, a2)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: a1, a2 ! Combined indices of Hamiltonian ia={ik,iu,io}
    ! Module in:
    ! integer(4) :: smap(hamsize, 5)      ! Mapping between combined and individual indices
    ! complex(8) :: sccli(no, nu, no, nu) ! Screened Coulomb interaction for current k kp combination
    ! real(8) :: ofac(hamsize)            ! Occupation factors
    ! Out:
    ! complex(8) :: scint  ! Screened Coulomb part of $H_{\alpha_1,\alpha_2}$
    !
    ! !DESCRIPTION:
    !   Given an two combined indices $\alpha = \{ \vec{k}_\alpha, o_\alpha, u_\alpha \}$
    !   this routine returns the associated screened Coulomb interaction part of the BSE-Hamiltonian\\
    !   $\text{scint}_{\alpha_1,\alpha_2} = -\sqrt{\left| f_{\vec{k}_{\alpha_1} o_{\alpha_1}} - f_{\vec{k}_{\alpha_1} u_{\alpha_1}} \right|} 
    !   W_{\vec{k}_{\alpha_1} o_{\alpha_1} u_{\alpha_1},\vec{k}_{\alpha_2} o_{\alpha_2} u_{\alpha_2}} 
    !   \sqrt{\left| f_{\vec{k}_{\alpha_2} o_{\alpha_2}} - f_{\vec{k}_{\alpha_2} u_{\alpha_2}} \right|}$.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      integer(4), intent(in) :: a1, a2

      integer(4) :: io1, iu1
      integer(4) :: io2, iu2

      ! Get individual indices
      io1 = smap(a1, 2) ! Relative
      iu1 = smap(a1, 3) !
      io2 = smap(a2, 2) ! Relative
      iu2 = smap(a2, 3) ! 

      scint = -sccli(io1, iu1, io2, iu2) * ofac(a1) * ofac(a2)
    end function scint
    !EOC

end subroutine setup_distributed_bse
!EOC
