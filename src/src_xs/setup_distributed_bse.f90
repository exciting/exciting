!BOP
! !ROUTINE: setup_distributed_bse
! !INTERFACE:
subroutine setup_distributed_bse(ham)
! !USES:
  use mod_constants
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
!   2d block cyclic distributed BSE-Hamiltonian matrix. Process 0 
!   reads {\tt EXCLI.OUT} and {\tt SCCLI.OUT} for each {\tt ikkp}
!   record and send the data block-wise to the responsible processes.
!
! !REVISION HISTORY:
!   Created. (Aurich)
!EOP
!BOC

  implicit none

  !! I/O
  type(dzmat), intent(inout) :: ham

  !! Local variables
  integer(4) :: ikkp
  integer(4) :: ik1, ik2

  complex(8), allocatable, dimension(:,:,:,:) :: excli, sccli
  complex(8), allocatable, dimension(:,:) :: excli_t, sccli_t
  logical, allocatable, dimension(:,:) :: lmap

  complex(8), allocatable, dimension(:,:) :: ebuff, sbuff

  integer(4) :: context, pr, pc
  integer(4) :: ig, jg, ii, jj, iblck, jblck
  integer(4) :: i, j, m, n, ib, jb
  integer(4) :: il, jl

  !write(*,*) "setupbse@", myprow, mypcol

  !!*****************************!!
  !! ONLY PROCESS 0 READS DATA   !!
  !! AND SENDS THE CORRESPONDING !!
  !! DATA CHUNKS TO THE OTHERS.  !!
  !!*****************************!!
  ! Allocate read in work arrays
  if(myprow == 0 .and. mypcol == 0) then 
    allocate(excli(nu,no,nu,no)) ! RR part of V
    allocate(sccli(nu,no,nu,no)) ! RR part of W
    allocate(excli_t(nou,nou))
    allocate(sccli_t(nou,nou))
    allocate(lmap(nou,nou))
  end if

  ! Block sizes
  iblck = ham%mblck
  jblck = ham%nblck 

  ! Context
  context = ham%context

  ! Send/receive buffer
  allocate(ebuff(iblck, jblck))
  allocate(sbuff(iblck, jblck))

  !write(6,'("P:",i3," x",i3," allocated work arrays.")') myprow, mypcol
  !write(6,'("P:",i3," x",i3," nkkp=",i8," nk=",i8," no=",i8," nu=",i8)')&
  ! & myprow, mypcol, nkkp, nk, no, nu

  ! Loop over ikkp blocks of the global 
  ! BSE Hamilton matrix.
  do ikkp = 1, nkkp

    ! Get k and kp for current kkp
    call kkpmap(ikkp, nk, ik1, ik2)

    !write(6,'("##### P:",i3," x",i3," ikkp=",i8," ik1=",i8," ik2=",i8)')&
    !  & myprow, mypcol, ikkp, ik1, ik2

    ! Position of ikkp block in global matrix
    ii = sum(kousize(1:ik1-1)) + 1
    jj = sum(kousize(1:ik2-1)) + 1

    !write(6,'("##### P:",i3," x",i3," ii=",i8," jj=",i8)')&
    !  & myprow, mypcol, ii, jj

    ! Size of sub-matrix
    m = kousize(ik1)
    n = kousize(ik2)

    !write(6,'("##### P:",i3," x",i3," m=",i8," n=",i8)')&
    !  & myprow, mypcol, m, n

    !!***********!!
    !! READ DATA !!
    !!***********!!
    if(myprow == 0 .and. mypcol == 0) then

      !write(6,'("***** P:",i3," x",i3," reading file.")') myprow, mypcol

      ! Read in screened coulomb interaction for ikkp
      select case(trim(input%xs%bse%bsetype))
        case('singlet', 'triplet')
          ! Read RR part of screened coulomb interaction W_{uoki,u'o'kj}
          call getbsemat('SCCLI.OUT', ikkp, nu, no, sccli)
      end select

      ! Read in exchange interaction for ikkp
      select case(trim(input%xs%bse%bsetype))
        case('RPA', 'singlet')
          ! Read RR part of exchange interaction v_{uoki,u'o'kj}
          call getbsemat('EXCLI.OUT', ikkp, nu, no, excli)
      end select

      !write(6,'("***** P:",i3," x",i3," Read for ikkp:",i8)') myprow, mypcol, ikkp

      !! Make map for which kou and k'o'u' to use in 
      !! the construction (may me a reduced set due to 
      !! partial occupations)
      lmap = matmul(reshape(kouflag(:,ik1),[nou, 1]), reshape(kouflag(:,ik2),[1, nou]))
      ! Shape interaction arrays more appropriately and select used kou combinations via map
      excli_t = zzero
      sccli_t = zzero
      excli_t(1:m,1:n) = &
        & reshape(pack(reshape(excli,[nou, nou]), lmap), [m, n]) 
      sccli_t(1:m,1:n) = &
        & reshape(pack(reshape(sccli,[nou, nou]), lmap), [m, n])

    end if

    !!******************!!
    !! SEND DATA CHUNKS !!
    !!******************!!
    ! Column index of global ikkp sub-matrix
    j = 1
    do while(j <= n)

      !write(6,'("> P(",i2,") x(",i2,") j=",i8)') myprow, mypcol, j

      ! Column index of global matrix
      jg = jj + j - 1

      !write(6,'("> P(",i2,") x(",i2,") jg=",i8)') myprow, mypcol, jg

      ! Calculate column block size
      if (j == 1) then 
        ! First column block size of global sub-matrix
        ! Adjust for possible truncation of first column block size 
        jb = jblck - mod(jj-1, jblck)
        jb = min(jb, n-j+1)
      else
        ! Adjust for possible truncation of last column block size 
        jb = min(jblck, n-j+1)
      end if

      !write(6,'("> P(",i2,") x(",i2,") jb=",i8)') myprow, mypcol, jb

      ! Row index of global sub-matrix
      i  = 1
      do while(i <= m)
        !write(6,'(">> P(",i2,") x(",i2,") i=",i8)') myprow, mypcol, i

        ! Row index of global matrix
        ig = ii + i - 1
        !write(6,'(">> P(",i2,") x(",i2,") ig=",i8)') myprow, mypcol, ig

        if( i == 1 ) then 
          ! First row block size of global sub-matrix
          ! Adjust for possible truncation of first row block size 
          ib = iblck - mod(ii-1, iblck)
          ib = min(ib, m-i+1)
        else
          ! Adjust for possible truncation of last row block size 
          ib = min(iblck, m-i+1)
        end if

        !! We have a sub block of the 
        !! global sub-matrix at coordinates i,j of size ib*jb
        !! that needs to be sent do one process only.
        ! Get process grid coordinates of responsible process.
        pr = indxg2p( ig, iblck, myprow, 0, nprow)
        pc = indxg2p( jg, jblck, mypcol, 0, npcol)
        ! Get position of ib*jb block in local matrix
        il = indxg2l( ig, iblck, pr, 0, nprow)
        jl = indxg2l( jg, jblck, pc, 0, npcol)

        ! No send needed, root is responsible for that block
        if( pr == 0 .and. pc == 0) then

          !write(6,'(">> Using ",i3," x",i3," subblock of submat", i8,&
          !  &" in root porcess",i3," x", i3)') ib, jb, ikkp, pr, pc

          !!********************!!
          !! PROCESS DATA BLOCK !!
          !!********************!!
          ! Assemble sub-block in local Hamilton matrix
          ham%za(il:il+ib-1, jl:jl+jb-1) =&
            & buildham(ig, jg, ib, jb, ofac(ig:ig+ib-1), ofac(jg:jg+jb-1),&
            & excli_t(i:i+ib-1, j:j+jb-1), sccli_t(i:i+ib-1, j:j+jb-1))

        else

          if( myprow == 0 .and. mypcol == 0) then 

            ! Prepare send packages
            ebuff(1:ib,1:jb) = excli_t(i:i+ib-1, j:j+jb-1)
            sbuff(1:ib,1:jb) = sccli_t(i:i+ib-1, j:j+jb-1)

            ! Send to target process
            !write(6,'("Sending ",i3," x",i3," subblock of submat", i8,&
            !  &" to porcess",i3," x", i3)') ib, jb, ikkp, pr, pc

            call zgesd2d(context, ib, jb, ebuff, iblck, pr, pc)
            call zgesd2d(context, ib, jb, sbuff, iblck, pr, pc)

          else if(myprow == pr .and. mypcol == pc) then

            ! Receive block
            !write(6,'("Receiving ",i3," x",i3," subblock of submat", i8,&
            !  &" form root, process",i3," x",i3)') ib, jb, ikkp, pr, pc
            call zgerv2d(context, ib, jb, ebuff, iblck, 0, 0)
            call zgerv2d(context, ib, jb, sbuff, iblck, 0, 0)

            !!********************!!
            !! PROCESS DATA BLOCK !!
            !!********************!!
            ! Assemble sub-block in local Hamilton matrix
            ham%za(il:il+ib-1, jl:jl+jb-1) =&
              & buildham(ig, jg, ib, jb, ofac(ig:ig+ib-1), ofac(jg:jg+jb-1),&
              & ebuff, sbuff)

          end if

        end if

        ! Next row block
        i = i + ib

      ! i while loop
      end do

      ! Next column block
      j = j + jb

    ! j while loop
    end do

  ! ikkp loop
  end do

  contains

    !BOP
    ! !ROUTINE: buildham
    ! !INTERFACE:
    function buildham(ig, jg, ib, jb, occ1, occ2, exc, scc)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: ig, jg ! Position of sub block in global matrix
    ! integer(4) :: ib, jb ! Sub block size
    ! real(8)    :: occ1(ib), occ2(jb)    ! Occupation factors
    ! complex(8) :: exc(ib, jb)           ! Exchange interaction 
    ! complex(8), optional :: scc(ib, jb) ! Screened Coulomb interaction 
    ! Out:
    ! complex(8) :: buildham(ib,jb)       ! Sub block of BSE-Hamiltonian
    !
    ! !DESCRIPTION:
    !   The function returns a sub block of the distributed BSE-Hamiltonian matrix:\\
    !   $H(i_g:ig+ib-1, j_g:j_g+jb-1)$ where each entry is computed according to \\
    !   $H(i, j) = E(i, j) + F(i) \left( - W(i, j) + 2*V(i, j) \right) F(j)$\\
    !   Only if the sub block contains diagonal elements of the matrix the kohn sham transition
    !   energies $E$ will be added. From the transition energies the gap energy is subtracted and
    !   the scissor is added. The exchange term $V$ is added optionally. \\
    !
    !   The matrix indices correspond to combined indices $\alpha$: \\
    !   Where $\alpha = \{ \vec{k}_\alpha, o_\alpha, u_\alpha \}$, so that: \\
    !   $F_{\alpha} = \sqrt{\left| f_{\vec{k}_{\alpha_1} o_{\alpha_1}} - f_{\vec{k}_{\alpha_1} u_{\alpha_1}} \right|}$ \\
    !   $V_{\alpha_1,\alpha_2} = V_{\vec{k}_{\alpha_1} o_{\alpha_1} u_{\alpha_1}, \vec{k}_{\alpha_2} o_{\alpha_2} u_{\alpha_2}}$ \\
    !   $W_{\alpha_1,\alpha_2} = W_{\vec{k}_{\alpha_1} o_{\alpha_1} u_{\alpha_1}, \vec{k}_{\alpha_2} o_{\alpha_2} u_{\alpha_2}}$ 
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
      integer(4), intent(in) :: ig, jg, ib, jb
      real(8), intent(in) :: occ1(ib), occ2(jb)
      complex(8), intent(in) :: scc(ib,jb)
      complex(8), intent(in), optional :: exc(ib,jb)
      complex(8) :: buildham(ib,jb)
      
      integer(4) :: r, c
      complex(8), parameter :: ztwo=(2.0d0,0.0d0) 
      
      do c = 1, jb
        do r = 1, ib
          if(present(exc)) then
            ! Singlet case with exchange interaction
            buildham(r, c) = occ1(r) * (ztwo * exc(r, c) - scc(r, c)) * occ2(c)
          else
            ! Triplet case without exchange interaction
            buildham(r, c) = -occ1(r) * scc(r, c) * occ2(c)
          end if
          ! Add KS transition energies
          if(ig+r-1 == jg+c-1) then 
            buildham(r, c) = buildham(r, c) + kstrans(ig+r-1)
          end if
        end do
      end do

    end function
    !EOC

    !BOP
    ! !ROUTINE: kstrans
    ! !INTERFACE:
    complex(8) function kstrans(a1)
    ! !USES:
      use modxs, only: bsed
      use modinput, only: input
      use mod_eigenvalue_occupancy, only: evalsv
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: a1   ! Combindex index of Hamiltonian ia={ik,iu,io}
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
    !   \epsilon_{\vec{k}_\alpha, o_\alpha}$. Subtracts the gap energy and 
    !   add a scissor shift.\\
    !   Needs the eigenvalues present in {\tt mod\_eigenvalue\_occupancy}.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
      implicit none

      integer(4), intent(in) :: a1

      integer(4) :: ik, ioabs, iuabs
      real(8) :: deval

      ! Get absolute state band indices form 
      ! combinded index
      iuabs = smap(a1, 1)
      ioabs = smap(a1, 2)
      ik = smap(a1, 3)
      ! Calculate ks energy difference for q=0
      deval = evalsv(iuabs,ik)-evalsv(ioabs,ik)+input%xs%scissor-egap-bsed
      kstrans = cmplx(deval, 0.0d0, 8)
    end function kstrans
    !EOC

end subroutine setup_distributed_bse
!EOC
