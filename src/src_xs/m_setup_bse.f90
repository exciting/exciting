module m_setup_bse
  use modmpi
  use modscl
  use modinput, only: input
  use mod_constants, only: zzero, zone
  use mod_eigenvalue_occupancy, only: evalsv
  use modbse
  use m_getunit
  use m_genfilname
  use m_putgetbsemat
      
  implicit none

  private

  public :: setup_bse, setup_distributed_bse

  contains

    !BOP
    ! !ROUTINE: setup_bse
    ! !INTERFACE:
    subroutine setup_bse(ham, iqmt, fcoup)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   integer(4) :: iqmt  ! Index of momentum transfer Q
    !   logical    :: fcoup ! If true, builds RA instead of RR block of BSE matrix 
    ! In/Out:
    !   complex(8) :: ham(:,:) ! RR or RA block of BSE-Hamiltonian matrix
    ! 
    ! !DESCRIPTION:
    !   The routine sets up the resonant-resonant or resonant-antiresonant block of
    !   the BSE-Hamiltonian matrix. The routine reads {\tt EXCLI.OUT} and
    !   {\tt SCCLI.OUT} ({\tt EXCLIC.OUT} and {\tt SCCLIC.OUT}).
    !
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !EOP
    !BOC

      !! I/O
      complex(8), intent(inout) :: ham(:, :)
      integer(4), intent(in) :: iqmt
      logical, intent(in) :: fcoup

      !! Local variables
      ! Indices
      integer(4) :: ikkp, ik, jk, iknr, jknr
      integer(4) :: ino, inu, jno, jnu, inou, jnou
      integer(4) :: ii, jj
      ! Work arrays
      complex(8), dimension(nou_bse_max, nou_bse_max) :: excli_t, sccli_t

      character(256) :: sfname
      character(256) :: efname

      ! Zero ham
      ham = zzero

      ! Select form which files W and V are to be read
      if(.not. fcoup) then
        call genfilname(basename=scclifbasename, iqmt=iqmt, filnam=sfname)
        call genfilname(basename=exclifbasename, iqmt=iqmt, filnam=efname)
      else
        call genfilname(basename=scclicfbasename, iqmt=iqmt, filnam=sfname)
        call genfilname(basename=exclicfbasename, iqmt=iqmt, filnam=efname)
      end if

      if(rank == 0) then 
        write(unitout, '("setup_bse: Reading form ", a, " and ", a)')&
          & trim(sfname), trim(efname)
      end if

      ! Set up kkp blocks of RR or RA Hamiltonian
      !! Note: If the Hamilton matrix 
      !! has the elements H^RR_{i,j}(qmt) and the indices enumerate the
      !! states according to
      !! i = {u1o1k1, u2o1k1, ..., uM(k1)o1k1,
      !!      uM(k1)o2k1, ..., uM(k1)oN(k1)k1, u1o1k2,
      !!      ..., uM(kO)oN(kO)kO} -> {1,...,\Prod_{i=1}^O M(ki)*N(ki)}
      !! then because of H^RR_{j,i} = H^RR*_{i,j} (or H^RA_{j,i} = H^RA_{i,j})
      !! only jk = ik,..,kO is needed.
      do ikkp = 1, nkkp_bse

        ! Get index ik jk form combination index ikkp
        call kkpmap(ikkp, nk_bse, ik, jk)

        ! Get global index iknr
        iknr = kmap_bse_rg(ik)
        ! Get ranges
        inu = koulims(2,iknr)-koulims(1,iknr)+1
        ino = koulims(4,iknr)-koulims(3,iknr)+1
        inou = inu*ino

        ! Get global index jknr
        jknr = kmap_bse_rg(jk)
        ! Get ranges
        jnu = koulims(2,jknr)-koulims(1,jknr)+1
        jno = koulims(4,jknr)-koulims(3,jknr)+1
        jnou = jnu*jno
        
        ! Read corresponding ikkp blocks of W and V from file
        select case(trim(input%xs%bse%bsetype))
          case('singlet', 'triplet')
            ! Read RR/RA part of screened coulomb interaction W_{iuioik,jujojk}(qmt)
            call b_getbsemat(trim(sfname), iqmt, ikkp, sccli_t(1:inou,1:jnou),&
              & input%xs%bse%blindread)
        end select

        select case(trim(input%xs%bse%bsetype))
          case('RPA', 'singlet')
            ! Read RR/RA part of exchange interaction v_{uoki,u'o'kj}
            call b_getbsemat(trim(efname), iqmt, ikkp, excli_t(1:inou,1:jnou),&
              & input%xs%bse%blindread)
        end select

        ! Position of ikkp block in global matrix
        ii = sum(kousize(1:iknr-1)) + 1
        jj = sum(kousize(1:jknr-1)) + 1

        !! RR only
        ! For blocks on the diagonal, add the KS transition
        ! energies to the diagonal of the block.
        if(.not. fcoup) then 
          if(iknr .eq. jknr) then
            call kstransdiag(ii, iknr, inou, ham(ii:ii+inou-1,jj:jj+jnou-1))
          end if
        end if
          
        !! RR and RA part
        ! Add correlation term and optionally exchange term
        ! (2* v_{iu io ik, jo ju jk}(qmt)
        ! - W_{iu io ik, jo ju jk})(qmt)
        ! * sqrt(abs(f_{io ik} - f_{ju jk+qmt}))
        ! * sqrt(abs(f_{jo jk} - f_{ju jk-qmt}))
        select case(trim(input%xs%bse%bsetype))
          case('RPA', 'singlet')
            call addint(inou, jnou, ham(ii:ii+inou-1,jj:jj+jnou-1),&
              & ofac(ii:ii+inou-1), ofac(jj:jj+jnou-1),&
              & sccli_t(1:inou,1:jnou), exc=excli_t(1:inou,1:jnou))
          case('triplet')
            call addint(inou, jnou, ham(ii:ii+inou-1,jj:jj+jnou-1),&
              & ofac(ii:ii+inou-1), ofac(jj:jj+jnou-1),&
              & sccli_t(1:inou,1:jnou))
        end select

      end do

      contains

        subroutine kstransdiag(ig, ik, m, hamblock)
          integer(4), intent(in) :: ig, ik, m
          complex(8), intent(out) :: hamblock(m,m)

          integer(4) :: io, iu, iou

          ! Calculate ks energy differences
          hamblock = zzero
          do iou = 1, m
            iu = smap(ig+iou-1, 1)
            io = smap(ig+iou-1, 2)
            ! de = e_{u, k+qmt} - e_{o, k} + scissor (RR)
            ! de = e_{u, k-qmt} - e_{o, k} + scissor (AA)
            ! Note: only qmt=0 supported
            hamblock(iou,iou) = cmplx(evalsv(iu,ik) - evalsv(io,ik) + evalshift,8)
          end do
        end subroutine kstransdiag

        subroutine addint(m, n, hamblock, oc1, oc2, scc, exc)
          integer(4), intent(in) :: m, n
          complex(8), intent(inout) :: hamblock(m,n)
          real(8), intent(in) :: oc1(m), oc2(n)
          complex(8), intent(in) :: scc(m,n)
          complex(8), intent(in), optional :: exc(m,n)
          
          integer(4) :: i, j
          integer(4) :: io1, iu1, iou1
          integer(4) :: io2, iu2, iou2
          complex(8) :: tmp

          if(present(exc)) then 
            do j= 1, n
              do i= 1, m
                hamblock(i,j) = hamblock(i,j)&
                  &+ oc1(i)*oc2(j) * (2.0d0 * exc(i,j) - scc(i,j))
              end do
            end do
          else
            do j= 1, n
              do i= 1, m
                hamblock(i,j) = hamblock(i,j) - oc1(i)*oc2(j) * scc(i,j)
              end do
            end do
          end if
        end subroutine addint

    end subroutine setup_bse
    !EOC

    !BOP
    ! !ROUTINE: setup_distributed_bse
    ! !INTERFACE:
    subroutine setup_distributed_bse(ham, iqmt, fcoup, binfo)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   integer(4) :: iqmt ! Index of momentum transfer Q
    !   logical :: fcoup   ! If true, builds RA instead of RR block of BSE matrix 
    !   type(blacsinfo) :: binfo ! Info type of the BLACS grid
    ! In/Out:
    !   type(dzmat) :: ham ! 2D block cyclic distributed RR or RA
    !                      ! block of BSE-Hamiltonian matrix
    ! 
    ! !DESCRIPTION:
    !   The routine sets up the content of the local array of the 
    !   2d block cyclic distributed resonant-resonant or resonant-antiresonant block of
    !   the BSE-Hamiltonian matrix. Process 0 reads {\tt EXCLI.OUT} and {\tt SCCLI.OUT}
    !   ({\tt EXCLIC.OUT} and {\tt SCCLIC.OUT}) for each {\tt ikkp}
    !   record and send the data block-wise to the responsible processes.
    !   If the matrix is not to be distributed the routine calls 
    !   {\tt setup\_bse} instead.
    !
    ! !REVISION HISTORY:
    !   Created. 2016 (Aurich)
    !EOP
    !BOC

      implicit none

      !! I/O
      type(dzmat), intent(inout) :: ham
      integer(4), intent(in) :: iqmt
      logical, intent(in) :: fcoup
      type(blacsinfo), intent(in) :: binfo

      !! Local variables
      integer(4) :: ikkp, ik, jk, iknr, jknr
      integer(4) :: ino, inu, jno, jnu, inou, jnou

      complex(8), allocatable, dimension(:,:) :: excli_t, sccli_t

      complex(8), allocatable, dimension(:,:) :: ebuff, sbuff

      integer(4) :: context, pr, pc
      integer(4) :: ig, jg, ii, jj, iblck, jblck
      integer(4) :: i, j, m, n, ib, jb
      integer(4) :: il, jl
      character(256) :: sfname, efname

      if(ham%isdistributed) then 

        ! Select form which files W and V are to be read
        if(.not. fcoup) then
          call genfilname(basename=scclifbasename, iqmt=iqmt, filnam=sfname)
          call genfilname(basename=exclifbasename, iqmt=iqmt, filnam=efname)
        else
          call genfilname(basename=scclicfbasename, iqmt=iqmt, filnam=sfname)
          call genfilname(basename=exclicfbasename, iqmt=iqmt, filnam=efname)
        end if

        if(mpiglobal%rank == 0) then 
          write(unitout, '("setup_distributed_bse: Reading form ", a, " and ", a)')&
            & trim(sfname), trim(efname)
        end if

#ifdef SCAL
        !!*****************************!!
        !! ONLY PROCESS 0 READS DATA   !!
        !! AND SENDS THE CORRESPONDING !!
        !! DATA CHUNKS TO THE OTHERS.  !!
        !!*****************************!!
        ! Allocate read in work arrays
        if(binfo%myprow == 0 .and. binfo%mypcol == 0) then 
          allocate(excli_t(nou_bse_max, nou_bse_max))
          allocate(sccli_t(nou_bse_max, nou_bse_max))
        end if

        ! Block sizes
        iblck = ham%mblck
        jblck = ham%nblck 

        ! Context
        context = ham%context
        if(context /= binfo%context) then 
          write(*, '("Error (setup_distributed_bse): Context mismatch:", 2i4)')&
            & context, binfo%context
          call terminate
        end if

        ! Send/receive buffer
        allocate(ebuff(iblck, jblck))
        allocate(sbuff(iblck, jblck))

        ! Loop over ikkp blocks of the global 
        ! BSE Hamilton matrix.
        do ikkp = 1, nkkp_bse

          ! Get index ik jk form combination index ikkp
          call kkpmap(ikkp, nk_bse, ik, jk)

          ! Get global index iknr
          iknr = kmap_bse_rg(ik)
          ! Get ranges
          inu = koulims(2,iknr)-koulims(1,iknr)+1
          ino = koulims(4,iknr)-koulims(3,iknr)+1

          ! Get global index jknr
          jknr = kmap_bse_rg(jk)
          ! Get ranges
          jnu = koulims(2,jknr)-koulims(1,jknr)+1
          jno = koulims(4,jknr)-koulims(3,jknr)+1

          ! Position of ikkp block in global matrix
          ii = sum(kousize(1:iknr-1)) + 1
          jj = sum(kousize(1:jknr-1)) + 1

          ! Size of sub-matrix
          inou = inu*ino
          jnou = jnu*jno

          !!***********!!
          !! READ DATA !!
          !!***********!!
          if(binfo%myprow == 0 .and. binfo%mypcol == 0) then

            ! Read in screened coulomb interaction for ikkp
            select case(trim(input%xs%bse%bsetype))
              case('singlet', 'triplet')
                ! Read RR/RA part of screened coulomb interaction W_{iuioik,jujojk}(qmt)
                call b_getbsemat(trim(sfname), iqmt, ikkp, sccli_t(1:inou,1:jnou),&
                  & input%xs%bse%blindread)
            end select

            ! Read in exchange interaction for ikkp
            select case(trim(input%xs%bse%bsetype))
              case('RPA', 'singlet')
                ! Read RR/RA part of exchange interaction v_{iuioik,jujojk}(qmt)
                call b_getbsemat(trim(efname), iqmt, ikkp, excli_t(1:inou,1:jnou),&
                  & input%xs%bse%blindread)
            end select

          end if

          !!******************!!
          !! SEND DATA CHUNKS !!
          !!******************!!
          ! Column index of global ikkp sub-matrix
          j = 1
          do while(j <= jnou)

            ! Column index of global matrix
            jg = jj + j - 1

            ! Calculate column block size
            if (j == 1) then 
              ! First column block size of global sub-matrix
              ! Adjust for possible truncation of first column block size 
              jb = jblck - mod(jj-1, jblck)
              jb = min(jb, jnou-j+1)
            else
              ! Adjust for possible truncation of last column block size 
              jb = min(jblck, jnou-j+1)
            end if

            !! We have a sub block of the 
            !! global sub-matrix at col j of colsize jb
            !! that needs to be sent do one process only.
            ! Get first process grid coordinate of responsible process.
            pc = indxg2p( jg, jblck, binfo%mypcol, 0, binfo%npcols)
            ! Get column position of ib*jb block in local matrix
            jl = indxg2l( jg, jblck, pc, 0, binfo%npcols)

            ! Row index of global sub-matrix
            i  = 1
            do while(i <= inou)

              ! Row index of global matrix
              ig = ii + i - 1

              if( i == 1 ) then 
                ! First row block size of global sub-matrix
                ! Adjust for possible truncation of first row block size 
                ib = iblck - mod(ii-1, iblck)
                ib = min(ib, inou-i+1)
              else
                ! Adjust for possible truncation of last row block size 
                ib = min(iblck, inou-i+1)
              end if

              !! We have a sub block of the 
              !! global sub-matrix at coordinates i,j of size ib*jb
              !! that needs to be sent do one process only.
              ! Get second process grid coordinate of responsible process.
              pr = indxg2p( ig, iblck, binfo%myprow, 0, binfo%nprows)
              ! Get row position of ib*jb block in local matrix
              il = indxg2l( ig, iblck, pr, 0, binfo%nprows)

              ! Root does data sending
              if( binfo%myprow == 0 .and. binfo%mypcol == 0) then 

                ! No send needed, root is responsible for that block
                if( pr == 0 .and. pc == 0) then

                  !!********************!!
                  !! PROCESS DATA BLOCK !!
                  !!********************!!
                  ! Assemble sub-block in local Hamilton matrix
                  call buildham(fcoup, ham%za(il:il+ib-1, jl:jl+jb-1),&
                    & ig, jg, ib, jb,&
                    & occ1=ofac(ig:ig+ib-1), occ2=ofac(jg:jg+jb-1),&
                    & scc=sccli_t(i:i+ib-1, j:j+jb-1), exc=excli_t(i:i+ib-1, j:j+jb-1))

                ! Send data 
                else

                  ! Prepare send packages
                  ebuff(1:ib,1:jb) = excli_t(i:i+ib-1, j:j+jb-1)
                  sbuff(1:ib,1:jb) = sccli_t(i:i+ib-1, j:j+jb-1)
                  call zgesd2d(context, ib, jb, ebuff, iblck, pr, pc)
                  call zgesd2d(context, ib, jb, sbuff, iblck, pr, pc)

                end if

              ! All others only receive
              else if(binfo%myprow == pr .and. binfo%mypcol == pc) then

                ! Receive block
                call zgerv2d(context, ib, jb, ebuff, iblck, 0, 0)
                call zgerv2d(context, ib, jb, sbuff, iblck, 0, 0)

                !!********************!!
                !! PROCESS DATA BLOCK !!
                !!********************!!
                ! Assemble sub-block in local Hamilton matrix
                call buildham(fcoup, ham%za(il:il+ib-1, jl:jl+jb-1),&
                  & ig, jg, ib, jb,&
                  & occ1=ofac(ig:ig+ib-1), occ2=ofac(jg:jg+jb-1),&
                  & scc=sbuff(1:ib,1:jb), exc=ebuff(1:ib,1:jb))

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
#else
        write(*,*) "setup_distributed_bse: Scalapack needed."
        call terminate
#endif

      else

        call setup_bse(ham%za, iqmt, fcoup)

      end if

    end subroutine setup_distributed_bse
    !EOC

    !BOP
    ! !ROUTINE: buildham
    ! !INTERFACE:
    subroutine buildham(fc, hamblck, ig, jg, ib, jb, occ1, occ2, scc, exc)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! logical :: fc ! Build RA instead of RR
    ! integer(4) :: ig, jg ! Position of sub block in global matrix
    ! integer(4) :: ib, jb ! Sub block size
    ! real(8)    :: occ1(ib), occ2(jb)    ! Occupation factors
    ! complex(8), optional :: scc(ib, jb) ! Screened Coulomb interaction 
    ! complex(8), optional :: exc(ib, jb) ! Exchange interaction 
    ! In/Out:
    ! complex(8) :: hamblck(ib,jb)       ! Sub block of BSE-Hamiltonian
    !
    ! !DESCRIPTION:
    !   The routine returns a sub block of the distributed BSE-Hamiltonian matrix:\\
    !   $H(i_g:ig+ib-1, j_g:j_g+jb-1)$ where each entry is computed according to \\
    !   $H(i, j) = E(i, j) + F(i) \left( - W(i, j) + 2*V(i, j) \right) F(j)$\\
    !   Only if the sub block contains diagonal elements of the matrix the kohn sham
    !   transition energies $E$ will be added (RR case only).
    !   From the transition energies the gap energy is subtracted and
    !   the scissor is added. The exchange term $V$ is added optionally. \\
    !
    !   The matrix indices correspond to combined indices $\alpha$: \\
    !   Where $\alpha = \{ \vec{k}_\alpha, o_\alpha, u_\alpha \}$, so that: \\
    !   $F_{\alpha} = \sqrt{ \left| f_{\vec{k}_{\alpha_1} o_{\alpha_1}} 
    !                        - f_{\vec{k}_{\alpha_1} u_{\alpha_1}} \right|}$ \\
    !   $V_{\alpha_1,\alpha_2} = 
    !    V_{\vec{k}_{\alpha_1} o_{\alpha_1} u_{\alpha_1},
    !        \vec{k}_{\alpha_2} o_{\alpha_2} u_{\alpha_2}}$ \\
    !   $W_{\alpha_1,\alpha_2} = 
    !    W_{\vec{k}_{\alpha_1} o_{\alpha_1} u_{\alpha_1},
    !        \vec{k}_{\alpha_2} o_{\alpha_2} u_{\alpha_2}}$ 
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
      integer(4), intent(in) :: ig, jg, ib, jb
      logical, intent(in) :: fc
      real(8), intent(in) :: occ1(ib), occ2(jb)
      complex(8), intent(in), optional :: scc(ib,jb)
      complex(8), intent(in), optional :: exc(ib,jb)
      complex(8), intent(inout) :: hamblck(ib,jb)
      
      integer(4) :: r, c
      complex(8), parameter :: ztwo=(2.0d0,0.0d0) 
      
      do c = 1, jb
        do r = 1, ib
          if(present(exc) .and. present(scc)) then
            ! Singlet case with exchange interaction
            hamblck(r, c) = occ1(r) * (ztwo * exc(r, c) - scc(r, c)) * occ2(c)
          else if(present(scc)) then
            ! Triplet case without exchange interaction
            hamblck(r, c) = -occ1(r) * scc(r, c) * occ2(c)
          end if
          ! Add KS transition energies
          if(.not. fc) then 
            if(ig+r-1 == jg+c-1) then 
              hamblck(r, c) = hamblck(r, c) + kstrans(ig+r-1)
            end if
          end if
        end do
      end do

    end subroutine buildham
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
      ! de = e_{u, k+qmt} - e_{o, k} + scissor (RR)
      ! de = e_{u, k-qmt} - e_{o, k} + scissor (AA)
      ! Note: only qmt=0 supported
      deval = evalsv(iuabs,ik)-evalsv(ioabs,ik)+evalshift
      kstrans = cmplx(deval, 0.0d0, 8)
    end function kstrans
    !EOC

end module m_setup_bse
