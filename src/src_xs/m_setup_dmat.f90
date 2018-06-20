module m_setup_dmat
  use modmpi
  use modscl
  use modinput
  use m_getunit
  use m_getpmat
  use mod_xsgrids
  use mod_constants
  use modbse
  use mod_eigenvalue_occupancy, only: evalsv
  use modxs, only: unitout, vkl0
  use m_writecmplxparts

  implicit none

  contains

    !BOP
    ! !ROUTINE: setup_dmat
    ! !INTERFACE:
    subroutine setup_dmat(dmat)
    ! !INPUT/OUTPUT PARAMETERS:
    ! Out:
    !   complex(8) :: dmat(hamsize,3)  ! Dipole operator matrix 
    ! 
    ! !DESCRIPTION:
    !   The routine generates the resonant Dipole matrix elements 
    !   $D^\text{rr}_{\alpha,j} = 
    !    \langle u_\alpha \vec{k}_\alpha | -\hat{r}_j | o_\alpha \vec{k}_\alpha \rangle
    !    =
    !    i \frac{\langle u_\alpha \vec{k}_\alpha | \hat{p}_j | o_\alpha \vec{k}_\alpha \rangle}
    !    {\epsilon_{u_\alpha, \vec{k}_\alpha}-\epsilon_{o_alpha \vec{k}_\alpha}}$
    !   
    !   Alpha is the combined index used in the BSE Hamiltonian.
    !
    ! !REVISION HISTORY:
    !   Created. (Aurich)
    !EOP
    !BOC
      complex(8), intent(out) :: dmat(hamsize,3)

      complex(8), allocatable :: pmuok(:,:,:,:), pmou_(:,:,:)
      real(8) :: t1, t0
      integer(4) :: io, iu, ioabs, iuabs, ik, iknr, i
      integer(4) :: ino, inu, ioabs1, iuabs1, ioabs2, iuabs2 
      integer(4) :: a1

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Building position operator matrix elements using momentum matrix elements  !
      ! and transition energies. If on top of GW, renormalize the p mat elements.  !
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Allocate momentum matrix slice needed.
      allocate(pmuok(3, nu_bse_max, no_bse_max, nk_bse))
      pmuok=zzero

      call timesec(t0)
      if(associated(input%gw)) then
        write(unitout, '("  Renormalizing momentum matrix elements with GW eigenvalues.")')
      end if
      ! Read in all needed momentum matrix elements
      do ik = 1, nk_bse
        iknr = kmap_bse_rg(ik)
        iuabs1 = koulims(1,iknr)
        iuabs2 = koulims(2,iknr)
        ioabs1 = koulims(3,iknr)
        ioabs2 = koulims(4,iknr)
        inu = iuabs2 - iuabs1 + 1
        ino = ioabs2 - ioabs1 + 1
        if (input%xs%bse%xas) then
          allocate(pmou_(3,ino,inu))
          call getpmatxas(iknr, vkl0,&
            & ioabs1, ioabs2, iuabs1, iuabs2,&
            & .true., 'PMAT_XS.OUT', pmou_(:,1:ino,1:inu))
          do i=1,3
            pmuok(i,:,:,ik)=transpose(conjg(pmou_(i,:,:)))
          end do
          deallocate(pmou_)
        else
          call getpmat(iknr, vkl0,&
            & iuabs1, iuabs2, ioabs1, ioabs2,&
            & .true., 'PMAT_XS.OUT', pmuok(:,1:inu,1:ino,ik))
        end if
        ! Din: Renormalise pm according to Del Sole PRB48, 11789(1993)
        ! \frac{v^\text{QP}_{okuk}}{E_uk - E_ok} \approx \frac{p^\text{LDA}_{okuk}}{e_uk - e_ok}
        !   In the case that we use the quasi-particle energies E but the LDA eigenfunctions:
        !   1) evalsv contains the E, while eval0 contains the e.
        !   2) The BSE diagonal contains \Delta E
        !   3) But the q->0 limit of <o|exp^{-iqr}|u> is still expressed in terms of
        !      \frac{p^\text{LDA}_{okuk}}{e_uk - e_ok}, instead of v and E.
        !   The following scales the LDA momentum matrix elements so that 3) is 
        !   true for calculations on top of LDA or GW (Eigenvalues only).
        if(associated(input%gw)) then
          !$OMP PARALLEL DO &
          !$OMP& COLLAPSE(2) &
          !$OMP& DEFAULT(SHARED), PRIVATE(io,iu)
          do io = 1, ino
            do iu = 1, inu
              pmuok(:,iu,io,ik) = pmuok(:,iu,io,ik)&
                &* (evalsv(iuabs1+iu-1,iknr) - evalsv(ioabs1+io-1,iknr))&
                &/ (eval0(iuabs1+iu-1,iknr) - eval0(ioabs1+io-1,iknr))
            end do
          end do
          !$OMP END PARALLEL DO
        end if 
      end do
      call timesec(t1)
      if (input%xs%BSE%outputlevelnumber == 1) &
        & write(unitout, '("    Time needed for reading pmat",f12.3,"s")') t1-t0

      call timesec(t0)
      !$OMP PARALLEL DO &
      !$OMP& DEFAULT(SHARED), PRIVATE(a1,iuabs,ioabs,iknr,iu,io,ik)
      do a1 = 1, hamsize

        ! Absolute indices
        iuabs = smap(1,a1)
        ioabs = smap(2,a1)
        iknr  = smap(3,a1)

        ! Relative indices
        iu = smap_rel(1,a1)
        io = smap_rel(2,a1)
        ik = smap_rel(3,a1)
        
        ! Build D-matrix from P-matrix
        ! \tilde{D}_{a,i} = 
        !   \sqrt{|f_{o_a,k_a}-f_{u_a,k_a}|} *
        !     i * P^j_{u_a,o_a,k_a} /(e_{u_a, k_a} - e_{o_a, k_a}) 
        ! Note: The scissor does not enter here, so subtract it again.
        dmat(a1, :) = zi * ofac(a1)&
          &*pmuok(:, iu, io, ik)/(de(a1)-sci)

      end do
      !$OMP END PARALLEL DO

      ! Momentum matrix elements no longer needed
      deallocate(pmuok)

      call timesec(t1)
      if (input%xs%BSE%outputlevelnumber == 1) &
        & write(unitout, '("    Time needed for constructing dmat",f12.3,"s")') t1-t0
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    end subroutine setup_dmat
    !EOC


    !BOP
    ! !ROUTINE: setup_dmat_dist
    ! !INTERFACE:
    subroutine setup_dmat_dist(dmat, binfo)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   type(blacsinfo) :: binfo ! Info type for BLACS grid
    ! In/Out:
    !   type(dzmat) :: dmat  ! 2D block cyclic distributed dipole operator matrix 
    ! 
    ! !DESCRIPTION:
    !   The routine sets up the content of the local array of the 
    !   2d block cyclic distributed dipole operator matrix. Process 0 
    !   reads {\tt PMATXS.OUT} for each {\tt ik}
    !   record and send the data block-wise to the responsible processes.
    !
    !   $D^\text{rr}_{\alpha,j} = 
    !    \langle u_\alpha \vec{k}_\alpha | -\hat{r}_j | o_\alpha \vec{k}_\alpha \rangle
    !    =
    !    i \frac{\langle u_\alpha \vec{k}_\alpha | \hat{p}_j | o_\alpha \vec{k}_\alpha \rangle}
    !    {\epsilon_{u_\alpha, \vec{k}_\alpha}-\epsilon_{o_alpha \vec{k}_\alpha}}$
    !   
    !   Alpha is the combined index used in the BSE Hamiltonian.
    !
    ! !REVISION HISTORY:
    !   Created. (Aurich)
    !EOP
    !BOC

      implicit none

      !! I/O
      type(blacsinfo), intent(in) :: binfo
      type(dzmat), intent(inout) :: dmat

      !! Local variables
      integer(4) :: ik, iknr
      integer(4) :: io, iu, ioabs, iuabs
      integer(4) :: ino, inu, ioabs1, iuabs1, ioabs2, iuabs2 

      complex(8), allocatable :: pmuo(:,:,:), pmuo_t(:,:), pbuff(:,:), pmou_(:,:,:)

      integer(4) :: context, pr, pc
      integer(4) :: ig, jg, ii, jj, iblck, jblck
      integer(4) :: i, j, m, n, ib, jb, iopt
      integer(4) :: il, jl

      if(dmat%isdistributed) then 

        ! Context
        context = dmat%context
        if(context /= binfo%context) then
          write(*,*) "Error (setup_dmat_dist): Context mismatch"
          call terminate
        end if

        ! Block sizes
        iblck = dmat%mblck
        jblck = dmat%nblck 

        ! Allocate work arrays
        if(binfo%myprow == 0 .and. binfo%mypcol == 0) then
          ! Array for reading in the momentum matrix elements
          ! form file for one k point ik.
          ! (maximal size over all ik)
          allocate(pmuo(3, nu_bse_max, no_bse_max))
          ! Array for the momentum matrix elements with combined indexing
          ! and only those transitions which have been selected at ik.
          ! (maximal size over all ik)
          allocate(pmuo_t(nou_bse_max, 3))
        end if

        ! Send/receive buffer
        allocate(pbuff(iblck, jblck))

        ! Loop over ik blocks of the global 
        ! dmat matrix
        ! Loop over participating k-points.
        do ik = 1, nk_bse

          ! Get index of non reduced k points
          ! in the set of all k points.
          iknr = kmap_bse_rg(ik)

          ! Get maximal and minimal unoccupied
          ! and occupied state index participating 
          ! at ik
          iuabs1 = koulims(1,iknr)
          iuabs2 = koulims(2,iknr)
          ioabs1 = koulims(3,iknr)
          ioabs2 = koulims(4,iknr)

          ! Get number of participating unoccupied 
          ! and occupied states (for reading in the
          ! momentum matrix elements)
          inu = iuabs2 - iuabs1 + 1
          ino = ioabs2 - ioabs1 + 1

          ! Position (coordinates of upper left corner)
          ! of ik block in global matrix
          ii = sum(kousize(1:iknr-1)) + 1
          jj = 1

          ! Size of sub-matrix
          m = kousize(iknr)
          n = 3

          !!***********************!!
          !! READ AND PREPARE DATA !! 
          !!***********************!!
          if(binfo%myprow == 0 .and. binfo%mypcol == 0) then

            ! Get momentum P_ou matrix for ik
            ! and those occupied and unoccupied states that 
            ! are participate at that k point. 
            if (input%xs%bse%xas) then
              ! Only pmou is written to file, but pmuo is needed
              allocate(pmou_(3,ino,inu))
              call getpmatxas(iknr, vkl0,&
                &  ioabs1, ioabs2, iuabs1, iuabs2,&
                & .true., 'PMAT_XS.OUT', pmou_(:,1:ino,1:inu))
              do i=1,3
                pmuo(i,:,:)=transpose(conjg(pmou_(i,:,:)))
              end do
              deallocate(pmou_)
            else
              call getpmat(iknr, vkl0,&
                & iuabs1, iuabs2, ioabs1, ioabs2,&
                & .true., 'PMAT_XS.OUT', pmuo(:,1:inu,1:ino))
            end if

            ! Din: Renormalise pm according to Del Sole PRB48, 11789(1993)
            ! \frac{v^\text{QP}_{okuk}}{E_uk - E_ok} \approx \frac{p^\text{LDA}_{okuk}}{e_uk - e_ok}
            !   In the case that we use the quasi-particle energies E but the LDA eigenfunctions:
            !   1) evalsv contains the E, while eval0 contains the e.
            !   2) The BSE diagonal contains \Delta E
            !   3) But the q->0 limit of <o|exp^{-iqr}|u> is still expressed in terms of
            !      \frac{p^\text{LDA}_{okuk}}{e_uk - e_ok}, instead of v and E.
            !   The following scales the LDA momentum matrix elements so that 3) is 
            !   true for calculations on top of LDA or GW (Eigenvalues only).
            if(associated(input%gw)) then
              !$OMP PARALLEL DO &
              !$OMP& COLLAPSE(2) &
              !$OMP& DEFAULT(SHARED), PRIVATE(io,iu)
              do io = 1, ino
                do iu = 1, inu
                  pmuo(:,iu,io) = pmuo(:,iu,io)&
                    &* (evalsv(iuabs1+iu-1,iknr) - evalsv(ioabs1+io-1,iknr))&
                    &/ (eval0(iuabs1+iu-1,iknr) - eval0(ioabs1+io-1,iknr))
                end do
              end do
              !$OMP END PARALLEL DO
            end if 

            ! Compact array and select participating transitions at ik
            !$OMP PARALLEL DO &
            !$OMP& COLLAPSE(2) &
            !$OMP& DEFAULT(SHARED), PRIVATE(i,ioabs,iuabs,iopt)
            do iopt = 1, 3
              do i = 1, m
                iuabs = smap(1, ii+i-1)-iuabs1+1
                ioabs = smap(2, ii+i-1)-ioabs1+1
                ! P_{iopt,io,iu} --> P_t{iou,iopt} iou= u1o1,u2o1,...
                pmuo_t(i, iopt) = pmuo(iopt, iuabs, ioabs)
              end do
            end do
            !$OMP END PARALLEL DO

          end if

          !!******************!!
          !! SEND DATA CHUNKS !!
          !!******************!!
          ! Column index of global ik sub-matrix
          j = 1
          do while(j <= n)

            ! Column index of global matrix
            jg = jj + j - 1

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
#ifdef SCAL
            !! We have a sub block of the 
            !! global sub-matrix at coordinates i,j of size ib*jb
            !! that needs to be sent to one process only.
            ! Get process grid coordinates (column) of responsible process.
            pc = indxg2p( jg, jblck, binfo%mypcol, 0, binfo%npcols)
            ! Get position of ib*jb block in local matrix
            jl = indxg2l( jg, jblck, pc, 0, binfo%npcols)
#else
            pc = 0
            jl = jg
#endif
            ! Row index of global sub-matrix
            i = 1
            do while(i <= m)

              ! Row index of global matrix
              ig = ii + i - 1

              if( i == 1 ) then 
                ! First row block size of global sub-matrix
                ! Adjust for possible truncation of first row block size 
                ib = iblck - mod(ii-1, iblck)
                ib = min(ib, m-i+1)
              else
                ! Adjust for possible truncation of last row block size 
                ib = min(iblck, m-i+1)
              end if
#ifdef SCAL
              !! We have a sub block of the 
              !! global sub-matrix at coordinates i,j of size ib*jb
              !! that needs to be sent do one process only.
              ! Get process grid coordinates (row) of responsible process.
              pr = indxg2p( ig, iblck, binfo%myprow, 0, binfo%nprows)
              ! Get position of ib*jb block in local matrix
              il = indxg2l( ig, iblck, pr, 0, binfo%nprows)
#else
              pr=0
              il=ig
#endif
              ! Root does data sending
              if( binfo%myprow == 0 .and. binfo%mypcol == 0) then 

                ! No send needed, root is responsible for that block
                if( pr == 0 .and. pc == 0) then

                  !!********************!!
                  !! PROCESS DATA BLOCK !!
                  !!********************!!
                  ! Assemble sub-block in local R matrix
                  dmat%za(il:il+ib-1, jl:jl+jb-1) =&
                    & builddmat(ib, jb, ofac(ig:ig+ib-1), de(ig:ig+ib-1)-sci,&
                    & pmuo_t(i:i+ib-1, j:j+jb-1))

                ! Send data
                else

                  ! Prepare send packages
                  pbuff(1:ib,1:jb) = pmuo_t(i:i+ib-1, j:j+jb-1)
#ifdef SCAL
                  ! Send to target process
                  call zgesd2d(context, ib, jb, pbuff, iblck, pr, pc)
#endif

                end if

              ! All others only receive
              else if(binfo%myprow == pr .and. binfo%mypcol == pc) then
#ifdef SCAL
                ! Receive block
                call zgerv2d(context, ib, jb, pbuff, iblck, 0, 0)
#endif
                !!********************!!
                !! PROCESS DATA BLOCK !!
                !!********************!!
                ! Assemble sub-block in local R matrix
                dmat%za(il:il+ib-1, jl:jl+jb-1) =&
                  & builddmat(ib, jb, ofac(ig:ig+ib-1), de(ig:ig+ib-1)-sci, pbuff(1:ib,1:jb))

              end if

              ! Next row block
              i = i + ib

            ! i while loop
            end do

            ! Next column block
            j = j + jb

          ! j while loop
          end do

        ! ik loop
        end do

      else

        call setup_dmat(dmat%za)

      end if

    end subroutine setup_dmat_dist
    !EOC

    !BOP
    ! !ROUTINE: builddmat
    ! !INTERFACE:
    function builddmat(ib, jb, occ, ediff, pmat)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: ib, jb        ! Sub block size
    ! real(8)    :: occ(ib)       ! Occupation factors
    ! real(8)    :: ediff(ib)     ! Energy difference e_u-e_o
    ! complex(8) :: pmat(ib, jb)  ! Puo(k) slice
    ! Out:
    ! complex(8) :: builddmat(ib,jb)  ! Sub block of D-mat
    !
    ! !DESCRIPTION:
    !   The function returns a sub block of the distributed Dipole operator matrix:\\
    !   $D(i_g:i_g+i_b-1, j_g:j_g+j_b-1)$ where each entry is computed according to \\
    !   $D(i, \text{opt}) = F(i) i P(i, \text{opt}) / E(i)$  \\
    !   $E$ are the Kohn-Sham transition energies, $F$ the occupation factors and $P$ 
    !   the momentum matrix elements.
    !   The index $i$ correspond to the combined index $\alpha$ and the index 
    !   $\text{opt}$ refers to a Cartesian direction.\\
    !   So with $\alpha = \{ u_\alpha, o_\alpha, \vec{k}_\alpha \}$ : \\
    !   $F_{\alpha} = \sqrt{\left| f_{\vec{k}_{\alpha} o_{\alpha}}
    !                   - f_{\vec{k}_{\alpha} u_{\alpha}} \right|}$, 
    !   $E_{\alpha} = \epsilon_{\vec{k}_{\alpha} u_{\alpha}}
    !                   - \epsilon_{\vec{k}_{\alpha} o_{\alpha}}$, 
    !   $P_{\alpha, \text{opt}} = P^\text{opt}_{u_{\alpha} o_{\alpha} \vec{k}_{\alpha}}$
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
      complex(8) :: builddmat(ib,jb)

      integer(4), intent(in) :: ib, jb
      real(8), intent(in) :: occ(ib), ediff(ib)
      complex(8), intent(in) :: pmat(ib,jb)
      
      integer(4) :: r, c
      
      !! Construct local 2d block cyclic dmat elements
      ! Build D-matrix from P-matrix 
      ! \tilde{D}_{u_{alpha},o_{alpha},k_{alpha}},j = 
      !   \sqrt{f_{o_{alpha},k_{alpha}}-f_{u_{alpha},k_{alpha}}} *
      !   i*P_{u_{alpha},o_{alpha},k_{alpha}},j /
      !     (e_{u_{alpha} k_{alpha}} - e_{o_{alpha} k_{alpha}})
      do c = 1, jb
        do r =1, ib
        builddmat(r, c) = occ(r)*zi*pmat(r, c)/ediff(r)
        end do
      end do

    end function builddmat
    !EOC

end module m_setup_dmat
