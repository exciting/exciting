module m_setup_rmat
  use modmpi
  use modscl
  use modinput
  use m_getunit
  use m_getpmat
  use mod_constants
  use modbse
  use mod_eigenvalue_occupancy, only: evalsv
  use modxs, only: unitout, vkl0
  use m_writecmplxparts

  implicit none

  contains

    !BOP
    ! !ROUTINE: setup_rmat
    ! !INTERFACE:
    subroutine setup_rmat(rmat)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   complex(8) :: rmat(hamsize,3)  ! Position operator matrix 
    ! 
    ! !DESCRIPTION:
    !   The routine generates the position matrix elements 
    !   $R_{\alpha,i} = 
    !    \langle v_\alpha \vec{k}_\alpha | \hat{r}_i | c_\alpha \vec{k}_\alpha \rangle
    !    =
    !    \frac{\langle v_\alpha \vec{k}_\alpha | \hat{p}_i | c_\alpha \vec{k}_\alpha \rangle}
    !    {\epsilon_{c_\alpha, \vec{k}_\alpha}-\epsilon_{v_alpha \vec{k}_\alpha}}$
    !   
    !   Alpha is the combined index used in the BSE Hamiltonian.
    !
    ! !REVISION HISTORY:
    !   Created. (Aurich)
    !EOP
    !BOC
      complex(8), intent(out) :: rmat(hamsize,3)

      complex(8), allocatable :: pmouk(:,:,:,:)
      real(8) :: t1, t0
      integer(4) :: io, iu, ioabs, iuabs, ik, iknr
      integer(4) :: ino, inu, ioabs1, iuabs1, ioabs2, iuabs2 
      integer(4) :: a1

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Building position operator matrix elements using momentum matrix elements  !
      ! and transition energies. If on top of GW, renormalize the p mat elements.  !
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Allocate momentum matrix slice needed.
      allocate(pmouk(3, no_bse_max, nu_bse_max, nk_bse))
      pmouk=zzero

      write(unitout, '("  Reading Pmat.")')
      call timesec(t0)
      if(associated(input%gw)) then
        write(unitout, '("  Also renormalizing Pmat with GW eigenvalues.")')
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
        call getpmat(iknr, vkl0,&
          & ioabs1, ioabs2, iuabs1, iuabs2,&
          & .true., 'PMAT_XS.OUT', pmouk(:,1:ino,1:inu,ik))
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
              pmouk(:,io,iu,ik) = pmouk(:,io,iu,ik)&
                &* (evalsv(iuabs1+iu-1,iknr)-evalsv(ioabs1+io-1,iknr))&
                &/ (eval0(iuabs1+iu-1,iknr)- eval0(ioabs1+io-1,iknr))
            end do
          end do
          !$OMP END PARALLEL DO
        end if 
      end do
      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0

      write(unitout, '("  Building Rmat.")')
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
        
        ! Build R-matrix from P-matrix
        ! \tilde{R}_{a,i} = 
        !   \sqrt{|f_{o_a,k_a}-f_{u_a,k_a}|} *
        !     P^i_{o_a,u_a,k_a} /(e_{u_a, k_a} - e_{o_a, k_a}) 
        ! Note: The scissor does not enter here, so subtract it again.
        rmat(a1, :) = ofac(a1)&
          &*pmouk(:, io, iu, ik)/(de(a1)-sci)

      end do
      !$OMP END PARALLEL DO
      if(rank==0) then 
        if(input%xs%bse%writeparts) then 
          do ik = 1, nk_bse
            call writecmplxparts('Pou_1', remat=dble(pmouk(1,:,:,ik)),&
              & immat=aimag(pmouk(1,:,:,ik)), ik1=ik)
            call writecmplxparts('Pou_2', remat=dble(pmouk(2,:,:,ik)),&
              & immat=aimag(pmouk(2,:,:,ik)), ik1=ik)
            call writecmplxparts('Pou_3', remat=dble(pmouk(3,:,:,ik)),&
              & immat=aimag(pmouk(3,:,:,ik)), ik1=ik)
          end do
          call writecmplxparts('Rmat', remat=dble(rmat), immat=aimag(rmat))
        end if
      end if
      ! Momentum matrix elements no longer needed
      deallocate(pmouk)

      call timesec(t1)
      write(unitout, '("    Time needed",f12.3,"s")') t1-t0
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    end subroutine setup_rmat
    !EOC


    !BOP
    ! !ROUTINE: setup_rmat_dist
    ! !INTERFACE:
    subroutine setup_rmat_dist(rmat, binfo)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    !   type(blacsinfo) :: binfo ! Info type for BLACS grid
    ! In/Out:
    !   type(dzmat) :: rmat  ! 2D block cyclic distributed position matrix 
    ! 
    ! !DESCRIPTION:
    !   The routine sets up the content of the local array of the 
    !   2d block cyclic distributed position operator matrix. Process 0 
    !   reads {\tt PMATXS.OUT} for each {\tt ik}
    !   record and send the data block-wise to the responsible processes.
    !
    ! !REVISION HISTORY:
    !   Created. (Aurich)
    !EOP
    !BOC

      implicit none

      !! I/O
      type(blacsinfo), intent(in) :: binfo
      type(dzmat), intent(inout) :: rmat

      !! Local variables
      integer(4) :: ik, iknr
      integer(4) :: io, iu, ioabs, iuabs
      integer(4) :: ino, inu, ioabs1, iuabs1, ioabs2, iuabs2 

      complex(8), allocatable :: pmou(:,:,:), pmou_t(:,:), pbuff(:,:)

      integer(4) :: context, pr, pc
      integer(4) :: ig, jg, ii, jj, iblck, jblck
      integer(4) :: i, j, m, n, ib, jb, iopt
      integer(4) :: il, jl

      ! Context
      context = rmat%context
      if(context /= binfo%context) then
        write(*,*) "Error (setup_rmat_dist): Context mismatch"
        call terminate
      end if

      ! Block sizes
      iblck = rmat%mblck
      jblck = rmat%nblck 

      ! Allocate work arrays
      if(binfo%myprow == 0 .and. binfo%mypcol == 0) then
        ! Array for reading in the momentum matrix elements
        ! form file for one k point ik.
        ! (maximal size over all ik)
        allocate(pmou(3, no_bse_max, nu_bse_max))
        ! Array for the momentum matrix elements with combined indexing
        ! and only those transitions which have been selected at ik.
        ! (maximal size over all ik)
        allocate(pmou_t(nou_bse_max, 3))
      end if

      ! Send/receive buffer
      allocate(pbuff(iblck, jblck))

      ! Loop over ik blocks of the global 
      ! rmat matrix
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
          call getpmat(iknr, vkl0,&
            & ioabs1, ioabs2, iuabs1, iuabs2,&
            & .true., 'PMAT_XS.OUT', pmou(:,1:ino,1:inu))

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
                pmou(:,io,iu) = pmou(:,io,iu)&
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
              pmou_t(i, iopt) = pmou(iopt, ioabs, iuabs)
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
                rmat%za(il:il+ib-1, jl:jl+jb-1) =&
                  & buildrmat(ib, jb, ofac(ig:ig+ib-1), de(ig:ig+ib-1)-sci,&
                  & pmou_t(i:i+ib-1, j:j+jb-1))

              ! Send data
              else

                ! Prepare send packages
                pbuff(1:ib,1:jb) = pmou_t(i:i+ib-1, j:j+jb-1)
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
              rmat%za(il:il+ib-1, jl:jl+jb-1) =&
                & buildrmat(ib, jb, ofac(ig:ig+ib-1), de(ig:ig+ib-1)-sci, pbuff(1:ib,1:jb))

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

    end subroutine setup_rmat_dist
    !EOC

    !BOP
    ! !ROUTINE: buildrmat
    ! !INTERFACE:
    function buildrmat(ib, jb, occ, ediff, pmat)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: ib, jb        ! Sub block size
    ! real(8)    :: occ(ib)       ! Occupation factors
    ! real(8)    :: ediff(ib)     ! Energy difference e_u-e_o
    ! complex(8) :: pmat(ib, jb)  ! Pou(k) slice
    ! Out:
    ! complex(8) :: buildrmat(ib,jb)  ! Sub block of R-mat
    !
    ! !DESCRIPTION:
    !   The function returns a sub block of the distributed position operator matrix:\\
    !   $R(i_g:i_g+i_b-1, j_g:j_g+j_b-1)$ where each entry is computed according to \\
    !   $R(i, \text{opt}) = F(i) P(i, \text{opt}) / E(i)$  \\
    !   $E$ are the Kohn-Sham transition energies, $F$ the occupation factors and $P$ 
    !   the momentum matrix elements.
    !   The index $i$ correspond to the combined index $\alpha$ and the index 
    !   $\text{opt}$ refers to a Cartesian direction.\\
    !   So with $\alpha = \{ u_\alpha, o_\alpha, \vec{k}_\alpha \}$ : \\
    !   $F_{\alpha} = \sqrt{\left| f_{\vec{k}_{\alpha} o_{\alpha}}
    !                   - f_{\vec{k}_{\alpha} u_{\alpha}} \right|}$, 
    !   $E_{\alpha} = \epsilon_{\vec{k}_{\alpha} u_{\alpha}}
    !                   - \epsilon_{\vec{k}_{\alpha} o_{\alpha}}$, 
    !   $P_{\alpha, \text{opt}} = P^\text{opt}_{o_{\alpha} u_{\alpha} \vec{k}_{\alpha}}$
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
      complex(8) :: buildrmat(ib,jb)

      integer(4), intent(in) :: ib, jb
      real(8), intent(in) :: occ(ib), ediff(ib)
      complex(8), intent(in) :: pmat(ib,jb)
      
      integer(4) :: r, c
      
      !! Construct local 2d block cyclic rmat elements
      ! Build R-matrix from P-matrix 
      ! \tilde{R}_{o_{alpha},u_{alpha},k_{alpha}},i = 
      !   \sqrt{f_{o_{alpha},k_{alpha}}-f_{u_{alpha},k_{alpha}}} *
      !   P_{o_{alpha},u_{alpha},k_{alpha}},i /
      !     (e_{u_{alpha} k_{alpha}} - e_{o_{alpha} k_{alpha}})
      do c = 1, jb
        do r =1, ib
        buildrmat(r, c) = occ(r)*pmat(r, c)/ediff(r)
        end do
      end do

    end function buildrmat
    !EOC

end module m_setup_rmat
