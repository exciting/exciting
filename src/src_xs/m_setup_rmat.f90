module m_setup_rmat
  use modmpi
  use modscl
  use modinput
  use m_getunit
  use m_getpmat
  use modbse
  use mod_kpoint, only: vkl
  use mod_eigenvalue_occupancy, only: evalsv

  implicit none

  contains

    !BOP
    ! !ROUTINE: setup_distributed_rmat
    ! !INTERFACE:
    subroutine setup_distributed_rmat(rmat, binfo)
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
      integer(4) :: i, j, m, n, ib, jb
      integer(4) :: il, jl

      ! Context
      context = rmat%context
      if(context /= binfo%context) then
        write(*,*) "Error (setup_distributed_rmat): Context mismatch"
        call terminate
      end if

      ! Block sizes
      iblck = rmat%mblck
      jblck = rmat%nblck 

      ! Allocate work arrays
      if(binfo%myprow == 0 .and. binfo%mypcol == 0) then
        allocate(pmou(3, no_bse_max, nu_bse_max))
        allocate(pmou_t(nou_bse_max, 3))
      end if

      ! Send/receive buffer
      allocate(pbuff(iblck, jblck))

      ! Loop over ik blocks of the global 
      ! rmat matrix
      do ik = 1, nk_bse

        iknr = kmap_bse_rg(ik)
        iuabs1 = koulims(1,iknr)
        iuabs2 = koulims(2,iknr)
        ioabs1 = koulims(3,iknr)
        ioabs2 = koulims(4,iknr)
        inu = iuabs2 - iuabs1 + 1
        ino = ioabs2 - ioabs1 + 1

        ! Position of ikkp block in global matrix
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
          call getpmat(iknr, vkl,&
            & ioabs1, ioabs2, iuabs1, iuabs2,&
            & .true., 'PMAT_XS.OUT', pmou(:,1:ino,1:inu))

          ! Din: Renormalise pm according to Del Sole PRB48, 11789(1993)
          ! P^\text{QP}_{okuk} = \frac{E_uk - E_ok}{e_uk - e_ok} P^\text{LDA}_{okuk}
          !   Where E are quasi-particle energies and e are KS energies.
          if(associated(input%gw)) then
            do io = 1, ino
              do iu = 1, inu
                pmou(:,io,iu) = pmou(:,io,iu)&
                  &* (evalsv(iuabs,iknr)-evalsv(ioabs,iknr))/(eval0(iuabs,iknr)&
                  &- eval0(ioabs,iknr))
              end do
            end do
          end if 

          ! Compact array
          do io = 1, ino
            do iu = 1, inu
              pmou_t(iu+(io-1)*ino, :) = pmou(:, io, iu)
            end do
          end do

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
          !! that needs to be sent do one process only.
          ! Get process grid coordinates of responsible process.
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
            ! Get process grid coordinates of responsible process.
            pr = indxg2p( ig, iblck, binfo%myprow, 0, binfo%nprows)
            ! Get position of ib*jb block in local matrix
            il = indxg2l( ig, iblck, pr, 0, binfo%nprows)
#else
            pr=0
            il=ig
#endif

            ! Root does data sending
            if ( binfo%myprow == 0 .and. binfo%mypcol == 0) then 

              ! No send needed, root is responsible for that block
              if( pr == 0 .and. pc == 0) then

                !!********************!!
                !! PROCESS DATA BLOCK !!
                !!********************!!
                ! Assemble sub-block in local R matrix
                rmat%za(il:il+ib-1, jl:jl+jb-1) =&
                  & buildrmat(ig, jg, ib, jb, ofac(ig:ig+ib-1), pmou_t(i:i+ib-1, j:j+jb-1))

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
                & buildrmat(ig, jg, ib, jb, ofac(ig:ig+ib-1), pbuff(1:ib,1:jb))

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

      contains

        !BOP
        ! !ROUTINE: buildrmat
        ! !INTERFACE:
        function buildrmat(ig, jg, ib, jb, occ, pmat)
        ! !INPUT/OUTPUT PARAMETERS:
        ! In:
        ! integer(4) :: ig, jg        ! Position of sub block in global matrix
        ! integer(4) :: ib, jb        ! Sub block size
        ! real(8)    :: occ(ib)       ! Occupation factors
        ! complex(8) :: pmat(ib, jb)  ! $P_{ou}(k)$ slice
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
        !   $P_{\alpha, \text{opt}} = P^\text{opt}_{o_{\alpha} u_{\alpha} \vec{k}_{\alpha}}
        !
        ! !REVISION HISTORY:
        !   Created 2016 (Aurich)
        !EOP
        !BOC
          integer(4), intent(in) :: ig, jg, ib, jb
          real(8), intent(in) :: occ(ib)
          complex(8), intent(inout), optional :: pmat(ib,jb)
          complex(8) :: buildrmat(ib,jb)
          
          real(8) :: evalfactors(ib)
          integer(4) :: iknr, ioabs, iuabs
          integer(4) :: ik, io, iu
          integer(4) :: r, c

          iknr = smap(ig, 3)

          do r = 1, ib
            
            ! State indices could be non
            ! continuous, due to filtering out of
            ! some transition
            iuabs = smap(ig+r-1, 1)
            ioabs = smap(ig+r-1, 2)

            ! Din: Renormalise pm according to Del Sole PRB48, 11789(1993)
            ! P^\text{QP}_{okuk} = \frac{E_uk - E_ok}{e_uk - e_ok} P^\text{LDA}_{okuk}
            !   Where E are quasi-particle energies and e are KS energies.
            if(associated(input%gw)) then
              evalfactors(r) = eval0(iuabs, iknr) - eval0(ioabs, iknr)
            else
              evalfactors(r) = evalsv(iuabs, iknr) - evalsv(ioabs, iknr)
            end if

          end do
          
          !! Construct local 2d block cyclic rmat elements
          ! Build R-matrix from P-matrix 
          ! \tilde{R}_{o_{alpha},u_{alpha},k_{alpha}},i = 
          !   \sqrt{f_{o_{alpha},k_{alpha}}-f_{u_{alpha},k_{alpha}}} *
          !   P_{o_{alpha},u_{alpha},k_{alpha}},i /
          !     (e_{u_{alpha} k_{alpha}} - e_{o_{alpha} k_{alpha}})
          do c = 1, jb
            buildrmat(:, c) = occ(:)*pmat(:, c)/evalfactors(:)
          end do

        end function buildrmat
        !EOC

    end subroutine setup_distributed_rmat
    !EOC

end module m_setup_rmat
