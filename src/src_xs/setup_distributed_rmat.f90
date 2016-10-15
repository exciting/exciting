!BOP
! !ROUTINE: setup_distributed_rmat
! !INTERFACE:
subroutine setup_distributed_rmat(rmat)
! !USES:
  use modmpi
  use modbse
  use modsclbse
  use modinput
  use m_getunit
  use m_getpmat
  use mod_kpoint, only: vkl
  use mod_eigenvalue_occupancy, only: evalsv
! !INPUT/OUTPUT PARAMETERS:
! In/Out:
! type(dzmat) :: rmat  ! 2D block cyclic distributed position matrix 
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
  type(dzmat), intent(inout) :: rmat

  !! Local variables
  integer(4) :: ik
  integer(4) :: io, iu, iorel, iurel

  complex(8), allocatable, dimension(:,:,:) :: pmou
  complex(8), allocatable, dimension(:,:) :: pmou_t
  logical, allocatable, dimension(:,:) :: lmap

  complex(8), allocatable, dimension(:,:) :: pbuff

  integer(4) :: context, pr, pc
  integer(4) :: ig, jg, ii, jj, iblck, jblck
  integer(4) :: i, j, m, n, ib, jb
  integer(4) :: il, jl

  ! Allocate work arrays
  if(myprow == 0 .and. mypcol == 0) then
    allocate(pmou(3, no, nu))
    allocate(pmou_t(nou, 3))
    allocate(lmap(nou, 3))
  end if

  ! Block sizes
  iblck = rmat%mblck
  jblck = rmat%nblck 

  ! Context
  context = rmat%context

  ! Send/receive buffer
  allocate(pbuff(iblck, jblck))

  ! Loop over ik blocks of the global 
  ! rmat matrix
  do ik = 1, nk

    ! Position of ikkp block in global matrix
    ii = sum(kousize(1:ik-1)) + 1
    jj = 1

    ! Size of sub-matrix
    m = kousize(ik)
    n = 3

    !!***********!!
    !! READ DATA !!
    !!***********!!
    if(myprow == 0 .and. mypcol == 0) then

      ! Get momentum matrix for ik
      call getpmat(ik, vkl,& 
        & bcouabs%il1, bcouabs%iu1,&
        & bcouabs%il2, bcouabs%iu2,&
        & .true., 'PMAT_XS.OUT', pmou)

      ! Make map
      lmap(:,1) = kouflag(:,ik)
      lmap(:,2) = kouflag(:,ik)
      lmap(:,3) = kouflag(:,ik)

      ! Shape momentum matrix arrays more appropriately and select used kou combinations
      pmou_t(1:m,1:n) = &
        & reshape(pack(transpose(reshape(pmou,[3, nou])), lmap), [m, n]) 

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
        jb = min(jb, n-j+1)
      end if

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
          ib = min(ib, m-i+1)
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

          !!********************!!
          !! PROCESS DATA BLOCK !!
          !!********************!!
          ! Assemble sub-block in local R matrix
          rmat%za(il:il+ib-1, jl:jl+jb-1) =&
            & buildrmat(ig, jg, ib, jb, ofac(ig:ig+ib-1), pmou_t(i:i+ib-1, j:j+jb-1))

        else

          if ( myprow == 0 .and. mypcol == 0) then 

            ! Prepare send packages
            pbuff(1:ib,1:jb) = pmou_t(i:i+ib-1, j:j+jb-1)

            ! Send to target process
            call zgesd2d(context, ib, jb, pbuff, iblck, pr, pc)

          else if(myprow == pr .and. mypcol == pc) then

            ! Receive block
            call zgerv2d(context, ib, jb, pbuff, iblck, 0, 0)

            !!********************!!
            !! PROCESS DATA BLOCK !!
            !!********************!!
            ! Assemble sub-block in local R matrix
            rmat%za(il:il+ib-1, jl:jl+jb-1) =&
              & buildrmat(ig, jg, ib, jb, ofac(ig:ig+ib-1), pbuff)

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

  ! ik loop
  end do

  contains

    !BOP
    ! !ROUTINE: buildrmat
    ! !INTERFACE:
    function buildrmat(ig, jg, ib, jb, occ, pmat)
    ! !USES:
      use mod_eigenvalue_occupancy, only: evalsv
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: ig, jg            ! Position of sub block in global matrix
    ! integer(4) :: ib, jb            ! Sub block size
    ! real(8)    :: occ(ib)          ! Occupation factors
    ! complex(8) :: pmat(ib, jb)      ! Exchange interaction 
    ! Out:
    ! complex(8) :: buildrmat(ib,jb)  ! Sub block of R-mat
    !
    ! !DESCRIPTION:
    !   The function returns a sub block of the distributed position operator matrix:\\
    !   $R(i_g:ig+ib-1, j_g:j_g+jb-1)$ where each entry is computed according to \\
    !   $R(i, o) = F(i) * P(i, o) / E(i)$  \\
    !   $E$ are the Kohn-Sham transition energies, $F$ the occupation factors and $P$ the momentum matrix elements.
    !   The index $i$ correspond to the combined indext $\alpha$ and the index $o$ refers to a Cartesian direction.\\
    !   So with $\alpha = \{ \vec{k}_\alpha, o_\alpha, u_\alpha \}$ : \\
    !   $F_{\alpha} = \sqrt{\left| f_{\vec{k}_{\alpha} o_{\alpha}} - f_{\vec{k}_{\alpha} u_{\alpha}} \right|}$ \\
    !   $E_{\alpha} = \epsilon_{\vec{k}_{\alpha} u_{\alpha}} - \epsilon_{\vec{k}_{\alpha} o_{\alpha}} \\
    !   $P_{\alpha, o} = P^o_{\vec{k}_{\alpha} o_{\alpha} u_{\alpha}}
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC
      integer(4), intent(in) :: ig, jg, ib, jb
      real(8), intent(in) :: occ(ib)
      complex(8), intent(inout), optional :: pmat(ib,jb)
      complex(8) :: buildrmat(ib,jb)
      
      integer(4) :: ik, io1, io2, iu1, iu2
      integer(4) :: r, c

      ik = smap(ig, 1)
      io1 = smap(ig, 4)
      iu1 = smap(ig, 5)
      io2 = io1+ib-1
      iu2 = iu1+ib-1

      ! Din: Renormalise pm according to Del Sole PRB48, 11789(1993)
      ! P^\text{QP}_{okuk} = \frac{E_uk - E_ok}{e_uk - e_ok} P^\text{LDA}_{okuk}
      !   Where E are quasi-particle energies and e are KS energies.
      if(associated(input%gw)) then
        do c = 1, jb
          pmat(:,c) = pmat(:, c)&
            &* (evalsv(iu1:iu2, ik) - evalsv(io1:io2, ik))&
            &/ (eval0(iu1:iu2, ik) - eval0(io1:io2, ik))
        end do
      end if 
      
      !! Construct local 2d block cyclic rmat elements
      ! Build complex conjugate R-matrix from p-matrix
      ! \tilde{R}^*_{u_{s1},o_{s1},k_{s1}},i = 
      ! (f_{o_{s1},k_{s1}}-f_{u_{s1},k_{s1}}) *
      !   P_{o_{s1},u_{s1},k_{s1}},i /(e_{o_{s1} k_{s1}} - e_{u_{s1} k_{s1}})
      do c = 1, jb
        buildrmat(:, c) = occ(:)*pmat(:, c)/(evalsv(io1:io2, ik) - evalsv(iu1:iu2, ik))
      end do

    end function buildrmat
    !EOC

end subroutine setup_distributed_rmat
!EOC
