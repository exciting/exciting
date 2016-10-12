!BOP
! !ROUTINE: setup_distributed_rmat
! !INTERFACE:
subroutine setup_distributed_rmat(rmat)
! !USES:
  use modmpi
  use modbse
  use modsclbse
  use mod_rrsring
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
!   2d block cyclic distributed position operator matrix.
!
! !REVISION HISTORY:
!   Created. (Aurich)
!EOP
!BOC
  implicit none

  !! I/O
  type(dzmat), intent(inout) :: rmat

  !! Local variables
  integer(4) :: ik, ik_loc, io, iu, iorel, iurel
  integer(4) :: id, iprc, ishift, a1, opt, i, j
  complex(8), allocatable, dimension(:,:,:) :: pmou

  integer(4) :: un
  logical :: verbose
  character(256) :: fnam
  verbose = .true.

  ! Allocate work arrays
  allocate(pmou(3, no, nu))

  ! Setup lists of k points to be 
  ! read, recived and sent by each process.
  write(*,*) "nk=", nk
  call setup_rrslists(nk)

  !! Test output
  if(verbose) then
    if(rank == 0) then
      write(*,'(" Number of data injections: ", i3, 1x, i3, 1x, i3)'), ndj, nall, rest
    end if
    write(*,'(" Rank", i3," Rank1d", i3 " Form rank:", i3, " To rank:", i3)') rank,&
      & mypcol1d_r, mysource, mytarget

    call getunit(un)
    fnam = ''
    write(fnam,*) mypcol1d_r
    fnam = "rese_rmat_"//trim(adjustl(fnam))
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
            write(*,'(" Rank:", i3, " id:", i4, " ishift:", i3,&
              & " RECEIVES DATA FROM FILE")') rank, id, ishift
          end if

          !!***********!!
          !! READ DATA !!
          !!***********!!
          ! Set ik do own list element
          ik = mysendlist(id, ishift)
          call getpmat(ik, vkl,& 
            & bcouabs%il1, bcouabs%iu1,&
            & bcouabs%il2, bcouabs%iu2,&
            & .true., 'PMAT_XS.OUT', pmou)

        ! Receive from previous process in ring
        else

          if(verbose) then
            write(*,'(" Rank:",i3, " id:", i4, " ishift:", i3,&
              &" RECEIVES DATA FROM RANK", i4)') rank, id, ishift, mysource
          end if

#ifdef SCAL
          !!**************!!
          !! RECEIVE DATA !!
          !!**************!!
          call zgerv2d(ictxt1d_r, 3*no, nu, pmou(1,1,1), 3*no, 0, mysource)
          ! Set ik do sources list element
          ik = myreceivelist(id, ishift)
#endif

        end if

        !!**************!!
        !! PROCESS DATA !!
        !!**************!!

        if(verbose) then
          write(*,'(" Rank:",i3, " id:", i4, " ishift:", i3,&
            & " PROCESSES DATA FOR IK", i4)') rank, id, ishift, ik
        end if
        ! Loop over local indices
        do i=1, rmat%nrows_loc
          do j=1, rmat%ncols_loc
            
#ifdef SCAL
            ! Get corresponding global indices
            a1  = indxl2g(i, rmat%mblck, myprow1d_c, 0, nprow1d_c)
            opt = indxl2g(j, rmat%mblck, mypcol1d_c, 0, npcol1d_c)
#else
            a1 = i
            opt = j
#endif
            ! Current PMAT data is for ik. Only write
            ! the local entries that correspond to global indices
            ! with that ik.
            ik_loc = smap(a1,1)
            if(ik_loc /= ik) then 
              cycle
            end if
            ! Get state indices w.r.t. the read PMAT block
            iorel = smap(a1, 2)
            iurel = smap(a1, 3)
            ! Get absolute state indices
            io = smap(a1, 4)
            iu = smap(a1, 5)

            ! Din: Renormalise pm according to Del Sole PRB48, 11789(1993)
            ! P^\text{QP}_{okuk} = \frac{E_uk - E_ok}{e_uk - e_ok} P^\text{LDA}_{okuk}
            !   Where E are quasi-particle energies and e are KS energies.
            if(associated(input%gw)) then
               pmou(opt,iorel,iurel) = pmou(opt,iorel,iurel)&
                 &* (evalsv(iu,ik)-evalsv(io,ik))/(eval0(iu,ik)- eval0(io,ik))
            end if 

            !! Construct local 2d block cyclic rmat elements
            ! Build complex conjugate R-matrix from p-matrix
            ! \tilde{R}^*_{u_{s1},o_{s1},k_{s1}},i = 
            ! (f_{o_{s1},k_{s1}}-f_{u_{s1},k_{s1}}) *
            !   P_{o_{s1},u_{s1},k_{s1}},i /(e_{o_{s1} k_{s1}} - e_{u_{s1} k_{s1}})
            rmat%za(i, j) = ofac(a1)&
              &* pmou(opt, iorel, iurel)/(evalsv(io, ik) - evalsv(iu, ik))

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
        call zgesd2d(ictxt1d_r, 3*no, nu, pmou(1,1,1), 3*no, 0, mytarget)
      end if
#endif
      
    end do datashift
  end do dataget

  ! Work arrays
  deallocate(myreceivelist)
  deallocate(mysendlist)
  deallocate(pmou)

end subroutine setup_distributed_rmat
!EOC
