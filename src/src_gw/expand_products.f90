
subroutine expand_products(ik,iq,nstart,nend,nsplit,mstart,mend,msplit,minm)
    use modmain, only : zzero
    use modgw,   only : mbsiz
    implicit none
    
    ! input variables
    integer, intent(in) :: ik
    integer, intent(in) :: iq
    integer, intent(in) :: nstart, nend
    integer, intent(in) :: nsplit         ! last index of valence states
    integer, intent(in) :: mstart ,mend
    integer, intent(in) :: msplit         ! last index of valence states
    complex(8), intent(out) :: minm(mbsiz,nstart:nend,mstart:mend)
    
    ! local variables
    integer :: cstart, cend
    integer :: iflag
    
    minm(:,:,:) = zzero
    
    if ((nsplit>0).and.(msplit<=0)) then

      !======================================================================
      ! When calculating the dielectric function:
      !
      ! n -> {1:nomax}{1:ncg}, so nsplit=nomax, when core states are treated
      ! m -> {numin:nstsv},       msplit is ignored
      !
      ! Important to specify msplit<=0 
      ! to ignore the self energy related part (below)
      ! 
      !======================================================================    
      if (nend<=nsplit) then
    
      ! calculate M^i_{nm}
      iflag = 1
      call expand_products_block(ik,iq,nstart,nend,mstart,mend, &
      &                          minm(:,nstart:nend,mstart:mend), &
      &                          iflag)
      
      else if (nstart>nsplit) then
    
        ! calculate M^i_{cm}
        iflag = 2
        cstart = nstart-nsplit
        cend = nend-nsplit
        call expand_products_block(ik,iq,cstart,cend,mstart,mend, &
        &                          minm(:,nstart:nend,mstart:mend), &
        &                          iflag)
    
      else
    
        ! calculate M^i_{nm}
        iflag = 1
        call expand_products_block(ik,iq,nstart,nsplit,mstart,mend, &
        &                          minm(:,nstart:nsplit,mstart:mend), &
        &                          iflag)
        ! calculate M^i_{cm}
        iflag = 2
        cstart = 1
        cend = nend-nsplit
        call expand_products_block(ik,iq,cstart,cend,mstart,mend, &
        &                          minm(:,nsplit+cstart:nend,mstart:mend), &
        &                          iflag)
      
      end if ! dielectric function part     
    
    else if ((nsplit<=0).and.(msplit>0)) then
      
      !======================================================================
      ! When calculating the self energy:
      !
      ! n -> {ibgw:nbgw},      nsplit is ignored
      ! m -> {1:nstsv}{1:ncg}, so msplit=nstsv, when core states are treated
      !
      ! Important to specify nsplit<=0 
      ! to ignore the dielectric function part (above)
      !
      !======================================================================    
      if (mend<=msplit) then
    
        ! calculate M^i_{nm}
        iflag = 1
        call expand_products_block(ik,iq,nstart,nend,mstart,mend, &
        &                          minm(:,nstart:nend,mstart:mend), &
        &                          iflag)
      
      else if (mstart>msplit) then
    
        ! calculate M^i_{nc}
        iflag = 3
        cstart = mstart-msplit
        cend = mend-msplit
        call expand_products_block(ik,iq,nstart,nend,cstart,cend, &
        &                          minm(:,nstart:nend,mstart:mend), &
        &                          iflag)
    
      else
    
        ! calculate M^i_{nm}
        iflag = 1
        call expand_products_block(ik,iq,nstart,nend,mstart,msplit, &
        &                          minm(:,nstart:nend,mstart:msplit), &
        &                          iflag)
        ! calculate M^i_{nc}
        iflag = 3
        cstart = 1
        cend = mend-msplit
        call expand_products_block(ik,iq,nstart,nend,cstart,cend, &
        &                          minm(:,nstart:nend,msplit+cstart:mend), &
        &                          iflag)
      
      end if ! self energy part
      
    else

      write(*,*) 'emergency stop (expand_product): wrong usage of the subroutine!'
      write(*,*) 'currently implemented options are:'
      write(*,*) "nsplit>0 and msplit<=0 used when calculating the dielectric function"
      write(*,*) "nsplit<=0 and msplit>0 used when calculating the self energy"
      stop
      
    end if
      
    return
    
contains

  !=============================================================================
  !
  !=============================================================================
    subroutine expand_products_block(ik,iq,nstart,nend,mstart,mend,minm,iflag)
        use modmain,    only : zone, zzero
        use modgw,      only : locmatsiz, matsiz, mbsiz, barc, fgw
        use mod_mpi_gw
        implicit none
        ! input variables
        integer(4), intent(in) :: ik    ! the index of the first k-point
        integer(4), intent(in) :: iq    ! the index of the q-point
        integer(4), intent(in) :: nstart, nend  ! range of n states
        integer(4), intent(in) :: mstart, mend  ! range of m states
        complex(8), intent(out):: minm(mbsiz,nstart:nend,mstart:mend)
        integer(4), intent(in) :: iflag ! 1-calcminm; 2-calcmicm; 3-calcminc
        ! local variables
        integer(4) :: ndim, mdim, nmdim
        complex(8), allocatable :: minm_(:,:,:)
#ifdef MPI        
        integer(4) :: mcount, msize, mtype
        integer(4) :: mfirst, mlast
        integer(4) :: mcounts(0:nmax_procs), mdispls(0:nmax_procs) 
        complex(8), allocatable :: minm_p(:,:,:)
#endif

        ndim = nend-nstart+1
        mdim = mend-mstart+1
        
        !if ((nproc_col==1).or.(mdim<nproc_col)) then
        if (.true.) then
          ! write(*,*) '(expand_products::expand_products_block) run in serial mode'
          
          nmdim = ndim*mdim
          select case(iflag)
          
            case(1)
              ! calculate v^{1/2}*M^i_{nm}
              allocate(minm_(matsiz,nstart:nend,mstart:mend))
              call calcminm2(ik,iq,nstart,nend,mstart,mend,minm_)
              call zgemm('c','n', &
              &          mbsiz,nmdim,matsiz, &
              &          zone, &
              &          barc,matsiz, &
              &          minm_,matsiz, &
              &          zzero,minm,mbsiz)
              deallocate(minm_)
              
            case(2)
              ! calculate v^{1/2}*M^i_{cm}
              allocate(minm_(locmatsiz,nstart:nend,mstart:mend))
              call calcmicm(ik,iq,nstart,nend,mstart,mend,minm_)
              call zgemm('c','n', &
              &          mbsiz,nmdim,locmatsiz, &
              &          zone, &
              &          barc,matsiz, &
              &          minm_,locmatsiz, &
              &          zzero,minm,mbsiz)
              deallocate(minm_)
              
            case(3)
              ! calculate v^{1/2}*M^i_{nc}
              allocate(minm_(locmatsiz,nstart:nend,mstart:mend))
              call calcminc(ik,iq,nstart,nend,mstart,mend,minm_)
              call zgemm('c','n', &
              &          mbsiz,nmdim,locmatsiz, &
              &          zone, &
              &          barc,matsiz, &
              &          minm_,locmatsiz, &
              &          zzero,minm,mbsiz)
              deallocate(minm_)
              
            case default
              write(*,*) 'ERROR(expand_products:expand_products_block):'
              write(*,*) 'Unknown iflag=', iflag
              
          end select
        
#ifdef MPI
        else

          ! mtype                : a derived type used for MPI send and receiv
          ! msize                : size of data block to be sent or received
          ! mcount               : size of each pieces in the unit of rstype
          ! mcounts(0:nproc_col) : the number of m-index for each process
          ! mdispls(0:nproc_col) : displacements        
          call mpi_set_range(nproc_col,  &
          &                  myrank_col, &
          &                  mdim,       &
          &                  mstart,     &
          &                  mfirst,     &
          &                  mlast,      &
          &                  mcounts,    &
          &                  mdispls)

          write(*,*) "(expand_products::expand_products_block): Parallel calculations of Minm"
          write(*,*) "myrank, myrank_row, myrank_col:", myrank, myrank_row, myrank_col
          write(*,*) "Range of m-index:", mfirst, mlast

          nmdim = ndim*(mlast-mfirst+1)
          allocate(minm_p(mbsiz,1:ndim,mfirst:mlast))
          select case(iflag)
          
            case(1)
              ! calculate M^i_{nm}
              allocate(minm_(matsiz,nstart:nend,mfirst:mlast))
              call calcminm2(ik,iq,nstart,nend,mfirst,mlast,minm_)
              ! Transform M^i_{nm} to the eigenvectors of the coulomb matrix
              call zgemm('c','n',mbsiz,nmdim,matsiz, &
              &          zone,barc,matsiz,minm_,matsiz, &
              &          zzero,minm_p,mbsiz)
              deallocate(minm_)
              
            case(2)
              ! calculate M^i_{cm}
              allocate(minm_(locmatsiz,nstart:nend,mfirst:mlast))
              call calcmicm(ik,iq,nstart,nend,mfirst,mlast,minm_)
              ! Transform M^i_{cm} to the eigenvectors of the coulomb matrix
              call zgemm('c','n',mbsiz,nmdim,locmatsiz, &
              &          zone,barc,matsiz,minm_,locmatsiz, &
              &          zzero,minm_p,mbsiz)
              deallocate(minm_)
              
            case(3)
              ! calculate M^i_{nc}
              allocate(minm_(locmatsiz,nstart:nend,mfirst:mlast))
              call calcminc(ik,iq,nstart,nend,mfirst,mlast,minm_)
              ! Transform M^i_{nc} to the eigenvectors of the coulomb matrix
              call zgemm('c','n',mbsiz,nmdim,locmatsiz, &
              &          zone,barc,matsiz,minm_,locmatsiz, &
              &          zzero,minm_p,mbsiz)
              deallocate(minm_)
              
            case default
              write(*,*) 'ERROR(expand_products::expand_products_block):'
              write(*,*) 'Unknown iflag=', iflag
              
          end select
          
          write(*,*) "(expand_products::expand_products_block): minm done, collect data"

          msize = mbsiz*ndim
          call MPI_type_contiguous(msize,MPI_DOUBLE_COMPLEX,mtype,ierr)
          call MPI_type_commit(mtype,ierr)

          mcount = mlast-mfirst+1
          call MPI_allgatherv(minm_p,  &
          &                   mcount,  &
          &                   mtype,   &
          &                   minm,    &
          &                   mcounts, &
          &                   mdispls, &
          &                   mtype,   &
          &                   mycomm_col, &
          &                   ierr)
          call MPI_type_free(mtype,ierr)
          deallocate(minm_p)
#endif
          
        end if ! parallel run 
       
      return    
    end subroutine
    
end subroutine
