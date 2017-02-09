!
subroutine calcmwm(nstart,nend,mstart,mend)
!
    use modgw,      only: freq, mwm, fgw
    use mod_mpi_gw
    implicit none
    
    ! input variables    
    integer(4), intent(in) :: nstart, nend
    integer(4), intent(in) :: mstart, mend
    ! local variables
    integer(4) :: iom
#ifdef MPI
    integer(4) :: ncount, nsize, ntype
    complex(8), allocatable :: mwm_p(:,:,:)
    complex(8), allocatable :: zbuf(:,:,:)
#endif
    
    if (nproc_col==1) then
    
      ! write(*,*) '(calcmwm): sequential loop over frequencies, rank=', myrank
      do iom = 1, freq%nomeg
        call calcmwm_block(iom,nstart,nend,mstart,mend, &
        &                  mwm(:,mstart:mend,iom))
      end do ! iom
      
#ifdef MPI
    else
    
      ! write(*,*) '(calcmwm): parallel loop over frequencies, rank=', myrank
      allocate(mwm_p(nstart:nend,mstart:mend,iomstart:iomend))
      do iom = iomstart, iomend
        call calcmwm_block(iom,nstart,nend,mstart,mend, &
        &                  mwm_p(:,:,iom))
      end do ! iom
      
      ! collect data from all processes in the same group
      ncount = iomend-iomstart+1
      nsize = (nend-nstart+1)*(mend-mstart+1)
      call MPI_Type_Contiguous(nsize,MPI_DOUBLE_COMPLEX,ntype,ierr)
      call MPI_Type_Commit(ntype,ierr)
      
      allocate(zbuf(nstart:nend,mstart:mend,1:freq%nomeg))
      call MPI_GatherV(mwm_p,ncount,ntype,zbuf,iomcnt,iomdsp,ntype,0,mycomm_col,ierr)
      if (myrank_col==0) then 
        do iom = 1, freq%nomeg
          mwm(:,mstart:mend,iom) = zbuf(:,:,iom) 
        end do 
      end if
      deallocate(zbuf)
      
      call MPI_type_free(ntype,ierr)
      
      deallocate(mwm_p)
   
#endif
    end if ! parallel
    
    return

contains

    !-------------------------------------------------------------------------------
    subroutine calcmwm_block(iom,nstart,nend,mstart,mend,xnm)
        use modinput
        use modmain, only: pi, zone, zzero
        use modgw,   only: vi, kqset, Gamma, singc1, singc2, mbsiz, &
        &                  minmmat, epsilon, epsh, epsw1, epsw2
        implicit none
        ! input variables
        integer(4), intent(in) :: iom
        integer(4), intent(in) :: nstart, nend
        integer(4), intent(in) :: mstart, mend
        complex(8), intent(out):: xnm(nstart:nend,mstart:mend)
        ! local variables
        integer(4) :: ie1, ie2, nmdim
        real(8)    :: wkq
        real(8)    :: vi4pi, coefs1, coefs2
        complex(8) :: wm(mbsiz,nstart:nend,mstart:mend)
        complex(8), external :: zdotu, zdotc
        external zhemm
        
        vi4pi = 4.d0*pi*vi
        coefs1 = singc1*sqrt(vi4pi)
        coefs2 = singc2*vi4pi
        
        ! q-point weight
        wkq = 1.d0/dble(kqset%nkpt)

        ! calculate \sum_{j} W^c_{ij}M^j_{nm}
        nmdim = (nend-nstart+1)*(mend-mstart+1)
        call zhemm('l','u',mbsiz,nmdim, &
        &          zone,epsilon(:,:,iom),mbsiz,minmmat,mbsiz, &
        &          zzero,wm,mbsiz)

        do ie2 = mstart, mend
          do ie1 = nstart, nend
            xnm(ie1,ie2) = wkq*zdotc(mbsiz,minmmat(:,ie1,ie2),1,wm(:,ie1,ie2),1)
            if ((Gamma).and.(ie1==ie2)) then
              xnm(ie1,ie2) = xnm(ie1,ie2) + &
              &  coefs2*epsh(iom,1,1) + &
              &  coefs1*(zdotu(mbsiz,minmmat(:,ie1,ie2),1,epsw2(:,iom,1),1) + &
              &          zdotc(mbsiz,minmmat(:,ie1,ie2),1,epsw1(:,iom,1),1))
            end if ! singular term
          end do
        end do

        return 
    end subroutine

end subroutine
