!BOP
!
!!ROUTINE: calcscgw0
!
!!INTERFACE:
!
subroutine calcscgw0 
!      
! !DESCRIPTION:
!
! This subroutine performs gw calculatons with more radial MPI parallelization 
! including parallelization with respect to both q-points and frequency mesh 
!
!!USES:
    use modinput
    use modmain,    only : zzero, nstsv, evalsv, efermi
    use modgw,      only : kset, kqset, freq,  &
    &                      nomax, numin, nbandsgw, ncg, &
    &                      ibgw, nbgw, time_selfc, time_io, fgw
    use mod_selfenergy, only : evalks, evalqp, eferqp, selfec, mwm  
    use mod_mpi_gw
    use m_getunit
      
!!LOCAL VARIABLES:
    implicit none
    integer(4) :: ikp, iq, mdim
    integer(4) :: isc, nscmax
    integer(4) :: fid
    character(120) :: fname_mwm
    real(8)    :: deltae
    real(8)    :: egk(kset%nkpt), egkold(kset%nkpt)
    real(8)    :: epsilon_sc
    real(8)    :: tstart, tend, t0, t1

    character(len=10), external :: int2str

!!REVISION HISTORY:
!
! Created 16.09.2008 by Hong Jiang
! Readjusted Dec. 2013 by DIN
!      
!EOP
!BOC
    call timesec(tstart)

    ! the loop is done only on root processes in each group  
#ifdef MPI
    if (myrank_col.ne.0) return
    mycomm = mycomm_row
    call mpi_set_range(nproc_row,myrank_row,kqset%nkpt,1,iqstart,iqend)
    write(*,*) '(calcscgw0): run in parallel'
    write(*,*) 'myrank_row, iqstart, iqend = ', myrank_row, iqstart, iqend
#else
    iqstart = 1
    iqend = kqset%nkpt
    write(*,*) '(calcscgw0): run sequentially'
#endif

    ! initialize band gap estimate
    egk(:) = evalqp(numin,:)-evalqp(nomax,:)
    
    ! disable energy shift in calcevalqp.f90
    input%gw%selfenergy%iopes=-1
    nscmax = 100
    epsilon_sc = 1.d-6
    do isc = 1, nscmax
    
      ! update evalsv on the root process
      if (myrank_row==0) then
        evalsv(ibgw:nbgw,:) = evalqp(ibgw:nbgw,:)-eferqp+efermi
      end if
        
#ifdef MPI
      ! Broadcast evalsv to all processes
      if (nproc_row>1) then
        call MPI_BCAST(evalsv,nbandsgw*kset%nkpt, &
        &              MPI_DOUBLE_PRECISION, &
        &              0, mycomm, ierr)
      end if
#endif

      ! combined (core+valence) number of states
      mdim = nstsv
      if (input%gw%coreflag=='all') mdim = mdim+ncg
      
      !====================
      ! sum over q-points
      !====================
      selfec = zzero
      do iq = iqstart, iqend
        
        !===========================================================
        ! Read Sum_ij{M^i*W^c_{ij}(k,q;\omega)*conjg(M^j)} from file
        !===========================================================
        allocate(mwm(ibgw:nbgw,1:mdim,1:freq%nomeg))
        fname_mwm = 'MWM'//'-q'//trim(int2str(iq))//'.OUT'
        call getunit(fid)
        open(fid,File=fname_mwm,Action='Read',Form='Unformatted')
        
        !================================
        ! loop over irreducible k-points
        !================================
        do ikp = 1, kset%nkpt
        
          ! read M*W*M matrix
          call timesec(t0)
          read(fid) mwm
          call timesec(t1)
          time_io = time_io+t1-t0
          
          !===========================================================
          ! Calculate the contribution to the correlation self-energy
          !===========================================================
          call calcselfc_freqconv(ikp,iq,mdim)
          
        end do ! ikp
        
        deallocate(mwm)
        close(fid)
        
      end do ! iq

#ifdef MPI
      if (nproc_row > 1) then
        call mpi_sum_array(0,selfec,nbandsgw,freq%nomeg,kset%nkpt,mycomm)
      end if
#endif

      if (myrank_row==0) then
      
        call calcevalqp

        egkold(:) = egk(:)
        egk(:) = evalqp(numin,:)-evalqp(nomax,:)
        
        deltae = abs(maxval(egk-egkold))
        write(fgw,'(a,i5,e12.4)') "(calcscgw0): isc, deltae", isc, deltae
        
      end if ! myrank

#ifdef MPI
      if (nproc_row>1) then 
        call MPI_BCAST(deltae,1,MPI_DOUBLE_PRECISION,0,mycomm,ierr)
      end if 
#endif

      if (deltae < epsilon_sc) exit
      
    enddo  !! loop over isc

    if (myrank_row==0) then
      if (isc > nscmax) then
        write(fgw,*) '(calcscgw0): WARNING convergence is not reached!'
      end if 
    end if
    
    ! timing
    call timesec(tend)
    time_selfc = time_selfc+tend-tstart
    
    return
end subroutine
!EOC      
