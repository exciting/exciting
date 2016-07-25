
subroutine task_emac_q()

    use modinput
    use modmain
    use modgw
    use mod_mpi_gw
    use m_getunit
      
    implicit none
    integer :: iq, iom, iop, fid, im
    real(8) :: r0, r1, v(3)
    complex(8) :: e0, e1, e2, zt1
    complex(8), allocatable :: emac_q(:,:)
    integer :: i, nq
    real(8) :: a(3), b(3), v1(3), v2(3), v3(3), t1
    integer, allocatable :: iq_path(:)
    real(8), allocatable :: vq_path(:)
#ifdef MPI
    integer :: ntype, ncnt, nsize
#endif
    
    if (rank==0) call boxmsg(fgw,'=','Calculate the q-dependent macroscopic dielectric function for omega=0')

    !==========================
    ! Perform initialization
    !==========================
    
    ! initialize local GW MPI environment
    call init_mpi_gw
    
    ! prepare GW global data
    call init_gw
    
    ! clear not used anymore global exciting variables
    call clean_gndstate
    
    ! occupancy dependent BZ integration weights
    call kintw
    
    !==============
    ! Set \omega=0
    !==============
    freq%nomeg = 1
    iomstart = 1
    iomend = freq%nomeg
    
    !=============================
    ! Setup a path for q-points
    !=============================
    if (.not.associated(input%gw%plot1d)) then
      write(*,*)
      write(*,*) 'ERROR(task_emac_q): Element plot1d is not specified!'
      write(*,*)
      call terminate
    else
      allocate(iq_path(kqset%nkpt))
      iq_path(:) = 0
      ! select a subset of q-points that belong to the specified q-path
      a(:) = input%gw%plot1d%path%pointarray(1)%point%coord(:)
      b(:) = input%gw%plot1d%path%pointarray(2)%point%coord(:)
      nq = 0
      do iq = 1, kqset%nkpt
        v1(:) = kqset%vql(:,iq)-a(:)
        v2(:) = kqset%vql(:,iq)-b(:)
        call r3cross(v1,v2,v3)
        t1  = v3(1)**2+v3(2)**2+v3(3)**2
        if (t1<1d-8) then
          nq = nq+1
          iq_path(nq) = iq
          if (rank==0) then
            write(*,*) 'i, iq_path: ', nq, iq_path(nq)
          end if
        end if
      end do ! iq
    end if
    
    if (nq<1) then
      write(*,*)
      write(*,*) 'ERROR(task_emac_q): No q-points are found along the specified path!'
      write(*,*)
      call terminate
    end if
    
    do i = 1, nq
      iq = iq_path(i)
      Gamma = gammapoint(kqset%vqc(:,iq))
      if (Gamma) exit
    end do
    if ((Gamma).and.(.not.input%gw%rpmat)) then
      !========================================================
      ! calculate momentum matrix elements and store to a file
      !========================================================
      call calcpmatgw
    else
      write(*,*)
      write(*,*) 'Info(task_emac_q): Momentum matrix elements read from files.'
      write(*,*)
    end if

    ! \eps_{00}^{-1}(q)
    allocate(emac_q(nq,iomstart:iomend))
    emac_q(:,1:freq%nomeg) = 0.d0
    
#ifdef MPI    
    iqstart = firstofset(rank,nq)
    iqend = lastofset(rank,nq)
#else
    iqstart = 1
    iqend = nq
#endif
    
    do i = iqstart, iqend
    
      iq = iq_path(i)
      Gamma = gammapoint(kqset%vqc(:,iq))
      
      write(*,*) '--rank ', rank, ' --iq=', iq, ' --Gamma=', Gamma

      !========================================
      ! Calculate interstitial basis functions
      !========================================
      matsiz = locmatsiz+Gqset%ngk(1,iq)
      call diagsgi(iq)
      call calcmpwipw(iq)
    
      !======================================
      ! Calculate the bare coulomb potential
      !======================================
      call calcbarcmb(iq)
    
      !========================================
      ! Set v-diagonal MB and reduce its size
      !========================================
      call setbarcev(input%gw%barecoul%barcevtol)
      call delete_coulomb_potential
    
      !===================================
      ! Calculate the dielectric function
      !===================================
    
      call init_dielectric_function(mbsiz,iomstart,iomend,Gamma)
    
      select case (trim(input%gw%scrcoul%scrtype))
    
        case('ppm','PPM')
          call calcepsilon_ppm(iq,iomstart,iomend)
        
        case default
          call calcepsilon(iq,iomstart,iomend)
          !================================
          ! Invert the dielectric function
          !================================
          call calcinveps(iomstart,iomend)
        
      end select
      
      !---------------------------------
      ! Averaged \epsilon_{00}^{-1}
      !---------------------------------
      if (Gamma) then
        !emac_q(iomstart:iomend,i) = epsh(iomstart:iomend,1,1)+zone
        emac_q(iom,i) = epsilon(1,1,iom)+zone
      else
        do iom = iomstart, iomend
          emac_q(iom,i) = epsilon(1,1,iom)+zone
        end do ! iom
      end if
    
      ! clean unused data
      call delete_dielectric_function(Gamma)
      if (allocated(fnm)) deallocate(fnm)
      if (allocated(mpwipw)) deallocate(mpwipw)
      if (allocated(barc)) deallocate(barc)
      
    end do ! i
    
#ifdef MPI
    call mpi_allgatherv_ifc(nq,freq%nomeg,zbuf=emac_q)
    call barrier
#endif
    
    ! generate q-path
    allocate(vq_path(nq))
    a(:) = input%gw%plot1d%path%pointarray(1)%point%coord(:)
    do i = 1, nq
      v1(:) = kqset%vql(:,iq_path(i))-a(:)
      call r3mv(bvec,v1,v2)
      vq_path(i) = v2(1)**2+v2(2)**2+v2(3)**2
    end do    
    
    if (myrank==0) then
      write(fgw,*)
      write(fgw,*) 'Info(task_emac_q): '
      write(fgw,*) '  Averaged macroscopic dielectric function with and without local field'
      write(fgw,*)
      write(fgw,*)'# frequency    <eps_{00}^{-1}>'
      write(fgw,*)
      call getunit(fid)
      open(fid, File='EPS_Q.OUT', Form='Formatted', Action='Write', Status='Replace')
      write(fid,*)'# q-point    <eps_{00}^{-1}>'
      do i = 1, nq
        write(fid,10) vq_path(i), emac_q(1:freq%nomeg,i)
        write(fgw,10) vq_path(i), emac_q(1:freq%nomeg,i)
      end do
      close(fid)
      10 format(7F12.4)
    end if ! myrank    
    
    deallocate(iq_path)
    deallocate(vq_path)
    deallocate(emac_q)
    
    return
end subroutine
