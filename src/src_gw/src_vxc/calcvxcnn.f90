!BOP
!
!!ROUTINE: calcvxcnn
!
!!INTERFACE:
!
subroutine calcvxcnn
!
!!DESCRIPTION:
!
! This subroutine calculates the diagonal matrix elements of
! the exchange correlation potential (only for valence states).

!!USES:
    use modinput
    use modmain, only: apwordmax, lmmaxapw, lmmaxvr, natmtot, nlomax, &
                       nstfv, nmatmax, nspecies, zzero, &
                       nlotot, natoms, vxcmt, vxcir, zone, &
                       ldapu, llu
    use mod_vxc
    use modgw
    use modmpi
    use m_getunit
    use mod_hdf5
    use mod_hybrids, only: hybridhf, vxnl
    use modxs, only: isreadstate0

!!LOCAL VARIABLES:
    implicit none
    integer(4) :: ikp, ik
    integer(4) :: ist, jst, i, j, k, l, ispn
    integer(4) :: ia, is
    integer(4) :: ngp, fid
    real(8)    :: tstart, tend
    complex(8) :: zsum, zt1
    complex(8), allocatable :: apwalm(:,:,:,:)
    complex(8), allocatable :: evecfv(:,:)
    complex(8), allocatable :: h(:)
    real(8)    :: t0, t1
    character(80) :: filext_save
    logical :: isreadstate0_save

!!EXTERNAL ROUTINES:
    complex(8), external :: zdotc

!!REVISION HISTORY:
!
! Created    8th. Aug. 2006 by RGA
! Revisited       May  2011 by DIN
!
!EOP
!BOC
    call timesec(tstart)

    if (hybridhf) then
      filext_save = filext
      isreadstate0_save = isreadstate0
      filext = '.OUT'
      isreadstate0 = .false.
      call readstate()
      filext = filext_save
      isreadstate0 = isreadstate0_save
      ! read the non-local potential used in hybrids
      call read_vxnl()
    end if

    ! Global array to store <n|Vxc|n>
    if (allocated(vxcnn)) deallocate(vxcnn)
    allocate(vxcnn(nstfv,kset%nkpt))
    vxcnn = zzero

    ! allocate exchange-correlation integral arrays
    if (allocated(vxcraa)) deallocate(vxcraa)
    allocate(vxcraa(apwordmax, &
                    0:input%groundstate%lmaxmat, &
                    apwordmax, &
                    0:input%groundstate%lmaxapw, &
                    0:lmmaxvr, &
                    natmtot))
    if (allocated(vxcrloa)) deallocate(vxcrloa)
    allocate(vxcrloa(nlomax, &
                     apwordmax, &
                     0:input%groundstate%lmaxmat, &
                     0:lmmaxvr, natmtot))
    if (allocated(vxcrlolo)) deallocate(vxcrlolo)
    allocate(vxcrlolo(nlomax, &
                      nlomax, &
                      0:lmmaxvr, &
                      natmtot))

    ! Calculate radial integrals
    call vxcrad

    ! Fourier transform the interstitial part of Vxc
    call genvxcig

    allocate(apwalm(Gkset%ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(evecfv(nmatmax,nstfv))
    allocate(h(nmatmax))

    do ikp = firstofset(rank,kset%nkpt), lastofset(rank,kset%nkpt)

      ik = kset%ikp2ik(ikp)
      ngp = Gkqset%ngk(1,ik)
      call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evecfv)
      call match(ngp, Gkqset%gkc(:,1,ik), Gkqset%tpgkc(:,:,1,ik), &
                 Gkqset%sfacgk(:,:,1,ik), apwalm)

      do i = ibgw, nbgw
        h(:) = zzero
        ! muffin-tin contributions
        do is = 1, nspecies
          do ia = 1, natoms(is)
            call vxcaa(is, ia, ngp, apwalm, evecfv(:,i), h)
            call vxcalo(is, ia, ngp, apwalm, evecfv(:,i), h)
            call vxclolo(is, ia, ngp, evecfv(:,i), h)
            !----------
            ! LDA+U
            !----------
            if ((ldapu /= 0) .and. (llu(is) >= 0)) then
              call vmat_ldapu(is, ia, ngp, apwalm, evecfv(:,i), vxcnn(i,ikp))
            end if
          end do
        end do
        ! interstitial contribution
        call vxcistl(ngp, Gkqset%igkig(:,1,ik), evecfv(:,i), h)
        vxcnn(i,ikp) = vxcnn(i,ikp) + zdotc(ngp+nlotot, evecfv(:,i), 1, h, 1)
      end do ! i

      if (hybridhf) then
        ! setup the hybrid Vxc
        do i = 1, nstfv
          vxcnn(i,ikp) = input%groundstate%Hybrid%excoeff*vxnl(i,i,ikp) + vxcnn(i,ikp)
        end do
      end if

    enddo ! ikp

    deallocate(h)
    deallocate(apwalm)
    deallocate(evecfv)
    deallocate(vxcraa)
    deallocate(vxcrloa)
    deallocate(vxcrlolo)
    if (hybridhf) deallocate(vxnl)

#ifdef MPI
    call mpi_allgatherv_ifc(kset%nkpt, nstfv, zbuf=vxcnn)
    call barrier
#endif

    if (rank==0) then
      ! print results into binary file VXCNN.OUT
      call write_vxcnn()
      ! and text file VXCNN.DAT
      call getunit(fid)
      open(fid,file='VXCNN.DAT',form='FORMATTED',status='UNKNOWN')
      do ikp = 1, kset%nkpt
        write(fid,'("ik=",i4,"    vkl=",3f8.4)') ikp, kset%vkl(:,ikp)
        do i = ibgw, nbgw
          write(fid,'(i4,2f16.6)') i, vxcnn(i,ikp)
        end do
        write(fid,*)
      end do
      close(fid)
    end if

    call timesec(tend)
    time_vxc = time_vxc + (tend-tstart)

end subroutine
!EOC