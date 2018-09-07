
module mod_vxc

! APW-APW exchange-correlation  integrals
      real(8), allocatable :: vxcraa(:,:,:,:,:,:)
      
! local-orbital-APW exchange-correlation  integrals
      real(8), allocatable :: vxcrloa(:,:,:,:,:)
      
! local-orbital-local-orbital exchange-correlation  integrals
      real(8), allocatable :: vxcrlolo(:,:,:,:)

! G-space interstitial exchange-correlation potential
      complex(8), allocatable :: vxcig(:)
      
! diagonal matrix elements of the exchange-correlation potential
      complex(8), allocatable :: vxcnn(:,:)

contains

!------------------------------------------------------------------
  subroutine write_vxcnn()
    use modgw, only : ibgw, nbgw, kset
    use m_getunit
    implicit none
    integer :: fid, ik, ib

    call getunit(fid)

    ! text format
    open(fid,file='VXCNN.DAT',form='FORMATTED',status='UNKNOWN')
    do ik = 1, kset%nkpt
      write(fid,'(a,3f16.8,4x,f16.8)') '# k-point: ', kset%vkl(:,ik), kset%wkpt(ik)
      do ib = ibgw, nbgw
        write(fid,'(i6,f18.6)') ib, dble(vxcnn(ib,ik))
      end do
      write(fid,*); write(fid,*)
    end do
    close(fid)

    ! fortran binary format
    open(fid, File='VXCNN.OUT', Form='UNFORMATTED', Status='UNKNOWN')
    write(fid) ibgw, nbgw, kset%nkpt, vxcnn(ibgw:nbgw,1:kset%nkpt)
    close(fid)

    return
  end subroutine

!------------------------------------------------------------------
  subroutine read_vxcnn()
    use modgw, only : ibgw, nbgw, kset
    use m_getunit
    implicit none
    integer :: ib, nb, nk
    integer :: fid

    call getunit(fid)
    open(fid,file='VXCNN.OUT',form='UNFORMATTED',status='UNKNOWN')
    read(fid) ib, nb, nk
    close(fid)

    if (nk /= kset%nkpt) then
      write(*,*)'ERROR(mod_vxc::read_vxcnn): Wrong number of k-points'
      write(*,*)'    nk=', nk, '    nkpt=', kset%nkpt
      stop
    end if

    if ((ib /= ibgw).or.(nb /= nbgw)) then
      write(*,*)'ERROR(mod_vxc::read_vxcnn): Different number of bands'
      write(*,*)'    ib=',   ib, '    nb=', nb
      write(*,*)'  ibgw=', ibgw, '  nbgw=', nbgw
      stop
    end if

    open(fid,file='VXCNN.OUT',form='UNFORMATTED',status='UNKNOWN')
    read(fid) ib, nb, nk, vxcnn(ibgw:nbgw,1:kset%nkpt)
    close(fid)
      
    return
  end subroutine
      
end module
