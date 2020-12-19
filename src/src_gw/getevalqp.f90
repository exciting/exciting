
subroutine getevalqp(fname, nkp2, kvecs2, eqp2)

  use modmpi
  use modinput
  use constants,                only: zzero
  use mod_kpoint,               only: nkpt, vkl
  use mod_eigenvalue_occupancy, only: nstsv
  use modgw,                    only: ibgw, nbgw, nkp1, kvecs1, eks1, eqp1, eferqp, eferks

  implicit none
  character(*), intent(in)    :: fname
  integer(4),   intent(in)    :: nkp2
  real(8),      intent(in)    :: kvecs2(3,nkp2)
  real(8),      intent(inout) :: eqp2(nstsv,nkp2)

  logical :: exist
  integer(4) :: ik, ib, nb, nk, nqp
  integer(4) :: recl
  real(8), allocatable :: eqp(:)
  complex(8), allocatable :: de1(:,:), de2(:,:)

  !-----------------------------------------------------------------------------
  ! Read the file
  !-----------------------------------------------------------------------------      
  inquire(File=fname, Exist=exist)
  if (.not.exist) then
    write(*,*)'ERROR(getevalqp): File ', trim(fname), ' does not exist!'
    stop
  end if
      
  inquire(IoLength=recl) nkp1, ibgw, nbgw
  open(70, File=trim(fname), Action='READ', Form='UNFORMATTED', &
       Access='DIRECT', Recl=recl)
  read(70, Rec=1) nkp1, ibgw, nbgw
  close(70)
  ! print*, nkp1, ibgw, nbgw
      
  allocate(kvecs1(1:3,nkp1))
  allocate(eqp1(ibgw:nbgw,nkp1))
  allocate(eks1(ibgw:nbgw,nkp1))
  
  inquire(IoLength=recl) nkp1, ibgw, nbgw, kvecs1(1:3,1), &
          eqp1(ibgw:nbgw,1), eks1(ibgw:nbgw,1), &
          eferqp, eferks
  
  open(70, File=trim(fname), Action='READ', Form='UNFORMATTED', &
       Access='DIRECT', Recl=recl)
  
  nqp = nbgw-ibgw+1
  allocate(eqp(nqp))

  do ik = 1, nkp1
    read(70, Rec=ik) nk, ib, nb, kvecs1(:,ik), &
         eqp1(ibgw:nbgw,ik), eks1(ibgw:nbgw,ik), &
         eferqp, eferks
    ! write(*,*) 'ik=', ik
    ! do ib = ibgw, nbgw
    !   write(*,*) ib, eqp1(ib,ik), eks1(ib,ik)
    ! end do
    ! write(*,*)
  end do
  close(70)

  !----------------------------------------------
  ! Special case of only one k-point (molecules)
  !----------------------------------------------
  if (nkp1==1) then
    if (nkp2==1) then
      do ib = ibgw, min(nbgw, nstsv)
        eqp2(ib,1) = eqp1(ib,1)
      end do
      deallocate(kvecs1, eqp1, eks1)
      return
    else
      write(*,*) 'ERROR(getevalqp):' 
      write(*,*) '  Interpolation is not possible!'
      write(*,*) '  EVALQP.OUT file contains data only for a single k-point.'
      stop
    end if
  end if 

  !-----------------------------------------------------------------------------
  ! Interpolate the energies
  !-----------------------------------------------------------------------------      
  allocate(de1(nkp1,ibgw:nbgw))
  do ik = 1, nkp1
    de1(ik,:) = cmplx(eqp1(ibgw:nbgw,ik)-eks1(ibgw:nbgw,ik), 0.d0, 8)
    ! write(*,*) 'ik=', ik
    ! do ib = ibgw, nbgw
    !   write(*,*) ib, de1(ik,ib)
    ! end do
    ! write(*,*)
  enddo

  allocate(de2(nkpt,ibgw:nbgw))
  de2(:,:) = zzero
  
  call fourintp(de1, nkp1, kvecs1, de2, nkp2, vkl, nbgw-ibgw+1)

  do ib = ibgw, min(nbgw,nstsv)
     do ik = 1, nkpt
        eqp2(ib,ik) = eqp2(ib,ik) + dble(de2(ik,ib)) - eferqp
     enddo 
  enddo

  deallocate(de1, de2)
  deallocate(kvecs1, eqp1, eks1)

  return
end subroutine
