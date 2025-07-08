
subroutine getevalqp(fname, nkp2, kvecs2, eqp2)
  use bandstructure,            only: fourintp
  use constants,                only: zzero
  use mod_bands,                only: nkp1, kvecs1, eks1, eqp1
  use mod_eigenvalue_occupancy, only: nstsv
  use mod_kpoint,               only: nkpt, vkl
  use modgw,                    only: ibgw, nbgw, eferqp, eferks
  use modmpi,                   only: terminate_if_false
  use precision,                only: dp, i32

  implicit none

  character(*), intent(in)    :: fname
  integer(i32), intent(in)    :: nkp2
  real(dp),     intent(in)    :: kvecs2(3, nkp2)
  real(dp),     intent(inout) :: eqp2(nstsv, nkp2)

  logical :: exist
  integer(i32) :: ik, ib, nb, nk, nqp, unit
  integer(i32) :: recl
  real(dp), allocatable :: eqp(:)
  complex(dp), allocatable :: de1(:,:), de2(:,:)

  !-----------------------------------------------------------------------------
  ! Read the file
  !-----------------------------------------------------------------------------      
  inquire(File=fname, Exist=exist)
  call terminate_if_false( exist, 'ERROR(getevalqp): File ' // trim( fname ) // ' does not exist!')
      
  inquire(IoLength=recl) nkp1, ibgw, nbgw
  open(newunit=unit, File=trim(fname), Action='READ', Form='UNFORMATTED', &
       Access='DIRECT', Recl=recl)
  read(unit, Rec=1) nkp1, ibgw, nbgw
  close(unit)
      
  allocate(kvecs1(1:3,nkp1))
  allocate(eqp1(ibgw:nbgw,nkp1))
  allocate(eks1(ibgw:nbgw,nkp1))
  
  inquire(IoLength=recl) nkp1, ibgw, nbgw, kvecs1(1:3,1), &
          eqp1(ibgw:nbgw,1), eks1(ibgw:nbgw,1), &
          eferqp, eferks
  
  open(newunit=unit, File=trim(fname), Action='READ', Form='UNFORMATTED', &
       Access='DIRECT', Recl=recl)
  
  nqp = nbgw-ibgw+1
  allocate(eqp(nqp))

  do ik = 1, nkp1
    read(unit, Rec=ik) nk, ib, nb, kvecs1(:,ik), &
         eqp1(ibgw:nbgw,ik), eks1(ibgw:nbgw,ik), &
         eferqp, eferks
  end do
  close(unit)

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
  enddo

  allocate(de2(nkpt,ibgw:nbgw))
  de2(:,:) = zzero
  
  call fourintp(de1, nkp1, kvecs1, de2, nkp2, vkl, nbgw-ibgw+1)

  do ib = ibgw, min(nbgw,nstsv)
     do ik = 1, nkpt
        eqp2(ib,ik) = eqp2(ib,ik) + dble(de2(ik,ib)) - eferqp
     enddo 
  enddo

  deallocate(kvecs1, eqp1, eks1)

end subroutine
