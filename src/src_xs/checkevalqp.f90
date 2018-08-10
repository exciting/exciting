subroutine checkevalqp(nkp2,kvecs2,eqp2)

  use modmpi
  use modinput
  use mod_constants, only: zzero
  use mod_kpoint, only: nkpt, vkl
  use mod_eigenvalue_occupancy, only: nstsv
  use modgw, only: ibgw, nbgw, nkp1
  use modbse, only: koulims

  implicit none
      
  integer, intent(in) :: nkp2
  real(8), intent(in) :: kvecs2(3,nkp2)
  real(8), intent(inout):: eqp2(nstsv,nkp2)

  logical       :: exist
  integer(4)    :: recl
  integer(4)    :: ioglobalmin, iuglobalmax
  character(30) :: fname
  !-----------------------------------------------------------------------------
  ! Get global band index limits
  !-----------------------------------------------------------------------------      
  iuglobalmax=maxval(koulims(2,:))
  ioglobalmin=minval(koulims(3,:))
  !-----------------------------------------------------------------------------
  ! Read the file
  !-----------------------------------------------------------------------------      

  fname='EVALQP.OUT'
  inquire(File=fname, Exist=exist)
  if (.not.exist) then
    write(*,*)'ERROR(getevalqp): File EVALQP.OUT does not exist!'
    stop
  end if
      
  inquire(IoLength=recl) nkp1, ibgw, nbgw
  open(70, File=fname, Action='READ', Form='UNFORMATTED', &
  &    Access='DIRECT', Recl=recl)
  read(70, Rec=1) nkp1, ibgw, nbgw
  close(70)
  !------------------------------
  ! Data-set consistency check
  !------------------------------
  if ((ibgw > ioglobalmin).or.(nbgw < iuglobalmax)) then
    write(*,*)
    write(*,*)'ERROR(getevalqp):'
    write(*,*)'  Quasiparticle band interval does not agree with BSE band interval!'
    write(*,*)'  Quasiparticle band interval:'
    write(*,*)'  [lowest band, highest band]=[', ibgw, nbgw,']'
    write(*,*)'  BSE band interval:'
    write(*,*)'  [lowest band, highest band]=[', ioglobalmin, iuglobalmax,']'
    write(*,*)
    call terminate
  end if

  !----------------------------------------------
  ! Special case of only one k-point (molecules)
  !----------------------------------------------
  if (nkp1==1) then
    if (nkp2 .ne. 1) then
      write(*,*) 'ERROR(getevalqp):' 
      write(*,*) '  Interpolation is not possible!'
      write(*,*) '  EVALQP.OUT file contains data only for a single k-point.'
      call terminate
    end if
  end if 

  return
end subroutine checkevalqp
