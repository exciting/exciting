subroutine checkevalqp(fname, nkp2, kvecs2, eqp2)

  use modmpi
  use modinput, only: input
  use constants, only: zzero
  use mod_kpoint, only: nkpt, vkl
  use mod_eigenvalue_occupancy, only: nstsv
  use modgw, only: ibgw, nbgw, nkp1
  use modbse, only: koulims

  implicit none
  character(*), intent(in) :: fname
  integer, intent(in) :: nkp2
  real(8), intent(in) :: kvecs2(3,nkp2)
  real(8), intent(inout):: eqp2(nstsv,nkp2)

  logical       :: exist, fxas
  integer(4)    :: recl
  integer(4)    :: ioglobalmin, iuglobalmax, iuglobalmin

  !-----------------------------------------------------------------------------
  ! Get global band index limits
  !-----------------------------------------------------------------------------      
  iuglobalmax=maxval(koulims(2,:))
  iuglobalmin=minval(koulims(2,:))
  ioglobalmin=minval(koulims(3,:))
  !-----------------------------------------------------------------------------
  ! Find out whether XANES calculation
  !-----------------------------------------------------------------------------      
  fxas=input%xs%bse%xas
  !-----------------------------------------------------------------------------
  ! Read the file
  !-----------------------------------------------------------------------------      
  inquire(File=trim(fname), Exist=exist)
  if (.not.exist) then
    write(*,*)'ERROR(getevalqp): File EVALQP.OUT does not exist!'
    stop
  end if
      
  inquire(IoLength=recl) nkp1, ibgw, nbgw
  open(70, File=trim(fname), Action='READ', Form='UNFORMATTED', &
  &    Access='DIRECT', Recl=recl)
  read(70, Rec=1) nkp1, ibgw, nbgw
  close(70)
  !------------------------------
  ! Data-set consistency check
  !------------------------------
  if (fxas) then
    if ((ibgw > iuglobalmin).or.(nbgw < iuglobalmax)) then
      write(*,*)
      write(*,*)'ERROR(getevalqp):'
      write(*,*)'  Quasiparticle band interval does not agree with BSE band interval!'
      write(*,*)'  Quasiparticle band interval:'
      write(*,*)'  [lowest band, highest band]=[', ibgw, nbgw,']'
      write(*,*)'  XANES-BSE band interval:'
      write(*,*)'  [lowest unoccupied  band, highest unoccupied band]=[', iuglobalmin, iuglobalmax,']'
      write(*,*)
      call terminate
    end if
  else
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
