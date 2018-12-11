
subroutine putevalqp()

  use modgw
  use m_getunit

  implicit none

  integer :: recl
  integer :: ik
  integer :: fid

  inquire(IoLength=recl) kset%nkpt, ibgw, nbgw, &
  &       kset%vkl(:,1), &
  &       evalqp(ibgw:nbgw,1), &
  &       evalks(ibgw:nbgw,1), &
  &       eferqp, efermi

  call getunit(fid)

  open(fid, File='EVALQP.OUT', Action='WRITE', Form='UNFORMATTED', &
  &    Access='DIRECT',status='REPLACE', Recl=recl)

  do ik = 1, kset%nkpt

    write(fid, Rec=ik) kset%nkpt, ibgw, nbgw, &
    &     kset%vkl(:,ik), &
    &     evalqp(ibgw:nbgw,ik), &
    &     evalks(ibgw:nbgw,ik), &
    &     eferqp, efermi
  end do ! ikp

  close(fid)

  return
end subroutine
