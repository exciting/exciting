subroutine eqpongrid
  Use modinput
  Use modmain
  use mod_wannier

  implicit none

  integer :: nv, ik, nk, ib, nb, iknr, is, iv, m, ix, iy, iz
  integer :: nkp, ibgw, nbgw
  real(8) :: s(3), v1(3), dt, vk(3)
  real(8) :: eferqp, eferks
  real(8), allocatable :: eval2(:,:), vvl(:,:), ongrid(:,:)
  real(8), allocatable :: kvecs(:,:), eqp(:,:), eks(:,:)
  logical :: exist
  character(30) :: fname
  integer(4)    :: recl

  !-----------------------------------------------------------------------------
  ! Read the file
  !-----------------------------------------------------------------------------      

  fname='EVALQP.OUT'
  inquire(File=fname, Exist=exist)
  if (.not.exist) then
    write(*,*)'ERROR(getevalqp): File EVALQP.OUT does not exist!'
    stop
  end if
      
  inquire(IoLength=recl) nkp, ibgw, nbgw
  open(70, File=fname, Action='READ', Form='UNFORMATTED', &
  &    Access='DIRECT', Recl=recl)
  read(70, Rec=1) nkp, ibgw, nbgw
  close(70)
      
  allocate(kvecs(1:3,nkp))
  allocate(eqp(ibgw:nbgw,nkp))
  allocate(eks(ibgw:nbgw,nkp))
  
  inquire(IoLength=recl) nkp, ibgw, nbgw, kvecs(1:3,1), &
  &       eqp(ibgw:nbgw,1), eks(ibgw:nbgw,1), &
  &       eferqp, eferks
  
  open(70, File=fname, Action='READ', Form='UNFORMATTED', &
  &    Access='DIRECT', Recl=recl)
  
  do ik = 1, nkp
    read(70, Rec=ik) nk, ib, nb, kvecs(:,ik), &
    &    eqp(ibgw:nbgw,ik), eks(ibgw:nbgw,ik), &
    &    eferqp, eferks
  end do
  close(70)
  do ik = 1, nkp
    call findkpt( kvecs( :, ik), is, iv)
    eks( ibgw:nbgw, iv) = eqp( ibgw:nbgw, ik)
  end do

  nv = size( input%properties%bandstructure%plot1d%path%pointarray)
  allocate( vvl( 3, nv))
  
  do ik = 1, nv
    vvl( :, ik) = input%properties%bandstructure%plot1d%path%pointarray( ik)%point%coord
  end do
  
  open( 50, file='BANDONGRID-QP.OUT', action='WRITE', form='FORMATTED')
  dt = 0.d0
  do ik = 1, nv-1
    v1 = vvl( :, ik+1) - vvl( :, ik)
    do iz = -wf_kset%ngridk(1), wf_kset%ngridk(1)
      do iy = -wf_kset%ngridk(1), wf_kset%ngridk(1)
        do ix = -wf_kset%ngridk(1), wf_kset%ngridk(1)
          vk(1) = dble( ix)/wf_kset%ngridk(1)
          vk(2) = dble( iy)/wf_kset%ngridk(2)
          vk(3) = dble( iz)/wf_kset%ngridk(3)
          s = -1.d0
          if( abs( v1( 1)) .ge. input%structure%epslat) s(1) = (vk(1) - vvl( 1, ik))/v1( 1)
          if( abs( v1( 2)) .ge. input%structure%epslat) s(2) = (vk(2) - vvl( 2, ik))/v1( 2)
          if( abs( v1( 3)) .ge. input%structure%epslat) s(3) = (vk(3) - vvl( 3, ik))/v1( 3)
          !write( *, '(3F13.6)') s
          if( (s(2) .ge. 0.d0) .and. (s(2) .le. 1.d0)) s(1) = s(2)
          if( (s(3) .ge. 0.d0) .and. (s(3) .le. 1.d0)) s(1) = s(3)
          if( (s(1) .ge. 0.d0) .and. (s(1) .le. 1.d0) .and. &
              (norm2( vvl( :, ik) + s(1)*v1 - vk) .lt. input%structure%epslat)) then
            call r3mv( bvec, vk - vvl( :, ik), s)
            call findkpt( vk, is, iv)
            do is = wf_fst, wf_lst
              write( 50, '(8G18.10,I)') dt + norm2( s), eks( is, iv)-eferqp, vk, vkl( :, iv), iv
            end do
            write( 50, *)
          end if
        end do
      end do
    end do
    call r3mv( bvec, v1, s)
    dt = dt + norm2( s)
  end do
  close( 50)
 
  Return
end subroutine eqpongrid
