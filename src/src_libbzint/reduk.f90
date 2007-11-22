!BOP
!
! !ROUTINE: reduk
!
! !INTERFACE: 
subroutine reduk(nsymt,divsh,weight)
  !     
  ! !DESCRIPTION:
  !
  ! This subroutine calculates a set of reduced k-points for the integration
  ! with the tetrahedron method. The k-points are reduced by the symmetry and 
  ! a weight is put to each irreducible k-point to tell how many points it 
  ! represents before the symmetry operation.

  !
  ! !USES:

  use kgen_internals
  !<sag>
  use control, only: tetraifc
  !</sag>

  implicit none
  !<sag>
  logical, allocatable :: done(:)
  !</sag>


  ! !INPUT PARAMETERS:

  integer(4), intent(in) :: nsymt     ! Number of symmetry operations
  integer(4), intent(in) :: divsh

  ! !OUTPUT PARAMETERS:

  integer(4), intent(out):: weight(*) ! weight of the 
  !                                           reducible k-point
  ! !LOCAL VARIABLES:

  integer(4), allocatable :: kpav(:)

  integer(4) :: i,j,i1,i2,i3
  integer(4) :: onek(3),nkp(3),kk(3),kpid,nkpid,ktp(3)
  integer(4) :: members,starr(48)

  integer(4), external :: jget
  integer(4), external :: idkp

  !EOP
  !BOC
  external jset

  allocate(kpav(div(1)*div(2)*div(3)/32+1))

  nirkp=0
  do i=1,div(1)*div(2)*div(3)
     call jset(kpav,i,0)
  enddo
  !      write(6,*)'------------------------------------------------------'
  !      write(6,*)'               reduk: begin'
  !      write(6,*)'------------------------------------------------------'

  !<sag>
  if (trim(tetraifc)=='wien2k') then

     ! original code
     do i1=0,div(1)-1
        kk(1)=divsh*i1+shift(1)
        onek(1)=i1
        do i2=0,div(2)-1
           kk(2)=divsh*i2+shift(2)
           onek(2)=i2
           do i3=0,div(3)-1
              kk(3)=divsh*i3+shift(3)
              onek(3)=i3
              kpid=idkp(onek)
              if(jget(kpav,kpid).eq.0) then
                 nirkp=nirkp+1
                 ikpid(kpid)=nirkp
                 do i=1,48
                    starr(i)=0
                 enddo
                 !              write(24,*)
                 !              write(24,'(i6,3i4)')divsh,div
                 !              write(24,'(3i4,"  ",3i4)')onek,kk
                 do i=1,nsymt
                    do j=1,3
                       nkp(j)=mod(iio(j,1,i)*kk(1)+iio(j,2,i)*kk(2)+&
                            &                 iio(j,3,i)*kk(3),divsh*div(j))
                       ktp(j)=nkp(j)
                       nkp(j)=nkp(j)+(1-isign(1,nkp(j)))*divsh*div(j)/2
                       nkp(j)=(nkp(j)-shift(j))/divsh
                    enddo
                    !               write(24,'(i6,3i4,"  ",3i4)')i,ktp,nkp
                    nkpid=idkp(nkp)
                    !                write(24,'(4x,2i4)')nirkp,nkpid
                    call jset(kpav,nkpid,1)
                    redkp(nkpid)=kpid
                    starr(i)=nkpid
                 enddo
                 members=0
                 do i=1,nsymt
                    if(starr(i).ne.0) then
                       members=members+1
                       do j=i+1,nsymt
                          if(starr(i).eq.starr(j)) starr(j)=0
                       enddo
                    endif
                 enddo
                 weight(nirkp)=members
              else
                 ikpid(kpid)=0
              endif
              ! find coords of reduced onek
              call coorskp(redkp(kpid),nkp)
           enddo
        enddo
     enddo

     ! end original code
  else if (trim(tetraifc)=='exciting') then
     ! new code

     allocate(done(div(1)*div(2)*div(3)))
     done(:)=.false.


     do i3=0,div(3)-1
        kk(3)=divsh*i3+shift(3)
        onek(3)=i3
        do i2=0,div(2)-1
           kk(2)=divsh*i2+shift(2)
           onek(2)=i2
           do i1=0,div(1)-1
              kk(1)=divsh*i1+shift(1)
              onek(1)=i1
              kpid=idkp(onek)
write(*,'(a,3i6,i9)') 'libbzint: knr,kpid',i1,i2,i3,kpid
              !---            if(jget(kpav,kpid).eq.0) then
              if (.not.done(kpid)) then !+++
                 nirkp=nirkp+1
                 ikpid(kpid)=nirkp
                 do i=1,48
                    starr(i)=0
                 enddo
                 !              write(24,*)
                 !              write(24,'(i6,3i4)')divsh,div
                 !              write(24,'(3i4,"  ",3i4)')onek,kk
                 do i=1,nsymt
                    do j=1,3
                       nkp(j)=mod(iio(j,1,i)*kk(1)+iio(j,2,i)*kk(2)+&
                            &                 iio(j,3,i)*kk(3),divsh*div(j))
                       ktp(j)=nkp(j)
                       nkp(j)=nkp(j)+(1-isign(1,nkp(j)))*divsh*div(j)/2
                       nkp(j)=(nkp(j)-shift(j))/divsh
                    enddo
                    !               write(24,'(i6,3i4,"  ",3i4)')i,ktp,nkp
                    nkpid=idkp(nkp)
                    !                write(24,'(4x,2i4)')nirkp,nkpid
                    !---                call jset(kpav,nkpid,1)
                    done(nkpid)=.true. !+++
                    redkp(nkpid)=kpid
                    starr(i)=nkpid
!write(*,'(a,i9,3i6,i9)') '  sym,nkp,nkpid',i,nkp,nkpid
                 enddo
                 members=0
                 do i=1,nsymt
                    if(starr(i).ne.0) then
                       members=members+1
                       do j=i+1,nsymt
                          if(starr(i).eq.starr(j)) starr(j)=0
                       enddo
                    endif
                 enddo
                 weight(nirkp)=members
              else
                 ikpid(kpid)=0
              endif
              ! find coords of reduced onek
              call coorskp(redkp(kpid),nkp)
           enddo
        enddo
     enddo

     ! end new code
  end if ! if (tetraifc)
  !</sag>

  deallocate(iio)
  !      write(6,*)'------------------------------------------------------'
  !      write(6,*)'               reduk: end'
  !      write(6,*)'------------------------------------------------------'
  return
end subroutine reduk
!EOC
